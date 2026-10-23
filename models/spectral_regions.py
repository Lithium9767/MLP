"""Spectral discovery and frozen nearest-neighbour transfer for window features.

Spectral clustering is transductive: it supplies labels for its training graph.
``predict`` is an explicit extension that transfers those labels in the frozen
PCA space. Vote fractions describe agreement among neighbours, not calibrated
probabilities, biological function accuracy, or detection of unfamiliar inputs.
"""

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from sklearn.cluster import SpectralClustering
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph


@dataclass(frozen=True)
class SpectralConfig:
    """All fit/transfer parameters; no search or automatic choice of K."""

    n_clusters: int = 12
    pca_components: int = 20
    graph_neighbors: int = 15
    transfer_neighbors: int = 15
    seed: int = 42

    def __post_init__(self):
        for name, value in asdict(self).items():
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise ValueError(f"{name} must be an integer")
            object.__setattr__(self, name, int(value))
        if self.n_clusters < 2 or self.pca_components < 1:
            raise ValueError("n_clusters >= 2 and pca_components >= 1 are required")
        if self.graph_neighbors < 2 or self.transfer_neighbors < 1:
            raise ValueError("graph_neighbors >= 2 and transfer_neighbors >= 1 are required")
        if not 0 <= self.seed <= 2**32 - 1:
            raise ValueError("seed must be in [0, 2**32 - 1]")


def _matrix(values, name, *, columns=None, allow_empty=False):
    array = np.asarray(values)
    if array.ndim != 2 or array.dtype.kind not in "fiu":
        raise ValueError(f"{name} must be a real numeric two-dimensional array")
    if array.shape[1] < 1 or (not allow_empty and array.shape[0] < 1):
        raise ValueError(f"{name} must have nonzero dimensions")
    if columns is not None and array.shape[1] != columns:
        raise ValueError(f"{name} feature dimension must equal {columns}")
    if not np.isfinite(array).all():
        raise ValueError(f"{name} must contain only finite values")
    with np.errstate(over="ignore"):
        result = np.array(array, dtype=np.float64, order="C", copy=True)
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must be finite in float64 precision")
    return result


@dataclass
class FrozenSpectralModel:
    """The fitted PCA and labelled training points, with no refitting in predict.

    Numeric arrays are copied and read-only. ``feature_schema`` identifies the
    meaning and order of input features and is checked on model loading.
    """

    pca_mean: np.ndarray
    components: np.ndarray
    training_scores: np.ndarray
    labels: np.ndarray
    config: SpectralConfig
    diagnostics: dict
    feature_schema: str = ""

    def __post_init__(self):
        if not isinstance(self.config, SpectralConfig):
            raise ValueError("config must be a SpectralConfig")
        mean = np.asarray(self.pca_mean)
        if mean.ndim != 1:
            raise ValueError("pca_mean must be one-dimensional")
        mean = _matrix(mean[None, :], "pca_mean")[0]
        components = _matrix(self.components, "components", columns=len(mean))
        if len(components) != self.config.pca_components:
            raise ValueError("components do not match pca_components")
        scores = _matrix(self.training_scores, "training_scores", columns=len(components))
        labels = np.asarray(self.labels)
        if labels.ndim != 1 or labels.dtype.kind not in "iu" or len(labels) != len(scores):
            raise ValueError("labels must be an integer vector matching training_scores")
        if not np.array_equal(np.unique(labels), np.arange(self.config.n_clusters)):
            raise ValueError("labels must contain exactly the configured cluster IDs")
        if len(scores) <= self.config.n_clusters:
            raise ValueError("training sample count must exceed n_clusters")
        if self.config.graph_neighbors > len(scores) or self.config.transfer_neighbors > len(scores):
            raise ValueError("neighbour counts must not exceed training sample count")
        if not isinstance(self.feature_schema, str) or not isinstance(self.diagnostics, dict):
            raise ValueError("feature_schema must be a string and diagnostics a dict")
        json.dumps(self.diagnostics, allow_nan=False)
        for name, value in (("pca_mean", mean), ("components", components),
                            ("training_scores", scores),
                            ("labels", np.array(labels, dtype=np.int64, copy=True))):
            value.setflags(write=False)
            setattr(self, name, value)

    def predict(self, features):
        """Transfer frozen labels with inverse-distance neighbour voting.

        Exact zero-distance neighbours vote alone, with equal weights. Otherwise
        weights are proportional to inverse distance. A tied cluster vote picks
        the lowest cluster ID. ``vote_margin`` is the winning share minus the
        runner-up share; a unanimous vote has margin 1. No rejection threshold or
        probability calibration is assumed. Empty batches are accepted.
        """
        values = _matrix(features, "features", columns=len(self.pca_mean), allow_empty=True)
        output = {
            "cluster": np.empty(len(values), dtype=np.int64),
            "vote_fraction": np.empty(len(values), dtype=np.float64),
            "vote_margin": np.empty(len(values), dtype=np.float64),
            "nearest_distance": np.empty(len(values), dtype=np.float64),
        }
        if not len(values):
            return output
        scores = (values - self.pca_mean) @ self.components.T
        if not np.isfinite(scores).all():
            raise ValueError("PCA projection produced nonfinite scores")
        tree = cKDTree(self.training_scores)
        distances, indices = tree.query(scores, k=self.config.transfer_neighbors, workers=1)
        distances = distances.reshape(len(scores), self.config.transfer_neighbors)
        indices = indices.reshape(len(scores), self.config.transfer_neighbors)
        if not np.isfinite(distances).all():
            raise ValueError("query distances overflowed float64 precision")
        for row, (distance, index) in enumerate(zip(distances, indices)):
            zero = distance == 0
            # Scaling by the nearest distance preserves inverse-distance ratios
            # and avoids overflowing when a positive distance is extremely small.
            weights = zero.astype(float) if zero.any() else distance[0] / distance
            votes = np.bincount(self.labels[index], weights=weights,
                                minlength=self.config.n_clusters)
            fractions = votes / votes.sum()
            winner = int(np.argmax(fractions))
            runner_up = np.partition(fractions, -2)[-2]
            output["cluster"][row] = winner
            output["vote_fraction"][row] = fractions[winner]
            output["vote_margin"][row] = fractions[winner] - runner_up
            output["nearest_distance"][row] = distance[0]
        return output


def fit_spectral(features, config=SpectralConfig()):
    """Fit full-SVD PCA, a symmetric neighbour graph, and spectral cluster_qr.

    Self-neighbours count toward ``graph_neighbors`` to reproduce the exploratory
    graph definition. The graph is (A + A.T)/2, with edge weights 0, 0.5, or 1.
    Structural/function annotations are not inputs to this function.
    """
    values = _matrix(features, "features")
    if not isinstance(config, SpectralConfig):
        raise ValueError("config must be a SpectralConfig")
    n_samples, n_features = values.shape
    if n_samples <= config.n_clusters:
        raise ValueError("training sample count must exceed n_clusters")
    if config.pca_components > min(n_samples - 1, n_features):
        raise ValueError("pca_components exceeds the available feature/sample dimensions")
    if config.graph_neighbors > n_samples or config.transfer_neighbors > n_samples:
        raise ValueError("neighbour counts must not exceed training sample count")
    if np.all(values == values[0]):
        raise ValueError("training features must have nonzero variance")
    pca = PCA(n_components=config.pca_components, svd_solver="full")
    scores = pca.fit_transform(values)
    graph = kneighbors_graph(scores, n_neighbors=config.graph_neighbors,
                             mode="connectivity", include_self=True, n_jobs=1)
    graph = ((graph + graph.T) * 0.5).tocsr()
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        labels = SpectralClustering(
            n_clusters=config.n_clusters, n_components=config.n_clusters,
            affinity="precomputed", assign_labels="cluster_qr",
            eigen_solver="arpack", random_state=config.seed,
        ).fit_predict(graph)
    diagnostics = {
        "n_training_windows": n_samples,
        "n_input_features": n_features,
        "graph_connected_components": int(connected_components(graph, directed=False)[0]),
        "pca_explained_variance_ratio": float(pca.explained_variance_ratio_.sum()),
        "warnings": [str(item.message) for item in captured],
        "transfer_method": "inverse_distance_neighbours_in_frozen_pca",
        "vote_fraction_is_probability": False,
    }
    return FrozenSpectralModel(pca.mean_, pca.components_, scores, labels, config, diagnostics)


def _fingerprint(arrays, metadata):
    digest = hashlib.sha256(json.dumps(metadata, sort_keys=True, ensure_ascii=False,
                                        allow_nan=False).encode("utf-8"))
    for name, array in sorted(arrays.items()):
        digest.update(name.encode("ascii"))
        digest.update(str(array.dtype).encode("ascii"))
        digest.update(str(array.shape).encode("ascii"))
        digest.update(np.ascontiguousarray(array).tobytes())
    return digest.hexdigest()


def save_model(model, path, feature_schema):
    """Save numeric arrays plus versioned JSON metadata; never serialize pickle.

    The caller supplies a nonempty schema describing the feature construction,
    not merely its dimension. The target's parent directory must already exist.
    """
    if not isinstance(feature_schema, str) or not feature_schema.strip():
        raise ValueError("feature_schema must be a nonempty string")
    arrays = {name: getattr(model, name) for name in
              ("pca_mean", "components", "training_scores", "labels")}
    # Revalidate before writing in case a caller reassigned a model attribute.
    FrozenSpectralModel(**arrays, config=model.config, diagnostics=model.diagnostics,
                        feature_schema=feature_schema)
    metadata = {"format_version": 1, "config": asdict(model.config),
                "feature_schema": feature_schema, "diagnostics": model.diagnostics}
    metadata["fingerprint_sha256"] = _fingerprint(arrays, metadata)
    with Path(path).open("wb") as handle:
        np.savez_compressed(handle, **arrays,
                            metadata=np.array(json.dumps(metadata, ensure_ascii=False,
                                                        sort_keys=True, allow_nan=False)))
    model.feature_schema = feature_schema


def load_model(path, expected_schema=None):
    """Load a validated model, rejecting tampering and optional schema mismatch."""
    if expected_schema is not None and not isinstance(expected_schema, str):
        raise ValueError("expected_schema must be a string or None")
    required = {"pca_mean", "components", "training_scores", "labels", "metadata"}
    with np.load(path, allow_pickle=False) as archive:
        if set(archive.files) != required:
            raise ValueError("model archive has unexpected or missing fields")
        scalar = archive["metadata"]
        if scalar.ndim != 0 or scalar.dtype.kind != "U":
            raise ValueError("model metadata must be a Unicode JSON scalar")
        metadata = json.loads(str(scalar))
        arrays = {name: archive[name].copy() for name in required - {"metadata"}}
    fields = {"format_version", "config", "feature_schema", "diagnostics", "fingerprint_sha256"}
    if not isinstance(metadata, dict) or set(metadata) != fields:
        raise ValueError("model metadata has unexpected or missing fields")
    if type(metadata["format_version"]) is not int or metadata["format_version"] != 1:
        raise ValueError("unsupported model format_version")
    fingerprint = metadata.pop("fingerprint_sha256", None)
    if fingerprint != _fingerprint(arrays, metadata):
        raise ValueError("model fingerprint does not match its contents")
    schema = metadata.get("feature_schema")
    if not isinstance(schema, str) or not schema.strip():
        raise ValueError("model feature_schema must be a nonempty string")
    if expected_schema is not None and schema != expected_schema:
        raise ValueError("model feature_schema does not match expected_schema")
    config_values = metadata["config"]
    if not isinstance(config_values, dict) or set(config_values) != set(asdict(SpectralConfig())):
        raise ValueError("model config has unexpected or missing fields")
    return FrozenSpectralModel(**arrays, config=SpectralConfig(**config_values),
                               diagnostics=metadata["diagnostics"], feature_schema=schema)
