"""Contracts for graph fitting, frozen label transfer, and safe persistence."""

import json
from unittest.mock import patch

import numpy as np
import pytest
from sklearn.datasets import make_blobs, make_moons
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits

from models.spectral_regions import (
    FrozenSpectralModel, SpectralConfig, fit_spectral, load_model, save_model,
)


def simple_model(neighbours=2):
    return FrozenSpectralModel(
        pca_mean=np.zeros(2), components=np.eye(2),
        training_scores=np.array([[-1, 0], [1, 0], [-2, 0], [2, 0]], float),
        labels=np.array([1, 0, 1, 0]),
        config=SpectralConfig(n_clusters=2, pca_components=2, graph_neighbors=2,
                              transfer_neighbors=neighbours),
        diagnostics={"purpose": "explicit vote fixture"},
    )


def test_spectral_resolves_nonconvex_moons():
    features, target = make_moons(n_samples=240, noise=0.02, random_state=3)
    with threadpool_limits(limits=1):
        model = fit_spectral(features, SpectralConfig(
            n_clusters=2, pca_components=2, graph_neighbors=10, transfer_neighbors=5))
    assert adjusted_rand_score(target, model.labels) > 0.98
    assert model.diagnostics["vote_fraction_is_probability"] is False


def test_explicit_k12_is_reproducible_and_roundtrips(tmp_path):
    features, target = make_blobs(n_samples=240, n_features=24, centers=12,
                                  cluster_std=0.15, random_state=7)
    config = SpectralConfig(graph_neighbors=21, transfer_neighbors=5)
    with threadpool_limits(limits=1):
        model = fit_spectral(features, config)
        repeated = fit_spectral(features, config)
    assert np.unique(model.labels).size == 12
    assert adjusted_rand_score(target, model.labels) == 1
    assert adjusted_rand_score(model.labels, repeated.labels) == 1
    path = tmp_path / "model.npz"
    save_model(model, path, "esm2-residue-mean-minus-hmm-anchor-mean/v1")
    restored = load_model(path, expected_schema=model.feature_schema)
    assert restored.feature_schema == model.feature_schema
    for field, actual in restored.predict(features[:15]).items():
        np.testing.assert_array_equal(actual, model.predict(features[:15])[field])
    with np.load(path, allow_pickle=False) as archive:
        assert all(archive[name].dtype.kind != "O" for name in archive.files)
        assert json.loads(str(archive["metadata"]))["config"]["n_clusters"] == 12


def test_transfer_does_not_fit_and_is_batch_invariant():
    model = simple_model()
    features = np.array([[-0.3, 0], [0.8, 0], [0, 0], [20, 0]], float)
    before = model.training_scores.copy()
    with patch("models.spectral_regions.PCA", side_effect=AssertionError("refit")), \
            patch("models.spectral_regions.SpectralClustering", side_effect=AssertionError("refit")), \
            patch("models.spectral_regions.kneighbors_graph", side_effect=AssertionError("refit")):
        batch = model.predict(features)
        singles = [model.predict(row[None, :]) for row in features]
    for key in batch:
        np.testing.assert_allclose(batch[key], [item[key][0] for item in singles])
    np.testing.assert_array_equal(before, model.training_scores)
    assert not model.training_scores.flags.writeable
    assert batch["nearest_distance"][-1] == 18
    # A faraway query still receives a label: no fabricated OOD rejection.
    assert batch["cluster"][-1] == 0


def test_zero_distance_and_deterministic_cluster_ties():
    model = simple_model()
    exact = model.predict(np.array([[-1, 0]], float))
    assert exact["cluster"][0] == 1
    assert exact["vote_fraction"][0] == 1
    assert exact["vote_margin"][0] == 1
    assert exact["nearest_distance"][0] == 0
    tied = model.predict(np.zeros((1, 2)))
    assert tied["cluster"][0] == 0
    assert tied["vote_fraction"][0] == 0.5
    assert tied["vote_margin"][0] == 0
    single = simple_model(neighbours=1).predict(np.array([[1, 0]], float))
    assert single["cluster"][0] == 0
    assert single["vote_fraction"][0] == 1


def test_zero_distance_neighbours_exclude_nonzero_votes():
    model = simple_model()
    model = FrozenSpectralModel(
        model.pca_mean, model.components,
        np.array([[0, 0], [0.001, 0], [-2, 0], [2, 0]], float),
        model.labels, model.config, {},
    )
    prediction = model.predict(np.zeros((1, 2)))
    assert prediction["cluster"][0] == 1
    assert prediction["vote_fraction"][0] == 1


def test_duplicate_zero_distance_neighbours_receive_equal_votes():
    model = simple_model(neighbours=3)
    model = FrozenSpectralModel(
        model.pca_mean, model.components,
        np.array([[0, 0], [0, 0], [-2, 0], [2, 0]], float),
        model.labels, model.config, {},
    )
    prediction = model.predict(np.zeros((1, 2)))
    assert prediction["cluster"][0] == 0
    assert prediction["vote_fraction"][0] == 0.5
    assert prediction["vote_margin"][0] == 0


@pytest.mark.parametrize("features", [
    np.array([[np.nan, 0]]), np.array([[np.inf, 0]]), np.ones((2, 3)),
    np.ones(2), np.array([["1", "2"]]), np.ones((2, 2), dtype=complex),
])
def test_transfer_rejects_invalid_features(features):
    with pytest.raises(ValueError):
        simple_model().predict(features)


def test_empty_batch_preserves_output_types():
    prediction = simple_model().predict(np.empty((0, 2)))
    assert all(value.shape == (0,) for value in prediction.values())
    assert prediction["cluster"].dtype == np.int64


@pytest.mark.parametrize("changes", [
    {"n_clusters": 1}, {"pca_components": 0}, {"graph_neighbors": 1},
    {"transfer_neighbors": 0}, {"seed": -1}, {"seed": 2**32},
    {"n_clusters": 12.0}, {"seed": True},
])
def test_config_rejects_invalid_values(changes):
    with pytest.raises(ValueError):
        SpectralConfig(**changes)


def test_fit_rejects_insufficient_dimensions_and_no_variance():
    for features in (np.ones((12, 24)), np.ones((24, 4)), np.ones((24, 24))):
        with pytest.raises(ValueError):
            fit_spectral(features)


def test_load_rejects_schema_mismatch_and_modified_numeric_data(tmp_path):
    path = tmp_path / "model.npz"
    save_model(simple_model(), path, "feature/v1")
    with pytest.raises(ValueError, match="expected_schema"):
        load_model(path, expected_schema="different/v1")
    with np.load(path, allow_pickle=False) as archive:
        changed = {name: archive[name].copy() for name in archive.files}
    changed["training_scores"][0, 0] += 1
    np.savez_compressed(path, **changed)
    with pytest.raises(ValueError, match="fingerprint"):
        load_model(path)


def test_model_rejects_invalid_label_contract():
    model = simple_model()
    for labels in (np.array([0, 0, 0, 0]), np.array([0, 1, 1, 2]),
                   np.array([0, 1]), np.array([0., 1., 1., 0.])):
        with pytest.raises(ValueError, match="labels"):
            FrozenSpectralModel(model.pca_mean, model.components,
                                model.training_scores, labels, model.config, {})


def test_model_rejects_nan_arrays_and_save_requires_schema(tmp_path):
    model = simple_model()
    mean = model.pca_mean.copy()
    mean[0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        FrozenSpectralModel(mean, model.components, model.training_scores,
                            model.labels, model.config, {})
    with pytest.raises(ValueError, match="feature_schema"):
        save_model(model, tmp_path / "invalid.npz", "")


def test_numpy_integer_config_roundtrips_as_plain_json(tmp_path):
    model = simple_model()
    model.config = SpectralConfig(n_clusters=np.int64(2), pca_components=np.int64(2),
                                   graph_neighbors=2, transfer_neighbors=2)
    path = tmp_path / "integer_config.npz"
    save_model(model, path, "feature/v1")
    assert load_model(path).config == model.config


def test_transfer_rejects_distance_overflow():
    with pytest.raises(ValueError, match="distances overflowed"):
        simple_model().predict(np.array([[1e308, 1e308]]))
