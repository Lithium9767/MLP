"""Statistics for one frozen transfer to a reserved, group-disjoint cohort.

Reference profiles test whether a structure-informed interpretation transfers.
They are not function truth: reference mapping also selected these windows.
All thresholds below are discovery quantiles and only flag diagnostics; no
prediction is rejected, relabelled, or fitted by this module.
"""
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree


REGIONS = {
    "N_coil_2_10": range(2, 11),
    "alpha1_11_21": range(11, 22),
    "alpha1_beta_linker_22_23": range(22, 24),
    "beta1_24_34": range(24, 35),
    "beta_turn_35_36": range(35, 37),
    "beta2_37_48": range(37, 49),
    "beta_alpha2_linker_49": range(49, 50),
    "alpha2_50_62": range(50, 63),
    "C_coil_63_66": range(63, 67),
}


def _prepare_rows(rows, split):
    """Read accepted positions only; profiles retain all ten denominators."""
    if not rows:
        raise ValueError(f"{split} rows must not be empty")
    profiles = np.zeros((len(rows), 67), dtype=float)
    keys, sequence_info = set(), {}
    for i, row in enumerate(rows):
        if row["split"] != split:
            raise ValueError(f"Expected split={split}")
        sid, group = row["internal_id"], row["homology_cluster"]
        start, end, length = (int(row[key]) for key in
                              ("raw_start", "raw_end", "sequence_length"))
        if (not sid or not group or start < 1 or end > length
                or end - start + 1 != 10 or (start - 1) % 5):
            raise ValueError("Expected valid 10-aa windows with step 5")
        key = (sid, start, end)
        if key in keys:
            raise ValueError("Duplicate window key")
        keys.add(key)
        info = (group, length, row.get("sequence_sha256"))
        if sequence_info.setdefault(sid, info) != info:
            raise ValueError("Conflicting sequence metadata")
        text = row["reference_positions"]
        if not isinstance(text, str):
            raise ValueError("reference_positions must be a semicolon-separated string")
        try:
            positions = [int(value) for value in text.split(";") if value]
        except ValueError as error:
            raise ValueError("Invalid accepted reference position") from error
        if (not 8 <= len(positions) <= 10 or len(set(positions)) != len(positions)
                or any(p < 2 or p > 66 for p in positions)):
            raise ValueError("Expected 8-10 unique accepted positions in modeled reference 2-66")
        if "accepted_mapping_count" in row and int(row["accepted_mapping_count"]) != len(positions):
            raise ValueError("Accepted mapping count and positions disagree")
        profiles[i, np.array(positions) - 1] = 0.1
        profiles[i, 66] = (10 - len(positions)) / 10
    return profiles, sequence_info


def _balanced_weights(rows, indices):
    """Give groups, their sequences, and each sequence's windows equal weight."""
    indices = np.asarray(indices, dtype=int)
    if not len(indices):
        return np.empty(0)
    counts = Counter(rows[i]["internal_id"] for i in indices)
    group_by_sequence = {rows[i]["internal_id"]: rows[i]["homology_cluster"] for i in indices}
    sequences_by_group = Counter(group_by_sequence.values())
    n_groups = len(sequences_by_group)
    return np.array([1 / (n_groups * sequences_by_group[rows[i]["homology_cluster"]]
                         * counts[rows[i]["internal_id"]]) for i in indices])


def _profiles(rows, reference, labels, k):
    output = np.full((k, reference.shape[1]), np.nan)
    for cluster in range(k):
        indices = np.flatnonzero(labels == cluster)
        if len(indices):
            # Keep floating summation and the seeded permutation statistic stable
            # when an otherwise identical input table arrives in another order.
            indices = np.array(sorted(indices, key=lambda i: (
                rows[i]["homology_cluster"], rows[i]["internal_id"],
                int(rows[i].get("raw_start", 0)), int(rows[i].get("raw_end", 0)))))
            output[cluster] = _balanced_weights(rows, indices) @ reference[indices]
    return output


def _cross_group_distances(scores, groups):
    """Find nearest discovery points outside each query's entire homology group.

    Query more neighbours only for unresolved rows, avoiding an N-by-N matrix.
    Identical scores across different groups may legitimately have distance 0.
    """
    if len(set(groups)) < 2:
        raise ValueError("Cross-group distances require at least two discovery groups")
    tree = cKDTree(scores)
    output = np.empty(len(scores), dtype=float)
    unresolved = np.arange(len(scores))
    count = min(32, len(scores))
    while len(unresolved):
        remaining = []
        for begin in range(0, len(unresolved), 512):
            query = unresolved[begin:begin + 512]
            distances, indices = tree.query(scores[query], k=count, workers=1)
            distances = distances.reshape(len(query), count)
            indices = indices.reshape(len(query), count)
            outside = groups[indices] != groups[query, None]
            found = outside.any(axis=1)
            first = outside.argmax(axis=1)
            output[query[found]] = distances[np.flatnonzero(found), first[found]]
            remaining.extend(query[~found].tolist())
        unresolved = np.array(remaining, dtype=int)
        if len(unresolved) and count == len(scores):
            raise ValueError("No cross-group neighbour found")
        count = min(count * 2, len(scores))
    return output


def _residual_ratio(features, model):
    centered = np.asarray(features, dtype=float) - model.pca_mean
    projected = centered @ model.components.T
    residual = centered - projected @ model.components
    energy = np.einsum("ij,ij->i", centered, centered)
    residual_energy = np.einsum("ij,ij->i", residual, residual)
    return np.divide(residual_energy, energy, out=np.zeros(len(features)), where=energy > 0)


def _quantiles(values):
    levels = (0, 0.25, 0.5, 0.75, 0.95, 1)
    return {f"q{int(q * 100):02d}": float(value) for q, value in
            zip(levels, np.quantile(values, levels, method="linear"))}


def _support(rows, labels, cluster):
    indices = np.flatnonzero(labels == cluster)
    lengths = {rows[i]["internal_id"]: int(rows[i]["sequence_length"]) for i in indices}
    return {
        "windows": len(indices), "sequences": len(lengths),
        "groups": len({rows[i]["homology_cluster"] for i in indices}),
        "short_sequences_le100": sum(length <= 100 for length in lengths.values()),
        "mid_sequences_101_130": sum(101 <= length <= 130 for length in lengths.values()),
        "long_sequences_ge131": sum(length >= 131 for length in lengths.values()),
    }


def _region_fractions(profile):
    return {**{name: float(profile[np.array(list(positions)) - 1].sum())
               for name, positions in REGIONS.items()}, "unmapped": float(profile[66])}


def _shuffle_test(rows, reference, labels, discovery_profiles, observed, seed, permutations):
    """A conditional position-association diagnostic, never function evidence.

    Permuting labels only within a sequence preserves its cluster counts and
    group/sequence participation. Canonical raw-window order makes the seeded
    test invariant to the input table's row ordering.
    """
    by_sequence = {}
    for i, row in enumerate(rows):
        by_sequence.setdefault(row["internal_id"], []).append(i)
    blocks = [np.array(sorted(by_sequence[sid], key=lambda i:
                             (int(rows[i]["raw_start"]), int(rows[i]["raw_end"]))))
              for sid in sorted(by_sequence)]
    rng = np.random.default_rng(seed)
    null_scores = []
    for _ in range(permutations):
        shuffled = labels.copy()
        for indices in blocks:
            shuffled[indices] = rng.permutation(labels[indices])
        profile = _profiles(rows, reference, shuffled, len(discovery_profiles))
        supported = np.isfinite(profile).all(axis=1)
        tv = abs(profile[supported] - discovery_profiles[supported]).sum(axis=1) / 2
        null_scores.append(float(np.mean(1 - tv)))
    return {
        "seed": seed, "permutations": permutations,
        "statistic": "mean_1_minus_profile_TV_over_validation_supported_clusters",
        "observed": observed, "null_scores": null_scores,
        "p_value_upper_tail_plus_one": ((1 + sum(value >= observed for value in null_scores))
                                        / (permutations + 1)) if permutations else None,
        "interpretation": "Conditional within-sequence position association, not function accuracy.",
        "preserved": "Each sequence's cluster counts, homology groups and cluster participation.",
    }


def summarize_validation(discovery_rows, discovery_features, discovery_labels,
                         validation_rows, validation_features, prediction_dict,
                         frozen_model, *, seed=42, permutations=100):
    """Return summary, per-cluster comparisons and validation-window diagnostics.

    Rows are selected 10-aa windows with accepted ``reference_positions`` in a
    semicolon-separated string. All upstream cohort/feature provenance checks
    remain the caller's responsibility. This function additionally binds the
    discovery row order to the frozen model and rejects group/identity leakage.
    Excluded windows and zero-eligible sequences are not inputs to this API, so
    the caller must separately report those denominators for all 519 proteins.
    """
    if (isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0
            or isinstance(permutations, bool) or not isinstance(permutations, (int, np.integer))
            or permutations < 0):
        raise ValueError("seed and permutations must be nonnegative integers")
    discovery_reference, discovery_info = _prepare_rows(discovery_rows, "discovery")
    validation_reference, validation_info = _prepare_rows(validation_rows, "validation")
    discovery_groups = {info[0] for info in discovery_info.values()}
    validation_groups = {info[0] for info in validation_info.values()}
    discovery_hashes = {info[2] for info in discovery_info.values() if info[2]}
    validation_hashes = {info[2] for info in validation_info.values() if info[2]}
    overlaps = {"sequence_ids": len(set(discovery_info) & set(validation_info)),
                "sequence_hashes": len(discovery_hashes & validation_hashes),
                "homology_groups": len(discovery_groups & validation_groups)}
    if any(overlaps.values()):
        raise ValueError(f"Discovery/validation leakage: {overlaps}")
    k = frozen_model.config.n_clusters
    features = []
    for name, values, rows in [("discovery", discovery_features, discovery_rows),
                                ("validation", validation_features, validation_rows)]:
        values = np.asarray(values, dtype=float)
        if (values.shape != (len(rows), len(frozen_model.pca_mean))
                or not np.isfinite(values).all()):
            raise ValueError(f"Invalid {name} feature matrix")
        features.append(values)
    discovery_features, validation_features = features
    discovery_labels = np.asarray(discovery_labels)
    if not np.array_equal(discovery_labels, frozen_model.labels):
        raise ValueError("Discovery labels differ from the frozen model numbering/order")
    projected = (discovery_features - frozen_model.pca_mean) @ frozen_model.components.T
    if (projected.shape != frozen_model.training_scores.shape
            or not np.allclose(projected, frozen_model.training_scores, atol=1e-9, rtol=1e-9)):
        raise ValueError("Discovery features/order do not match the frozen PCA training scores")
    prediction = {name: np.asarray(prediction_dict[name]) for name in
                  ("cluster", "vote_fraction", "vote_margin", "nearest_distance")}
    if any(array.shape != (len(validation_rows),) or not np.isfinite(array).all()
           for array in prediction.values()):
        raise ValueError("Invalid prediction vectors")
    labels = prediction["cluster"]
    if (labels.dtype.kind not in "iu" or np.any(labels < 0) or np.any(labels >= k)
            or np.any(prediction["nearest_distance"] < 0)
            or any(np.any((prediction[name] < 0) | (prediction[name] > 1))
                   for name in ("vote_fraction", "vote_margin"))):
        raise ValueError("Prediction labels/diagnostics outside their valid ranges")
    discovery_profiles = _profiles(discovery_rows, discovery_reference, discovery_labels, k)
    validation_profiles = _profiles(validation_rows, validation_reference, labels, k)
    supported = np.isfinite(validation_profiles).all(axis=1)
    tv = abs(validation_profiles - discovery_profiles).sum(axis=1) / 2
    clusters = []
    for cluster in range(k):
        record = {"cluster": cluster, "reference_profile_TV": float(tv[cluster]) if supported[cluster] else None,
                  "validation_supported": bool(supported[cluster])}
        for prefix, rows, assigned, profile in [
                ("discovery", discovery_rows, discovery_labels, discovery_profiles[cluster]),
                ("validation", validation_rows, labels, validation_profiles[cluster])]:
            record.update({f"{prefix}_{name}": value
                           for name, value in _support(rows, assigned, cluster).items()})
            record.update({f"{prefix}_fraction_{name}": value if np.isfinite(value) else None
                           for name, value in _region_fractions(profile).items()})
            record[f"{prefix}_reference_profile_1_66_and_unmapped"] = (
                profile.tolist() if np.isfinite(profile).all() else None)
        clusters.append(record)
    training_groups = np.array([row["homology_cluster"] for row in discovery_rows])
    distances = _cross_group_distances(frozen_model.training_scores, training_groups)
    training_residual = _residual_ratio(discovery_features, frozen_model)
    validation_residual = _residual_ratio(validation_features, frozen_model)
    distance_threshold = float(np.quantile(distances, 0.95, method="linear"))
    residual_threshold = float(np.quantile(training_residual, 0.95, method="linear"))
    far = prediction["nearest_distance"] > distance_threshold
    high_residual = validation_residual > residual_threshold
    weights = _balanced_weights(validation_rows, np.arange(len(validation_rows)))
    diagnostics = []
    for i, row in enumerate(validation_rows):
        diagnostics.append({"internal_id": row["internal_id"], "homology_cluster": row["homology_cluster"],
                            "split": "validation", "raw_start": int(row["raw_start"]),
                            "raw_end": int(row["raw_end"]), "candidate_cluster": int(labels[i]),
                            **{name: float(prediction[name][i]) for name in
                               ("vote_fraction", "vote_margin", "nearest_distance")},
                            "pca_residual_energy_ratio": float(validation_residual[i]),
                            "distance_above_discovery_cross_group_p95": bool(far[i]),
                            "residual_above_discovery_p95": bool(high_residual[i]),
                            "diagnostic_flags_are_rejections": False})
    observed = float(np.mean(1 - tv[supported]))
    summary = {
        "operation": "frozen_reserved_cohort_transfer_statistics",
        "discovery_selected_windows": len(discovery_rows),
        "validation_selected_windows": len(validation_rows),
        "discovery_selected_sequences": len(discovery_info),
        "validation_selected_sequences": len(validation_info),
        "discovery_selected_groups": len(discovery_groups),
        "validation_selected_groups": len(validation_groups),
        "cohort_intersections": overlaps,
        "configured_clusters": k, "validation_supported_clusters": int(supported.sum()),
        "profile_weighting": "Within cluster: equal groups, equal sequences within group, equal windows within sequence.",
        "profile_denominator": "All 10 raw residues; reference positions 1-66 plus unmapped.",
        "macro_reference_profile_TV": float(np.mean(tv[supported])),
        "macro_reference_profile_overlap": observed,
        "distance_diagnostic": {
            "threshold": distance_threshold, "threshold_source": "Discovery nearest outside entire homology group; linear 95th percentile.",
            "discovery_cross_group_distance_quantiles": _quantiles(distances),
            "validation_nearest_distance_quantiles": _quantiles(prediction["nearest_distance"]),
            "validation_fraction_above_threshold": float(far.mean()),
            "validation_group_balanced_fraction_above_threshold": float(weights @ far),
        },
        "pca_residual_diagnostic": {
            "definition": "sum((centered - projected@components)^2) / sum(centered^2); zero for zero energy.",
            "threshold": residual_threshold, "threshold_source": "Discovery residual ratio; linear 95th percentile.",
            "discovery_quantiles": _quantiles(training_residual),
            "validation_quantiles": _quantiles(validation_residual),
            "validation_fraction_above_threshold": float(high_residual.mean()),
            "validation_group_balanced_fraction_above_threshold": float(weights @ high_residual),
        },
        "vote_fraction_quantiles": _quantiles(prediction["vote_fraction"]),
        "vote_margin_quantiles": _quantiles(prediction["vote_margin"]),
        "within_sequence_position_shuffle": _shuffle_test(
            validation_rows, validation_reference, labels, discovery_profiles, observed, int(seed), int(permutations)),
        "validation_used": True, "model_refitted": False,
        "diagnostic_flags_are_rejections": False, "measured_function_accuracy": None,
        "interpretation": "Reserved-group structural-profile replication under a structure-informed protocol; no biological function truth is available.",
        "excluded_window_denominators_included": False,
    }
    return summary, clusters, diagnostics
