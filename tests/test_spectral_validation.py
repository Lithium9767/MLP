"""Validate independent-group statistics rather than an invented function score."""
from unittest.mock import patch

import numpy as np
import pytest

from models.spectral_regions import FrozenSpectralModel, SpectralConfig
from models.spectral_validation import (
    _balanced_weights, _cross_group_distances, _profiles, _residual_ratio,
    summarize_validation,
)


def fixture():
    def row(sid, group, start, positions, split):
        return {"internal_id": sid, "homology_cluster": group, "split": split,
                "raw_start": start, "raw_end": start + 9, "sequence_length": 80,
                "reference_positions": ";".join(map(str, positions)), "accepted_mapping_count": 10}
    discovery = [row("d1", "g1", 1, range(11, 21), "discovery"),
                 row("d1", "g1", 6, range(25, 35), "discovery"),
                 row("d2", "g2", 1, range(11, 21), "discovery"),
                 row("d2", "g2", 6, range(25, 35), "discovery")]
    validation = [row("v1", "h1", 1, range(11, 21), "validation"),
                  row("v1", "h1", 6, range(25, 35), "validation"),
                  row("v2", "h2", 1, range(11, 21), "validation"),
                  row("v2", "h2", 6, range(25, 35), "validation")]
    train = np.array([[-1, 0, 0], [1, 0, 0], [-2, 0, .1], [2, 0, .1]], float)
    values = np.array([[-1.1, 0, 0], [1.1, 0, 0], [-2.1, 0, .5], [2.1, 0, .5]])
    labels = np.array([0, 1, 0, 1])
    model = FrozenSpectralModel(np.zeros(3), np.eye(3)[:2], train[:, :2], labels,
                                SpectralConfig(n_clusters=2, pca_components=2,
                                               graph_neighbors=2, transfer_neighbors=1), {})
    return discovery, train, labels, validation, values, model.predict(values), model


def test_replication_is_frozen_and_diagnostics_are_not_accuracy_or_rejection():
    args = fixture()
    before = args[-1].training_scores.copy()
    with patch("models.spectral_regions.PCA", side_effect=AssertionError("refit")), \
            patch("models.spectral_regions.SpectralClustering", side_effect=AssertionError("refit")):
        summary, clusters, diagnostics = summarize_validation(*args, permutations=20)
    assert summary["macro_reference_profile_TV"] == pytest.approx(0)
    assert summary["validation_supported_clusters"] == 2
    assert summary["cohort_intersections"] == {"sequence_ids": 0, "sequence_hashes": 0, "homology_groups": 0}
    assert summary["measured_function_accuracy"] is None
    assert summary["model_refitted"] is False
    assert summary["validation_used"] is True
    assert summary["within_sequence_position_shuffle"]["p_value_upper_tail_plus_one"] >= 1 / 21
    assert len(clusters) == 2 and len(diagnostics) == 4
    assert diagnostics[-1]["residual_above_discovery_p95"] is True
    assert all(not row["diagnostic_flags_are_rejections"] for row in diagnostics)
    assert clusters[0]["validation_fraction_alpha1_11_21"] == pytest.approx(1)
    np.testing.assert_array_equal(before, args[-1].training_scores)


def test_group_then_sequence_then_window_weighting_and_replication_invariance():
    rows = [{"internal_id": "a", "homology_cluster": "A"},
            {"internal_id": "a", "homology_cluster": "A"},
            {"internal_id": "b", "homology_cluster": "A"},
            {"internal_id": "c", "homology_cluster": "B"}]
    np.testing.assert_allclose(_balanced_weights(rows, np.arange(4)), [.125, .125, .25, .5])
    reference = np.array([[1., 0], [1, 0], [1, 0], [0, 1]])
    actual = _profiles(rows, reference, np.zeros(4, int), 1)
    np.testing.assert_allclose(actual, [[.5, .5]])
    # A repeated sequence with identical profile in group A must not increase A's weight.
    repeated = rows + [{"internal_id": "a_copy", "homology_cluster": "A"}]
    result = _profiles(repeated, np.vstack([reference, [1, 0]]), np.zeros(5, int), 1)
    np.testing.assert_allclose(result, actual)


def test_training_distance_excludes_whole_group_and_escalates_past_32():
    scores = np.arange(81, dtype=float)[:, None]
    groups = np.array(["A"] * 80 + ["B"])
    distances = _cross_group_distances(scores, groups)
    np.testing.assert_allclose(distances[:80], 80 - np.arange(80))
    assert distances[-1] == 1
    with pytest.raises(ValueError, match="at least two"):
        _cross_group_distances(scores, np.array(["A"] * 81))
    # Zero distance is legitimate if the equal points belong to different groups.
    np.testing.assert_array_equal(_cross_group_distances(np.zeros((2, 1)), np.array(["A", "B"])), [0, 0])


def test_zero_energy_residual_is_finite_zero():
    model = fixture()[-1]
    np.testing.assert_array_equal(_residual_ratio(np.zeros((2, 3)), model), [0, 0])


def test_validation_shift_cannot_move_discovery_thresholds_or_reject_labels():
    args = list(fixture())
    before = summarize_validation(*args, permutations=0)[0]
    args[4] = args[4] + np.array([100., 0., 100.])
    args[5] = args[-1].predict(args[4])
    after, _, diagnostics = summarize_validation(*args, permutations=0)
    for key in ("distance_diagnostic", "pca_residual_diagnostic"):
        assert before[key]["threshold"] == after[key]["threshold"]
        assert after[key]["validation_fraction_above_threshold"] == 1
    assert [row["candidate_cluster"] for row in diagnostics] == args[5]["cluster"].tolist()
    assert all(row["distance_above_discovery_cross_group_p95"] for row in diagnostics)


def test_validation_row_order_cannot_change_seeded_null():
    args = list(fixture())
    original = summarize_validation(*args, permutations=25)[0]
    order = np.array([3, 0, 2, 1])
    args[3] = [args[3][i] for i in order]
    args[4] = args[4][order]
    args[5] = {name: values[order] for name, values in args[5].items()}
    reordered = summarize_validation(*args, permutations=25)[0]
    assert original["within_sequence_position_shuffle"] == reordered["within_sequence_position_shuffle"]
    assert original["macro_reference_profile_TV"] == reordered["macro_reference_profile_TV"]


def test_missing_validation_cluster_is_unsupported_not_a_zero_distance_match():
    args = list(fixture())
    args[5]["cluster"] = np.zeros(4, dtype=int)
    summary, clusters, _ = summarize_validation(*args, permutations=0)
    assert summary["validation_supported_clusters"] == 1
    assert clusters[1]["validation_supported"] is False
    assert clusters[1]["reference_profile_TV"] is None
    assert clusters[1]["validation_fraction_unmapped"] is None
    assert summary["within_sequence_position_shuffle"]["p_value_upper_tail_plus_one"] is None


def test_unmapped_fraction_uses_all_ten_residues():
    args = list(fixture())
    args[3][0]["reference_positions"] = ";".join(map(str, range(11, 19)))
    args[3][0]["accepted_mapping_count"] = 8
    _, clusters, _ = summarize_validation(*args, permutations=0)
    assert clusters[0]["validation_fraction_unmapped"] == pytest.approx(.1)
    assert clusters[0]["validation_fraction_alpha1_11_21"] == pytest.approx(.9)
    assert clusters[0]["reference_profile_TV"] == pytest.approx(.1)


@pytest.mark.parametrize("kind", ["groups", "ids", "hashes", "labels", "features", "split", "duplicates", "tentative"])
def test_reject_leakage_misalignment_and_unaccepted_positions(kind):
    args = list(fixture())
    if kind == "groups":
        for row in args[3]: row["homology_cluster"] = "g1"
    elif kind == "ids":
        for row in args[3]: row["internal_id"] = "d1" if row["internal_id"] == "v1" else "d2"
    elif kind == "hashes":
        for row in args[0] + args[3]: row["sequence_sha256"] = "same"
    elif kind == "labels":
        args[2] = 1 - args[2]
    elif kind == "features":
        args[1] = args[1][::-1]
    elif kind == "split":
        args[3][0]["split"] = "new"
    elif kind == "duplicates":
        args[3][1] = dict(args[3][0])
    else:
        args[3][0]["reference_positions"] = ""
        args[3][0]["tentative_reference_positions"] = "11;12;13;14;15;16;17;18;19;20"
    with pytest.raises(ValueError):
        summarize_validation(*args, permutations=0)
