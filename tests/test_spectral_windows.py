"""Input identity, frozen correspondence and exact feature arithmetic contracts."""

import csv
from unittest.mock import patch

import numpy as np
import pytest

from bioinformatics.m2_pf00741 import sequence_digest
from features.spectral_windows import (
    accepted_reference_mapping, build_window_features, load_hmm_anchors,
    load_primary_records,
)


REFERENCE = "MSIQKSTNSSSLAEVIDRILDKGIVIDAFARVSVVGIEILTIEARVVIASVDTWLRYAEAVGLLRDDVEENGLPERSNSSEGQPRFSI"


def record(sid="a", sequence=None, split="discovery", group="group_a"):
    if sequence is None:
        sequence = {"a": "M", "b": "A", "c": "C"}.get(sid, "G") + REFERENCE[1:]
    return dict(internal_id=sid, sequence_id="accession_" + sid, sequence=sequence,
                sequence_length=str(len(sequence)), sequence_sha256=sequence_digest(sequence),
                sequence_qc_eligible="true", primary_analysis_eligible="true",
                analysis_cohort="primary", split=split, homology_cluster=group)


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def source_files(tmp_path, records):
    fields = ("internal_id", "sequence_id", "analysis_cohort", "primary_analysis_eligible",
              "homology_cluster", "split")
    return (write_csv(tmp_path / "metadata.csv", records),
            write_csv(tmp_path / "split.csv", [{k: r[k] for k in fields} for r in records]))


def coordinate_rows(row):
    return [dict(internal_id=row["internal_id"], sequence_sha256=row["sequence_sha256"],
                 split=row["split"], analysis_cohort="primary", raw_position=str(h + 10),
                 residue=row["sequence"][h + 9], hmm_match_state=str(h)) for h in range(1, 40)]


def test_primary_records_join_sort_and_filter(tmp_path):
    a, b = record(), record("b", split="validation", group="other")
    nonprimary = record("c", group="third")
    nonprimary.update(primary_analysis_eligible="false", analysis_cohort="sensitivity")
    metadata, splits = source_files(tmp_path, [b, nonprimary, a])
    assert [r["internal_id"] for r in load_primary_records(metadata, splits, "discovery")] == ["a"]
    assert load_primary_records(metadata, splits, "validation")[0]["homology_cluster"] == "other"
    with pytest.raises(ValueError, match="Unknown split"):
        load_primary_records(metadata, splits, "test")


@pytest.mark.parametrize("field,value", [
    ("sequence_sha256", "bad"), ("sequence_length", "2"), ("sequence", "X" * 88),
    ("sequence_qc_eligible", "false"), ("analysis_cohort", "sensitivity"),
    ("sequence_id", ""), ("primary_analysis_eligible", "True"),
])
def test_bad_primary_identity_rejected(tmp_path, field, value):
    row = record()
    row[field] = value
    metadata, splits = source_files(tmp_path, [row])
    with pytest.raises(ValueError):
        load_primary_records(metadata, splits, "discovery")


def test_duplicate_metadata_and_cross_split_homology_rejected(tmp_path):
    a, b = record(), record("b", split="validation")
    metadata, splits = source_files(tmp_path, [a, a])
    with pytest.raises(ValueError, match="duplicate"):
        load_primary_records(metadata, splits, "discovery")
    metadata, splits = source_files(tmp_path, [a, b])
    with pytest.raises(ValueError, match="crosses"):
        load_primary_records(metadata, splits, "discovery")


@pytest.mark.parametrize("duplicate", ["sequence_id", "sequence"])
def test_duplicate_accession_or_sequence_under_different_ids_rejected(tmp_path, duplicate):
    a, b = record(), record("b", group="other")
    if duplicate == "sequence_id":
        b["sequence_id"] = a["sequence_id"]
    else:
        b["sequence"] = a["sequence"]
        b["sequence_sha256"] = a["sequence_sha256"]
    metadata, splits = source_files(tmp_path, [a, b])
    with pytest.raises(ValueError, match="duplicate"):
        load_primary_records(metadata, splits, "discovery")


def test_split_metadata_fields_and_coverage_are_checked(tmp_path):
    a, b = record(), record("b", group="other")
    metadata, splits = source_files(tmp_path, [a, b])
    with splits.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["sequence_id"] = "incorrect_accession"
    write_csv(splits, rows)
    with pytest.raises(ValueError, match="disagree"):
        load_primary_records(metadata, splits, "discovery")
    write_csv(splits, rows[1:])
    with pytest.raises(ValueError, match="cover"):
        load_primary_records(metadata, splits, "discovery")


def test_actual_hmm_anchors_preserved_without_imputation(tmp_path):
    row = record()
    coordinates = coordinate_rows(row)
    del coordinates[5]
    anchors = load_hmm_anchors(write_csv(tmp_path / "map.csv", coordinates), [row])["a"]
    assert len(anchors) == 38
    assert 6 not in anchors
    assert anchors[7] == 17


@pytest.mark.parametrize("field,value", [
    ("sequence_sha256", "bad"), ("split", "validation"), ("raw_position", "0"),
    ("residue", "X"), ("hmm_match_state", "40"), ("analysis_cohort", "sensitivity"),
])
def test_bad_hmm_coordinates_rejected(tmp_path, field, value):
    row = record()
    coordinates = coordinate_rows(row)
    coordinates[0][field] = value
    with pytest.raises(ValueError):
        load_hmm_anchors(write_csv(tmp_path / "map.csv", coordinates), [row])


def test_duplicate_hmm_states_positions_and_insufficient_support_rejected(tmp_path):
    row = record()
    rows = coordinate_rows(row)
    write_csv(tmp_path / "map.csv", rows + [rows[0]])
    with pytest.raises(ValueError, match="Duplicate"):
        load_hmm_anchors(tmp_path / "map.csv", [row])
    rows[1]["hmm_match_state"] = "1"
    with pytest.raises(ValueError, match="Duplicate HMM state"):
        load_hmm_anchors(write_csv(tmp_path / "map.csv", rows), [row])
    rows[1]["hmm_match_state"] = "2"
    with pytest.raises(ValueError, match="at least 20"):
        load_hmm_anchors(write_csv(tmp_path / "map.csv", rows[:10]), [row])


def test_identical_reference_extensions_and_missing_core_state():
    anchors = {h: h + 10 for h in range(1, 40) if h != 5}
    mapping, quality = accepted_reference_mapping(REFERENCE, anchors, REFERENCE)
    assert mapping == {q: q for q in range(2, 67) if q != 15}
    assert quality["N"]["passes_descriptive_sequence_support"]
    assert quality["alpha2"]["identity"] == 1
    assert quality["C_context"]["mapped_residues"] == 17
    assert 1 not in mapping and 67 not in mapping


def test_failed_extensions_do_not_fill_missing_core_states():
    anchors = {h: h + 10 for h in range(1, 40) if h != 5}
    empty = ({}, {"max_optimal_ties": 33, "truncated_ties": True, "scores": []})
    with patch("features.spectral_windows._extension_maps", return_value=empty):
        mapping, quality = accepted_reference_mapping(REFERENCE, anchors, REFERENCE)
    assert mapping == {q: h + 10 for h, q in anchors.items()}
    assert not quality["alpha2"]["passes_descriptive_sequence_support"]
    assert 15 not in mapping


def test_window_order_selection_and_exact_actual_anchor_background():
    row = record()
    anchors = {h: h + 10 for h in range(1, 40) if h != 5}
    mapping = {q: h + 10 for h, q in anchors.items()}
    # Different coefficients per channel expose accidental axis or whole-protein means.
    embedding = np.arange(88, dtype=np.float32)[:, None] * np.arange(1, 481, dtype=np.float32)
    selected, features, all_rows = build_window_features(
        [row], {"a": embedding}, {"a": anchors}, {"a": mapping})
    assert len(all_rows) == 16
    assert [r["raw_start"] for r in all_rows] == list(range(1, 77, 5))
    assert [r["raw_start"] for r in selected] == list(range(11, 42, 5))
    background = embedding[np.array(list(anchors.values())) - 1].mean(axis=0)
    expected = np.stack([embedding[r["raw_start"] - 1:r["raw_end"]].mean(axis=0) - background
                         for r in selected])
    np.testing.assert_array_equal(features, expected)
    assert features.dtype == np.float32
    assert selected[0]["reference_positions"] == "11;12;13;14;16;17;18;19;20"
    assert selected[0]["fragment"] == REFERENCE[10:20]
    assert not all_rows[0]["selected"]


@pytest.mark.parametrize("corruption", ["shape", "nonfinite", "dtype", "key", "hash", "mapping", "duplicate"])
def test_bad_feature_inputs_fail_before_prediction(corruption):
    row = record()
    anchors = {h: h + 10 for h in range(1, 40)}
    embedding = np.zeros((88, 480), np.float32)
    embeddings, mappings, records = {"a": embedding}, {"a": {q: q for q in range(11, 50)}}, [row]
    if corruption == "shape":
        embeddings["a"] = embedding[:-1]
    elif corruption == "nonfinite":
        embeddings["a"][0, 0] = np.nan
    elif corruption == "dtype":
        embeddings["a"] = embedding.astype(int)
    elif corruption == "key":
        embeddings = {}
    elif corruption == "hash":
        row["sequence_sha256"] = "bad"
    elif corruption == "mapping":
        mappings["a"][11] = 12
    elif corruption == "duplicate":
        records.append(row)
    with pytest.raises(ValueError):
        build_window_features(records, embeddings, {"a": anchors}, mappings)


def test_empty_selected_windows_have_explicit_shape():
    row = record()
    anchors = {h: h + 10 for h in range(1, 40)}
    selected, features, all_rows = build_window_features(
        [row], {"a": np.zeros((88, 480), np.float32)}, {"a": anchors},
        {"a": {q: q for q in range(11, 50)}}, width=88, min_accepted=88)
    assert selected == [] and features.shape == (0, 480)
    assert len(all_rows) == 1 and not all_rows[0]["selected"]
