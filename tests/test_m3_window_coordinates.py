"""M3 C coordinate handoff: synthetic edge cases, never discovery evidence."""

import csv
from copy import deepcopy

import pytest

from bioinformatics.m2_pf00741 import sequence_digest
from bioinformatics.m3_window_coordinates import CoordinateIndex, load_coordinate_index


SEQUENCE = "ACDEFGHIKL"


def record(identifier="GVPA_000001", *, split="discovery", primary=True):
    return {"internal_id": identifier, "sequence": SEQUENCE,
            "sequence_sha256": sequence_digest(SEQUENCE), "split": split,
            "analysis_cohort": "primary" if primary else "sensitivity_partial",
            "primary_analysis_eligible": str(primary).lower()}


def status(source, *, kind="accepted", domains=1, mapped=5, states=4):
    return {"internal_id": source["internal_id"], "sequence_sha256": source["sequence_sha256"],
            "split": source["split"], "analysis_cohort": source["analysis_cohort"],
            "is_partial_or_conflict": str(source["primary_analysis_eligible"] == "false"),
            "status": kind, "n_ga_domains": str(domains),
            "n_mapped_residues": str(mapped), "n_match_states": str(states)}


def coordinate(source, position, state, *, domain=1):
    return {"internal_id": source["internal_id"], "sequence_sha256": source["sequence_sha256"],
            "split": source["split"], "analysis_cohort": source["analysis_cohort"],
            "domain_number": str(domain), "raw_position": str(position),
            "residue": SEQUENCE[position - 1],
            "hmm_match_state": "" if state is None else str(state),
            "is_insertion": str(state is None)}


def structure(*, unmodeled_state=4):
    rows = []
    for state in range(1, 9):
        modeled = state != unmodeled_state
        rows.append({"deposited_position": str(state),
                     "pdb_residue_number": str(100 + state) if modeled else "",
                     "pdb_insertion_code": "", "residue": "A",
                     "modeled": str(modeled),
                     "secondary_structure": "H" if modeled else "",
                     "hmm_match_state": str(state),
                     "is_hmm_insertion": "False", "hmm_aligned": "True"})
    return rows


def fixture_rows():
    source = record()
    # HMM state 3 is deleted between natural residues 3 and 5; residue 4 is
    # a true alignment insertion and has no HMM or fabricated PDB coordinate.
    coords = [coordinate(source, pos, state) for pos, state in
              ((2, 1), (3, 2), (4, None), (5, 4), (6, 5))]
    return [source], [status(source)], coords, structure()


def test_window_maps_insertion_internal_deletion_and_unmodeled_reference():
    index = CoordinateIndex.from_rows(*fixture_rows(), hmm_length=8)
    mapped = index.map_window("GVPA_000001", 2, 6)
    assert mapped["start_1based"] == 2 and mapped["end_1based_inclusive"] == 6
    assert mapped["window_length"] == 5
    assert mapped["n_match_residues"] == 4
    assert mapped["n_insertion_residues"] == 1
    assert mapped["match_fraction"] == pytest.approx(0.8)
    assert mapped["coverage_status"] == "mapped"
    assert mapped["hmm_match_states"] == [1, 2, 4, 5]
    assert mapped["hmm_deleted_states_internal"] == [3]
    insertion = mapped["residues"][2]
    assert insertion["raw_position"] == 4
    assert insertion["mapping_status"] == "insertion"
    assert insertion["hmm_match_state"] is None
    assert insertion["reference_7r1c"] is None
    assert insertion["reference_status"] == "no_match_state"
    unmodeled = mapped["residues"][3]
    assert unmodeled["hmm_match_state"] == 4
    assert unmodeled["reference_status"] == "reference_unmodeled"
    assert unmodeled["reference_7r1c"]["deposited_position"] == 4
    assert unmodeled["reference_7r1c"]["pdb_residue_number"] is None


def test_unaligned_low_coverage_and_partial_are_not_force_mapped():
    source, statuses, coords, reference = fixture_rows()
    sensitivity = record("GVPA_000002", primary=False)
    source.append(sensitivity)
    statuses.append(status(sensitivity, kind="no_hit", domains=0, mapped=0, states=0))
    index = CoordinateIndex.from_rows(source, statuses, coords, reference, hmm_length=8)
    left = index.map_window("GVPA_000001", 1, 4)
    assert left["coverage_status"] == "low_match_coverage"
    assert left["n_unaligned_residues"] == 1
    assert left["residues"][0]["mapping_status"] == "unaligned"
    assert left["residues"][0]["reference_7r1c"] is None
    assert left["hmm_deleted_states_internal"] == []
    no_hit = index.map_window(sensitivity["internal_id"], 1, 5)
    assert no_hit["coverage_status"] == "unmapped"
    assert no_hit["primary_analysis_eligible"] is False
    assert no_hit["analysis_cohort"] == "sensitivity_partial"
    assert no_hit["hmm_match_states"] == []
    assert all(row["reference_7r1c"] is None for row in no_hit["residues"])


def test_failed_m2_alignment_retains_no_coordinates_even_with_reported_ga_domain():
    source = record("GVPA_000003")
    index = CoordinateIndex.from_rows([source], [status(source, kind="failed", domains=1,
                                                         mapped=0, states=0)], [], structure(),
                                      hmm_length=8)
    result = index.map_window(source["internal_id"], 1, 3)
    assert result["m2_hmm_status"] == "failed"
    assert result["coverage_status"] == "unmapped"
    assert result["hmm_match_states"] == []


def test_multi_domain_overlap_is_explicitly_ambiguous():
    source = record()
    coords = [coordinate(source, 2, 1, domain=1), coordinate(source, 3, 2, domain=1),
              coordinate(source, 3, 5, domain=2), coordinate(source, 4, 6, domain=2)]
    index = CoordinateIndex.from_rows([source], [status(source, kind="multi_hit", domains=2,
                                                         mapped=3, states=4)], coords,
                                      structure(), hmm_length=8)
    result = index.map_window(source["internal_id"], 2, 4)
    assert result["coverage_status"] == "ambiguous_multi_domain"
    assert result["n_ambiguous_residues"] == 1
    assert result["residues"][1]["hmm_match_state"] is None
    assert result["residues"][1]["reference_status"] == "ambiguous_multi_domain"
    assert len(result["residues"][1]["coordinate_options"]) == 2
    assert result["hmm_match_states"] == [1, 6]


def test_disjoint_domains_in_one_window_are_not_presented_as_one_mapping():
    source = record()
    coords = [coordinate(source, 2, 1, domain=1), coordinate(source, 3, 2, domain=1),
              coordinate(source, 5, 5, domain=2), coordinate(source, 6, 6, domain=2)]
    index = CoordinateIndex.from_rows([source], [status(source, kind="multi_hit", domains=2,
                                                         mapped=4, states=4)], coords,
                                      structure(), hmm_length=8)
    result = index.map_window(source["internal_id"], 2, 6)
    assert result["coverage_status"] == "multi_domain_window"
    assert result["n_domains"] == 2 and result["domain_numbers"] == [1, 2]
    assert result["n_unaligned_residues"] == 1
    assert result["hmm_deleted_states_internal"] == []


@pytest.mark.parametrize("change,reason", [
    (lambda records, statuses, coords, reference: coords[0].update(residue="W"), "residue differs"),
    (lambda records, statuses, coords, reference: coords.append(deepcopy(coords[0])), "Repeated domain/position"),
    (lambda records, statuses, coords, reference: coords[0].update(sequence_sha256="bad"), "identity differs"),
    (lambda records, statuses, coords, reference: statuses[0].update(n_mapped_residues="4"), "counts disagree"),
    (lambda records, statuses, coords, reference: statuses[0].update(internal_id="unknown"), "Unknown or duplicate status"),
    (lambda records, statuses, coords, reference: coords[2].update(is_insertion="False"), "Insertion must"),
    (lambda records, statuses, coords, reference: reference[3].update(pdb_residue_number="104"), "Unmodeled"),
    (lambda records, statuses, coords, reference: [row.update(domain_number="2") for row in coords], "domain numbering has gaps"),
])
def test_corrupted_m2_or_reference_rows_are_rejected(change, reason):
    rows = fixture_rows()
    change(*rows)
    with pytest.raises(ValueError, match=reason):
        CoordinateIndex.from_rows(*rows, hmm_length=8)


def test_unknown_or_out_of_bounds_window_rejected():
    index = CoordinateIndex.from_rows(*fixture_rows(), hmm_length=8)
    with pytest.raises(KeyError, match="Unknown internal_id"):
        index.map_window("not-in-input", 1, 2)
    for start, end in ((0, 1), (2, 11), (5, 4), (True, 2)):
        with pytest.raises(ValueError, match="Window"):
            index.map_window("GVPA_000001", start, end)
    with pytest.raises(ValueError, match="min_match_fraction"):
        index.map_window("GVPA_000001", 1, 2, min_match_fraction=1.5)


def test_csv_loader_uses_m2_tables_without_rerunning_hmm(tmp_path):
    records, statuses, coords, reference = fixture_rows()
    paths = [tmp_path / name for name in ("status.csv", "coordinates.csv", "7r1c.csv")]
    for path, rows in zip(paths, (statuses, coords, reference)):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    index = load_coordinate_index(records, *paths, hmm_length=8)
    assert index.map_window(records[0]["internal_id"], 2, 6)["hmm_deleted_states_internal"] == [3]
