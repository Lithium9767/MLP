import pytest

from scripts.m3_discovery_coordinate_grid import discovery_windows, PUBLIC_RESIDUE_COLUMNS


def test_discovery_grid_isolated_and_no_padded_tail():
    records = {
        "primary": {"split": "discovery", "primary_analysis_eligible": True,
                    "sequence": "A" * 37},
        "short": {"split": "discovery", "primary_analysis_eligible": True,
                  "sequence": "A" * 29},
        "sensitivity": {"split": "discovery", "primary_analysis_eligible": False,
                        "sequence": "A" * 40},
        "validation": {"split": "validation", "primary_analysis_eligible": True,
                       "sequence": "A" * 40},
    }
    rows = list(discovery_windows(records, 30, 5))
    assert [(row["internal_id"], row["start_1based"], row["end_1based_inclusive"])
            for row in rows] == [("primary", 1, 30), ("primary", 6, 35)]
    assert len({row["window_id"] for row in rows}) == 2


@pytest.mark.parametrize("width,step", [(0, 5), (30, 0), (-1, 5)])
def test_invalid_grid_parameters(width, step):
    with pytest.raises(ValueError, match="positive"):
        list(discovery_windows({}, width, step))


def test_public_residue_table_does_not_redistribute_sequence_letters():
    assert "residue" not in PUBLIC_RESIDUE_COLUMNS
    assert "7r1c_reference_residue" not in PUBLIC_RESIDUE_COLUMNS
    assert {"raw_position", "hmm_match_state", "7r1c_modeled"} <= set(PUBLIC_RESIDUE_COLUMNS)
