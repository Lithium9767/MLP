from copy import deepcopy

import pytest

from scripts.m3_candidate_coordinate_handoff import validate_window


def fixture():
    record = {"split": "discovery", "primary_analysis_eligible": True,
              "sequence_sha256": "abc", "sequence": "A" * 40}
    row = {"split": "discovery", "sequence_sha256": "abc", "raw_start": "1",
           "raw_end": "30", "window_length": "30", "sequence_length": "40",
           "hmm_states": "", "mapped_fraction": "0", "probability": "0", "cluster": "-1"}
    mapped = {"hmm_match_states_json": "[]", "match_fraction": 0}
    return row, record, mapped


def test_unmapped_noise_is_retained_without_placeholder_state():
    validate_window(*fixture())


@pytest.mark.parametrize("field,value", [
    ("split", "validation"), ("sequence_sha256", "wrong"),
    ("raw_end", "31"), ("hmm_states", "1"),
    ("mapped_fraction", "0.5"), ("probability", "nan"), ("cluster", "4")])
def test_reject_incompatible_d_window(field, value):
    row, record, mapped = fixture()
    row[field] = value
    with pytest.raises(ValueError):
        validate_window(row, record, mapped)


def test_reject_sensitivity_or_validation_input():
    for change in ({"primary_analysis_eligible": False}, {"split": "validation"}):
        row, record, mapped = fixture()
        record.update(change)
        with pytest.raises(ValueError, match="discovery primary"):
            validate_window(row, record, mapped)
