"""Synthetic checks for the M3 C handoff CLI; no discovery/validation run."""

import csv
from pathlib import Path
from types import SimpleNamespace

import pytest

from bioinformatics.m2_pf00741 import digest
from scripts.m3_coordinate_handoff import (
    DEFAULTS, PINNED_M2_RECEIPT_SHA256, ROOT, read_windows,
    require_validation_freeze,
)


def test_pinned_m2_receipts_match_committed_sources():
    for name, expected in PINNED_M2_RECEIPT_SHA256.items():
        assert digest(ROOT / DEFAULTS[name]) == expected


def _window_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("window_id", "internal_id",
                                                     "start_1based", "end_1based_inclusive"))
        writer.writeheader()
        writer.writerows(rows)


def _fake_index():
    return SimpleNamespace(records={
        "DISCOVERY_PRIMARY": {"sequence": "A" * 40, "split": "discovery",
                              "primary_analysis_eligible": True},
        "VALIDATION_PRIMARY": {"sequence": "A" * 40, "split": "validation",
                               "primary_analysis_eligible": True},
        "DISCOVERY_PARTIAL": {"sequence": "A" * 40, "split": "discovery",
                              "primary_analysis_eligible": False},
    })


def test_validation_requires_existing_nonempty_freeze_record(tmp_path):
    freeze = tmp_path / "candidate_freeze.md"
    with pytest.raises(ValueError, match="requires A's"):
        require_validation_freeze("validation", None)
    freeze.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="requires A's"):
        require_validation_freeze("validation", freeze)
    freeze.write_text("Frozen by A before validation.\n", encoding="utf-8")
    with pytest.raises(ValueError, match="front matter"):
        require_validation_freeze("validation", freeze)
    freeze.write_text(
        "---\nstatus: not_frozen\nfrozen_at_utc: 2026-09-29T10:00:00Z\n"
        f"candidate_manifest_sha256: {'a' * 64}\n"
        f"evaluation_rules_sha256: {'b' * 64}\n---\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="not frozen"):
        require_validation_freeze("validation", freeze)
    freeze.write_text(freeze.read_text(encoding="utf-8").replace("not_frozen", "frozen"),
                      encoding="utf-8")
    assert len(require_validation_freeze("validation", freeze)) == 64
    with pytest.raises(ValueError, match="only valid in validation"):
        require_validation_freeze("discovery", freeze)


def test_window_csv_requires_explicit_coordinates_and_frozen_queue(tmp_path):
    source = tmp_path / "windows.csv"
    _window_csv(source, [{"window_id": "w1", "internal_id": "DISCOVERY_PRIMARY",
                          "start_1based": "1", "end_1based_inclusive": "30"}])
    result = read_windows(source, _fake_index(), split="discovery", cohort="primary")
    assert result[0]["start_1based"] == 1
    assert result[0]["end_1based_inclusive"] == 30
    with pytest.raises(ValueError, match="belongs to discovery"):
        read_windows(source, _fake_index(), split="validation", cohort="primary")
    _window_csv(source, [{"window_id": "w2", "internal_id": "DISCOVERY_PARTIAL",
                          "start_1based": "1", "end_1based_inclusive": "30"}])
    with pytest.raises(ValueError, match="outside the primary queue"):
        read_windows(source, _fake_index(), split="discovery", cohort="primary")
    assert read_windows(source, _fake_index(), split="discovery", cohort="sensitivity")


@pytest.mark.parametrize("rows,error", [
    ([{"window_id": "w1", "internal_id": "UNKNOWN", "start_1based": "1",
       "end_1based_inclusive": "30"}], "Unknown window"),
    ([{"window_id": "w1", "internal_id": "DISCOVERY_PRIMARY", "start_1based": "0",
       "end_1based_inclusive": "30"}], "Out-of-range"),
    ([{"window_id": "w1", "internal_id": "DISCOVERY_PRIMARY", "start_1based": "1.0",
       "end_1based_inclusive": "30"}], "Noninteger"),
    ([{"window_id": "w1", "internal_id": "DISCOVERY_PRIMARY", "start_1based": "1",
       "end_1based_inclusive": "30"},
      {"window_id": "w2", "internal_id": "DISCOVERY_PRIMARY", "start_1based": "1",
       "end_1based_inclusive": "30"}], "Duplicate window coordinates"),
])
def test_window_csv_rejects_unknown_malformed_or_repeated_rows(tmp_path, rows, error):
    source = tmp_path / "windows.csv"
    _window_csv(source, rows)
    with pytest.raises(ValueError, match=error):
        read_windows(source, _fake_index(), split="discovery", cohort="primary")
