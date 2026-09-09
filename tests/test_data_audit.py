"""Handoff integrity and label-state checks using only synthetic sequences."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from bioinformatics.data_audit import audit_data


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "bioinformatics" / "examples"
FIELDS = ["sequence_id", "sequence", "sequence_length", "sequence_hash", "valid_residues",
          "accession", "database", "source_url", "retrieval_date", "organism", "taxonomy_id",
          "lineage", "label", "label_definition", "label_source", "evidence_level",
          "similarity_cluster", "split", "exclusion_reason"]


def row(identifier, sequence, **overrides):
    result = {field: "" for field in FIELDS}
    result.update(sequence_id=identifier, sequence=sequence, sequence_length=str(len(sequence)),
                  sequence_hash=hashlib.sha256(sequence.upper().encode("ascii")).hexdigest(),
                  valid_residues="true", label_definition="pending_prediction_target_confirmation",
                  evidence_level="unlabelled_sequence")
    result.update(overrides)
    return result


def write_csv(path, rows, fields=FIELDS):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def inputs(tmp_path, rows=None, fasta=None, split_rows=None):
    if rows is None:
        rows = [row("a", "ACD"), row("b", "ADE")]
    if fasta is None:
        fasta = "".join(f'>{item["sequence_id"]}\n{item["sequence"]}\n'
                        for item in rows if not item["exclusion_reason"])
    fasta_path, metadata_path = tmp_path / "sequences.fasta", tmp_path / "metadata.csv"
    fasta_path.write_text(fasta, encoding="utf-8")
    write_csv(metadata_path, rows)
    split_path = None
    if split_rows is not None:
        split_path = tmp_path / "split.csv"
        write_csv(split_path, split_rows, ["sequence_id", "similarity_cluster", "split"])
    return fasta_path, metadata_path, split_path


def codes(report, kind="errors"):
    return {item["code"] for item in report[kind]}


def test_unlabelled_is_valid_for_description_not_supervision(tmp_path):
    report = audit_data(*inputs(tmp_path))
    assert report["format_valid"]
    assert report["descriptive_sequence_analysis_ready"]
    assert report["labels"]["status"] == "unlabelled"
    assert report["labels"]["unknown_records"] == 2
    assert report["labels"]["value_counts"] == {}
    assert not report["labels"]["supervised_training_approved"]
    assert report["source_completeness"]["retrieval_date"]["missing_included"] == 2
    assert report["split"]["status"] == "not_provided"


def test_excluded_rows_are_not_required_in_fasta(tmp_path):
    rows = [row("a", "ACD"), row("excluded", "", valid_residues="false", exclusion_reason="empty_sequence")]
    report = audit_data(*inputs(tmp_path, rows))
    assert report["format_valid"]
    assert report["summary"]["excluded_records"] == 1
    assert report["summary"]["fasta_records"] == 1


def test_excluded_only_rows_cannot_appear_in_fasta(tmp_path):
    rows = [row("a", "ACD"), row("x", "ADE", exclusion_reason="other")]
    assert "fasta_extra_ids" in codes(audit_data(*inputs(tmp_path, rows, ">a\nACD\n>x\nADE\n")))


def test_duplicate_ids_not_silently_collapsed(tmp_path):
    report = audit_data(*inputs(tmp_path, [row("same", "ACD"), row("same", "ADE")]))
    assert {"duplicate_fasta_ids", "duplicate_metadata_ids"} <= codes(report)
    assert not report["format_valid"]


def test_fasta_and_metadata_must_match_ids_and_complete_sequence(tmp_path):
    report = audit_data(*inputs(tmp_path, fasta=">a\nACE\n>unexpected\nADE\n"))
    assert {"fasta_missing_ids", "fasta_extra_ids", "sequence_mismatch"} <= codes(report)


@pytest.mark.parametrize("change,expected", [
    ({"sequence_length": "4"}, "sequence_length_mismatch"),
    ({"sequence_length": "3.0"}, "sequence_length_mismatch"),
    ({"sequence_hash": "0" * 64}, "sequence_hash_mismatch"),
    ({"valid_residues": "false"}, "invalid_residue_flag"),
    ({"sequence": "A-C"}, "invalid_metadata_sequence"),
    ({"sequence": "ß"}, "invalid_metadata_sequence"),
])
def test_metadata_integrity_failures(tmp_path, change, expected):
    record = row("a", "ACD")
    record.update(change)
    report = audit_data(*inputs(tmp_path, [record], ">a\nACD\n"))
    assert expected in codes(report)


def test_ascii_case_normalization_is_consistent(tmp_path):
    report = audit_data(*inputs(tmp_path, [row("a", "acd")], ">a\nACD\n"))
    assert report["format_valid"]


@pytest.mark.parametrize("fasta,expected", [
    ("ACD", "fasta_parse_error"), (">\nACD", "fasta_parse_error"),
    (">a\nß", "invalid_fasta_sequence"), (">a\nAC*", "invalid_fasta_sequence"),
    (">a\nA C", "invalid_fasta_sequence"), ("", "empty_fasta"),
])
def test_bad_fasta_is_reported(tmp_path, fasta, expected):
    assert expected in codes(audit_data(*inputs(tmp_path, fasta=fasta)))


def test_partial_labels_and_zero_label_are_not_missing(tmp_path):
    rows = [row("a", "ACD", label="0", label_definition="synthetic_test_target",
                label_source="synthetic_test_only", evidence_level="synthetic"), row("b", "ADE")]
    report = audit_data(*inputs(tmp_path, rows))
    assert report["format_valid"]
    assert report["labels"]["status"] == "partially_labelled"
    assert report["labels"]["value_counts"] == {"0": 1}
    assert report["labels"]["unknown_records"] == 1


def test_missing_evidence_and_conflicting_labels_require_review(tmp_path):
    rows = [row("a", "ACD", label="0"), row("b", "ACD", label="1")]
    report = audit_data(*inputs(tmp_path, rows))
    assert report["format_valid"]  # Evidence is a readiness state, not file corruption.
    assert report["labels"]["status"] == "evidence_review_required"
    assert {"label_evidence_missing", "sequence_label_conflict", "retained_exact_duplicates"} <= codes(report, "warnings")
    assert report["labels"]["conflicts"][0]["labels"] == ["0", "1"]


def test_different_target_definitions_require_review(tmp_path):
    rows = [row("a", "ACD", label="0", label_definition="target_A", label_source="test", evidence_level="test"),
            row("b", "ADE", label="1", label_definition="target_B", label_source="test", evidence_level="test")]
    assert "label_definition_conflict" in codes(audit_data(*inputs(tmp_path, rows)), "warnings")


def test_valid_split_and_distributions(tmp_path):
    rows = [row("a", "ACD", organism="synthetic_one"), row("b", "ADEG", organism="synthetic_two")]
    split_rows = [dict(sequence_id="a", similarity_cluster="group_a", split="train"),
                  dict(sequence_id="b", similarity_cluster="group_b", split="test")]
    report = audit_data(*inputs(tmp_path, rows, split_rows=split_rows))
    assert report["format_valid"]
    assert report["split"]["status"] == "assignment_checks_passed"
    assert report["distributions"]["length_aa"] == {3: 1, 4: 1}
    assert report["distributions"]["organism"] == {"synthetic_one": 1, "synthetic_two": 1}
    assert report["distributions"]["label_unknown_records"] == 2


@pytest.mark.parametrize("split_rows,expected", [
    ([dict(sequence_id="a", similarity_cluster="g", split="train"),
      dict(sequence_id="a", similarity_cluster="h", split="test")], {"duplicate_split_ids", "split_missing_ids", "exact_sequence_cross_split"}),
    ([dict(sequence_id="a", similarity_cluster="g", split="train"),
      dict(sequence_id="b", similarity_cluster="g", split="test")], {"cluster_cross_split"}),
    ([dict(sequence_id="a", similarity_cluster="", split="validation"),
      dict(sequence_id="extra", similarity_cluster="x", split="test")], {"missing_similarity_cluster", "invalid_split_name", "split_missing_ids", "split_extra_ids"}),
])
def test_split_rejects_leakage_and_bad_coverage(tmp_path, split_rows, expected):
    report = audit_data(*inputs(tmp_path, split_rows=split_rows))
    assert expected <= codes(report)
    assert report["split"]["status"] == "invalid"


def test_exact_duplicate_sequences_cannot_cross_split_even_in_different_clusters(tmp_path):
    rows = [row("a", "ACD"), row("b", "ACD")]
    split_rows = [dict(sequence_id="a", similarity_cluster="g", split="train"),
                  dict(sequence_id="b", similarity_cluster="h", split="test")]
    assert "exact_sequence_cross_split" in codes(audit_data(*inputs(tmp_path, rows, split_rows=split_rows)))


def test_split_in_metadata_is_checked_and_external_conflicts_rejected(tmp_path):
    rows = [row("a", "ACD", similarity_cluster="g", split="train"),
            row("b", "ADE", similarity_cluster="g", split="test")]
    assert "cluster_cross_split" in codes(audit_data(*inputs(tmp_path, rows)))
    split_rows = [dict(sequence_id="a", similarity_cluster="g", split="train"),
                  dict(sequence_id="b", similarity_cluster="h", split="test")]
    assert "split_metadata_conflict" in codes(audit_data(*inputs(tmp_path, rows, split_rows=split_rows)))


@pytest.mark.parametrize("metadata", [
    "sequence_id,sequence_id\na,a\n",
    "sequence_id\na\n",
    "sequence_id,sequence,sequence_length,sequence_hash,exclusion_reason\na,ACD,3\n",
])
def test_malformed_csv_reported(tmp_path, metadata):
    fasta_path, metadata_path, _ = inputs(tmp_path)
    metadata_path.write_text(metadata, encoding="utf-8")
    assert "metadata_parse_error" in codes(audit_data(fasta_path, metadata_path))


def test_synthetic_fixture_cli_hashes_and_no_overwrite(tmp_path):
    out = tmp_path / "audit with spaces"
    command = [sys.executable, "-m", "bioinformatics.data_audit", "--fasta",
               str(EXAMPLES / "synthetic_original.fasta"), "--metadata",
               str(EXAMPLES / "synthetic_metadata.csv"), "--out-dir", str(out)]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    report = json.loads((out / "audit.json").read_text(encoding="utf-8"))
    assert report["summary"]["included_records"] == 4
    assert report["labels"]["status"] == "unlabelled"
    assert report["inputs"]["fasta"]["sha256"] == hashlib.sha256((EXAMPLES / "synthetic_original.fasta").read_bytes()).hexdigest()
    before = (out / "audit.json").read_bytes()
    assert subprocess.run(command, cwd=ROOT, capture_output=True).returncode == 2
    assert (out / "audit.json").read_bytes() == before


def test_cli_invalid_input_still_writes_report(tmp_path):
    fasta, metadata, _ = inputs(tmp_path, fasta=">a\nAC*\n")
    out = tmp_path / "invalid_audit"
    result = subprocess.run([sys.executable, "-m", "bioinformatics.data_audit", "--fasta", str(fasta),
                             "--metadata", str(metadata), "--out-dir", str(out)], cwd=ROOT, capture_output=True)
    assert result.returncode == 2
    report = json.loads((out / "audit.json").read_text(encoding="utf-8"))
    assert not report["format_valid"]
    assert "invalid_fasta_sequence" in codes(report)


def test_missing_file_returns_report(tmp_path):
    _, metadata, _ = inputs(tmp_path)
    report = audit_data(tmp_path / "missing.fa", metadata)
    assert {"input_unreadable", "fasta_parse_error"} <= codes(report)
