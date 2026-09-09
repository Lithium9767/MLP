"""Pfam TSV integrity, inclusive coordinate and unknown-status regression checks."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from bioinformatics.msa_analysis import read_fasta
from bioinformatics.pfam import (
    read_interproscan_tsv, residue_annotations, sequence_coverage, validate_accession,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "bioinformatics" / "examples"


def make_row(sequence="ACDE", key="a", start="1", end="4", signature="PF00741"):
    return [key, hashlib.md5(sequence.encode("ascii")).hexdigest(), str(len(sequence)),
            "Pfam", signature, "SYNTHETIC mock hit; not a biological annotation",
            start, end, "1e-10", "T", "09-09-2026", "-", "-"]


def parse_rows(tmp_path, rows, originals=None):
    path = tmp_path / "annotations.tsv"
    path.write_text("\n".join("\t".join(row) for row in rows) + "\n", encoding="utf-8")
    return read_interproscan_tsv(path, {"a": "ACDE"} if originals is None else originals, "PF00741")


def test_fixture_keeps_overlaps_versioned_accession_and_uncertain_sequences():
    originals = read_fasta(EXAMPLES / "synthetic_original.fasta", aligned=False)
    hits, audit = read_interproscan_tsv(EXAMPLES / "synthetic_interproscan.tsv", originals, "PF00741")
    assert len(hits) == 2
    assert hits[1]["signature_accession"] == "PF00741.1"
    assert hits[1]["pfam_accession"] == "PF00741"
    assert audit["n_input_rows"] == 4
    assert audit["n_filtered_rows"] == 2
    assert audit["selected_hit_reported"] == ["toy_1"]
    assert audit["no_selected_hit_uncertain"] == ["toy_2", "toy_3"]
    assert audit["not_reported_uncertain"] == ["toy_4"]
    assert audit["all_sequences_reported_hit_fraction"] == 0.25
    rows = list(residue_annotations(originals, hits, audit, "PF00741"))
    toy1 = [row for row in rows if row["sequence_id"] == "toy_1"]
    assert [row["residue_position"] for row in toy1] == list(range(1, 18))
    assert "".join(row["residue"] for row in toy1) == originals["toy_1"]
    assert toy1[0]["in_selected_hit"] == "true"
    assert toy1[3]["n_covering_hits"] == 2
    assert toy1[3]["hit_ids"] == "tsv-line-1|tsv-line-2"
    assert toy1[15]["in_selected_hit"] == "true"
    assert toy1[16]["in_selected_hit"] == "false"
    assert toy1[16]["n_covering_hits"] == 0
    assert toy1[16]["residue_fraction"] == 1
    assert toy1[0]["residue_fraction"] == pytest.approx(1 / 17)
    for row in rows:
        if row["sequence_id"] != "toy_1":
            assert row["in_selected_hit"] == ""
            assert row["n_covering_hits"] is None
            assert row["annotation_status"].endswith("uncertain")
    coverage = list(sequence_coverage(originals, hits, audit, "PF00741"))
    assert coverage[0]["covered_residue_count"] == 16  # overlap counted once
    assert coverage[0]["reported_interval_coverage_fraction"] == pytest.approx(16 / 17)
    assert coverage[1]["covered_residue_count"] is None
    assert coverage[3]["reported_interval_coverage_fraction"] is None


def test_single_residue_hit_at_final_position_and_duplicates_are_retained(tmp_path):
    row = make_row(start="4", end="4")
    hits, audit = parse_rows(tmp_path, [row, row])
    assert len(hits) == 2
    annotations = list(residue_annotations({"a": "ACDE"}, hits, audit, "PF00741"))
    assert [row["in_selected_hit"] for row in annotations] == ["false", "false", "false", "true"]
    assert annotations[-1]["n_covering_hits"] == 2
    assert len(set(hit["hit_id"] for hit in hits)) == 2
    coverage = list(sequence_coverage({"a": "ACDE"}, hits, audit, "PF00741"))[0]
    assert coverage["covered_residue_count"] == 1
    assert coverage["reported_interval_coverage_fraction"] == 0.25


def test_empty_report_is_unknown_for_every_sequence(tmp_path):
    hits, audit = parse_rows(tmp_path, [])
    assert not hits
    assert audit["n_input_rows"] == 0
    assert audit["not_reported_uncertain"] == ["a"]
    assert audit["all_sequences_reported_hit_fraction"] == 0
    annotations = list(residue_annotations({"a": "ACDE"}, hits, audit, "PF00741"))
    assert all(row["in_selected_hit"] == "" for row in annotations)
    assert all(row["annotation_status"] == "not_reported_uncertain" for row in annotations)


@pytest.mark.parametrize("count", [11, 13, 14, 15])
def test_documented_column_counts_preserve_optional_text(tmp_path, count):
    row = make_row()[:11] if count == 11 else make_row()
    row += ["GO:0005515(InterPro)", "Reactome:R-SYNTHETIC"][:max(0, count - 13)]
    hits, _ = parse_rows(tmp_path, [row])
    assert json.loads(hits[0]["optional_columns_json"]) == row[13:]
    assert hits[0]["interpro_accession"] == "-"


@pytest.mark.parametrize("column,value,match", [
    (0, "missing", "unknown sequence_id"),
    (1, "0" * 32, "MD5 mismatch"), (1, "bad", "invalid sequence MD5"),
    (2, "3", "length mismatch"), (2, "4.0", "positive integer"), (2, "0", "positive integer"),
    (3, "-", "invalid analysis"), (4, "PF0074", "invalid Pfam"),
    (6, "0", "positive integer"), (6, "-1", "positive integer"),
    (6, "2.0", "positive integer"), (6, "5", "start <= end"),
    (7, "5", "sequence length"), (7, "0", "positive integer"),
    (8, "nan", "invalid match score"), (8, "Infinity", "invalid match score"),
    (8, "-0.01", "invalid match score"), (8, "1_0", "invalid match score"),
    (9, "F", "status T"), (10, "31-02-2026", "invalid run date"),
    (10, "2026-09-09", "DD-MM-YYYY"), (11, "IPRxyz", "invalid InterPro"),
    (12, "orphan description", "description without accession"),
    (5, "", "empty TSV field"), (0, " a", "whitespace"),
])
def test_rejects_malformed_fields(tmp_path, column, value, match):
    row = make_row()
    row[column] = value
    with pytest.raises(ValueError, match=match):
        parse_rows(tmp_path, [row])


@pytest.mark.parametrize("count", [1, 10, 12, 16])
def test_rejects_wrong_field_count(tmp_path, count):
    row = (make_row() + ["-"] * 3)[:count]
    with pytest.raises(ValueError, match="TSV fields"):
        parse_rows(tmp_path, [row])


def test_filtered_rows_still_require_identity_and_coordinate_validation(tmp_path):
    row = make_row(signature="PF00001")
    row[1] = "0" * 32
    with pytest.raises(ValueError, match="MD5 mismatch"):
        parse_rows(tmp_path, [row])
    row = make_row(signature="G3DSA:1.10.10.10", end="5")
    row[3] = "Gene3D"
    with pytest.raises(ValueError, match="sequence length"):
        parse_rows(tmp_path, [row])


def test_uppercase_md5_tiny_evalue_and_interpro_metadata_are_preserved(tmp_path):
    row = make_row()
    row[1] = row[1].upper()
    row[8] = "1e-1000"
    row[11:13] = ["IPR000001", "SYNTHETIC InterPro description"]
    hits, _ = parse_rows(tmp_path, [row])
    assert hits[0]["score"] == "1e-1000"
    assert hits[0]["sequence_md5"] == row[1].lower()
    assert hits[0]["interpro_accession"] == "IPR000001"


@pytest.mark.parametrize("accession", ["", "PF741", "PF00741.1", "pf00741", "PF00741,PF01132"])
def test_requires_explicit_canonical_pfam_accession(accession):
    with pytest.raises(ValueError, match="Pfam accession"):
        validate_accession(accession)


def command(out):
    return [sys.executable, "-m", "bioinformatics.pfam",
            "--original", str(EXAMPLES / "synthetic_original.fasta"),
            "--interproscan-tsv", str(EXAMPLES / "synthetic_interproscan.tsv"),
            "--pfam-accession", "PF00741", "--out-dir", str(out),
            "--run-id", "test-pfam-synthetic", "--dataset-version", "synthetic-v1",
            "--purpose", "synthetic", "--interproscan-version", "synthetic-not-run",
            "--pfam-version", "synthetic-not-scanned",
            "--source-description", "Hand-constructed synthetic fixture; not a biological result"]


def test_cli_hashes_audit_coverage_and_no_overwrite(tmp_path):
    out = tmp_path / "pfam output with spaces"
    result = subprocess.run(command(out), cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["n_sequences"] == 4
    assert manifest["n_selected_hits"] == 2
    assert manifest["purpose"] == "synthetic"
    assert manifest["interproscan_version"] == "synthetic-not-run"
    assert manifest["pfam_version"] == "synthetic-not-scanned"
    assert manifest["join_keys"] == ["sequence_id", "residue_position"]
    for name, digest in manifest["outputs"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == digest
    assert manifest["interproscan_tsv"]["sha256"] == hashlib.sha256(
        (EXAMPLES / "synthetic_interproscan.tsv").read_bytes()).hexdigest()
    with (out / "residue_annotations.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 69
    assert all(row["in_selected_hit"] == "" for row in rows if row["sequence_id"] == "toy_4")
    with (out / "per_sequence_coverage.csv").open(encoding="utf-8", newline="") as handle:
        coverage = list(csv.DictReader(handle))
    assert coverage[0]["covered_residue_count"] == "16"
    assert coverage[-1]["covered_residue_count"] == ""
    before = {p.name: p.read_bytes() for p in out.iterdir()}
    again = subprocess.run(command(out), cwd=ROOT, capture_output=True, text=True)
    assert again.returncode == 2
    assert {p.name: p.read_bytes() for p in out.iterdir()} == before


def test_cli_bad_input_fails_before_creating_output(tmp_path):
    out = tmp_path / "invalid output"
    bad_input = tmp_path / "bad.tsv"
    bad_input.write_text("bad\trow\n", encoding="utf-8")
    args = command(out)
    args[args.index("--interproscan-tsv") + 1] = str(bad_input)
    result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 2
    assert "TSV fields" in result.stderr
    assert not out.exists()


def test_cli_has_no_default_family(tmp_path):
    args = command(tmp_path / "no family")
    index = args.index("--pfam-accession")
    del args[index:index + 2]
    result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 2
    assert "--pfam-accession" in result.stderr
