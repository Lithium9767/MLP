"""A3M state semantics, residue round trips, fail-closed validation and CLI checks."""

from argparse import Namespace
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from bioinformatics.a3m import coordinate_rows, read_a3m, run, sequence_qc, validate_originals
from bioinformatics.msa_analysis import coordinate_rows as match_only_coordinates
from bioinformatics.msa_analysis import read_fasta


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "bioinformatics" / "examples"
FIXTURE = EXAMPLES / "synthetic.a3m"
ORIGINAL = EXAMPLES / "synthetic_a3m_original.fasta"


def test_match_states_preserve_inserted_residues_and_audit_annotations():
    alignment = read_a3m(FIXTURE)
    assert alignment.match_alignment == {"query": "AC-EF", "homolog_one": "-CDEF", "homolog_two": "A-DE-"}
    assert alignment.full_sequences == {"query": "MACDEEFG", "homolog_one": "CTDEF", "homolog_two": "AQXDE"}
    assert alignment.match_columns == 5
    assert [row["sequence_id"] for row in alignment.excluded_records] == ["ss_pred", "ss_conf", "sa_pred", "aa_cons"]
    assert alignment.excluded_records[1]["content_sha256"] == hashlib.sha256(b"98765").hexdigest()
    assert alignment.comments[0]["text"].startswith("# SYNTHETIC")
    validate_originals(alignment, read_fasta(ORIGINAL, aligned=False))


def test_insertion_anchors_ranks_and_actual_original_positions():
    rows = [row for row in coordinate_rows(read_a3m(FIXTURE)) if row["sequence_id"] == "query"]
    assert [(row["original_position"], row["residue"], row["match_position"],
             row["insertion_anchor"], row["insertion_rank"]) for row in rows] == [
        (1, "M", None, 0, 1), (2, "A", 1, None, None), (3, "C", 2, None, None),
        (4, "D", None, 3, 1), (5, "E", None, 3, 2),
        (6, "E", 4, None, None), (7, "F", 5, None, None), (8, "G", None, 5, 1),
    ]
    assert [row["match_residue_position"] for row in rows] == [None, 1, 2, None, None, 3, 4, None]


def test_every_full_residue_roundtrips_and_match_only_bridge_is_explicit():
    alignment = read_a3m(FIXTURE)
    rows = list(coordinate_rows(alignment))
    for key, full in alignment.full_sequences.items():
        selected = [row for row in rows if row["sequence_id"] == key]
        assert "".join(row["residue"] for row in selected) == full
        assert [row["original_position"] for row in selected] == list(range(1, len(full) + 1))
        match_rows = [row for row in selected if row["state"] == "match"]
        old_rows = [row for row in match_only_coordinates({key: alignment.match_alignment[key]})
                    if row["residue_position"] is not None]
        assert [(r["match_position"], r["match_residue_position"], r["residue"]) for r in match_rows] == [
            (r["msa_position"], r["residue_position"], r["residue"]) for r in old_rows
        ]
        assert all(row["match_position"] is None and row["match_residue_position"] is None
                   for row in selected if row["state"] == "insert")
    # This is the critical coordinate trap: residue 3 of the match-only query
    # is residue 6 of the reconstructed full sequence, not original residue 3.
    query_e = next(r for r in rows if r["sequence_id"] == "query" and r["match_residue_position"] == 3)
    assert query_e["original_position"] == 6


def test_optional_dot_padding_does_not_change_coordinates(tmp_path):
    plain = tmp_path / "plain.a3m"
    padded = tmp_path / "padded.a3m"
    plain.write_text(">a\nmAC-d.eEFg\n", encoding="utf-8")
    padded.write_text(">a\n..m.AC..-d...eE.F..g.\n", encoding="utf-8")
    left, right = read_a3m(plain), read_a3m(padded)
    assert left.match_alignment == right.match_alignment
    assert left.full_sequences == right.full_sequences
    assert list(coordinate_rows(left)) == list(coordinate_rows(right))


def test_bom_wrapped_records_comments_and_first_header_token(tmp_path):
    path = tmp_path / "wrapped.a3m"
    path.write_text("\ufeff# title\n>sp|P12345|A description\nACd\n# inside record\nE-\n>b\nACe\nEF\n", encoding="utf-8")
    alignment = read_a3m(path)
    assert alignment.sequences == {"sp|P12345|A": "ACdE-", "b": "ACeEF"}
    assert alignment.headers["sp|P12345|A"] == "sp|P12345|A description"
    assert len(alignment.comments) == 2


@pytest.mark.parametrize("contents,match", [
    ("", "no protein sequences"), ("ACD", "before A3M"), (">\nACD", "empty A3M identifier"),
    (">a\nAC\n>a description\nAC", "duplicate"),
    (">ss_conf\n99\n>ss_conf\n98\n>a\nAC", "duplicate"),
    (">a>bad\nAC", "invalid A3M identifier"), (">a\x00\nAC", "invalid A3M identifier"),
    (">a\n>b\nAC", "empty A3M sequence"), (">a\nac", "no match columns"),
    (">a\n---", "no match residues"), (">a\na-b-c", "no match residues"),
    (">a\nACdE\n>b\nACEF", "match column count"),
    (">a\nAC*", "invalid A3M"), (">a\nA C", "invalid A3M"),
    (">a\nAC1", "invalid A3M"), (">a\nACé", "invalid A3M"),
    (">ss_conf\n99", "no protein sequences"),
])
def test_malformed_inputs_are_rejected(tmp_path, contents, match):
    path = tmp_path / "bad.a3m"
    path.write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError, match=match):
        read_a3m(path)


@pytest.mark.parametrize("mutation,match", [("missing", "identifiers differ"), ("extra", "identifiers differ"),
                                             ("match_only", "including inserts"), ("substitution", "including inserts")])
def test_original_validation_never_confuses_match_only_with_full_sequence(mutation, match):
    alignment = read_a3m(FIXTURE)
    originals = dict(alignment.full_sequences)
    if mutation == "missing":
        originals.pop("query")
    elif mutation == "extra":
        originals["extra"] = "ACD"
    elif mutation == "match_only":
        originals["query"] = alignment.match_alignment["query"].replace("-", "")
    else:
        originals["query"] = "MACDEEFW"
    with pytest.raises(ValueError, match=match):
        validate_originals(alignment, originals)


def test_sequence_quality_counts_include_insertions_and_ambiguous_residues():
    rows = list(sequence_qc(read_a3m(FIXTURE)))
    assert rows[0] == dict(sequence_id="query", match_columns=5, full_length=8, match_residues=4,
                           inserted_residues=4, match_gaps=1, insert_gaps=1, ambiguous_residues=0)
    assert rows[2]["ambiguous_residues"] == 1
    assert sum(row["full_length"] for row in rows) == 18
    assert sum(row["inserted_residues"] for row in rows) == 7


def command(out):
    return [sys.executable, "-m", "bioinformatics.a3m", "--input", str(FIXTURE),
            "--original", str(ORIGINAL), "--out-dir", str(out), "--purpose", "synthetic"]


def test_cli_outputs_are_reconstructable_hashed_and_never_overwritten(tmp_path):
    out = tmp_path / "conversion with spaces"
    result = subprocess.run(command(out), cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["n_sequences"] == 3
    assert manifest["match_columns"] == 5
    assert manifest["n_full_residues"] == 18
    assert manifest["n_inserted_residues"] == 7
    assert manifest["n_excluded_records"] == 4
    assert manifest["synthetic"] and manifest["original_verified"]
    assert manifest["input"]["sha256"] == hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
    assert manifest["original"]["sha256"] == hashlib.sha256(ORIGINAL.read_bytes()).hexdigest()
    for name, digest in manifest["outputs"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == digest
    assert read_fasta(out / "match_alignment.fasta") == read_a3m(FIXTURE).match_alignment
    assert read_fasta(out / "full_sequences.fasta", aligned=False) == read_a3m(FIXTURE).full_sequences
    with (out / "coordinate_map.csv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 18
    assert rows[0]["match_position"] == ""
    assert rows[0]["insertion_anchor"] == "0"
    before = {path.name: path.read_bytes() for path in out.iterdir()}
    repeated = subprocess.run(command(out), cwd=ROOT, capture_output=True, text=True)
    assert repeated.returncode == 2
    assert {path.name: path.read_bytes() for path in out.iterdir()} == before


def test_no_original_does_not_claim_verification_or_biological_validation(tmp_path):
    manifest = run(Namespace(input=FIXTURE, original=None, out_dir=tmp_path / "unchecked", purpose="synthetic"))
    assert not manifest["original_verified"]
    assert manifest["original"] is None
    assert any("untrimmed" in note for note in manifest["limitations"])


def test_cli_original_mismatch_fails_before_creating_outputs(tmp_path):
    wrong_original = tmp_path / "wrong.fasta"
    wrong_original.write_text(">query\nACEF\n>homolog_one\nCTDEF\n>homolog_two\nAQXDE\n", encoding="utf-8")
    out = tmp_path / "invalid output"
    args = command(out)
    args[args.index("--original") + 1] = str(wrong_original)
    result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 2
    assert "including inserts" in result.stderr
    assert not out.exists()


def test_cli_malformed_a3m_fails_before_creating_outputs(tmp_path):
    path = tmp_path / "bad.a3m"
    path.write_text(">a\nAC\n>b\nAc\n", encoding="utf-8")
    out = tmp_path / "invalid alignment"
    args = command(out)
    args[args.index("--input") + 1] = str(path)
    result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 2
    assert "match column count" in result.stderr
    assert not out.exists()
