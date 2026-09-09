"""Mathematical, coordinate, malformed-input and end-to-end checks; no network."""

import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import pytest

from bioinformatics.msa_analysis import (
    CANONICAL, column_scores, coordinate_rows, read_fasta, validate_originals,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "bioinformatics" / "examples"


def test_known_entropies():
    assert list(column_scores({"a": "A", "b": "A"}))[0]["entropy_bits"] == 0
    half = list(column_scores({"a": "A", "b": "C"}))[0]
    assert half["entropy_bits"] == 1
    assert half["consensus"] == "AC"
    assert half["conservation"] == pytest.approx(1 - 1 / math.log2(20))
    uniform = list(column_scores({r: r for r in CANONICAL}))[0]
    assert uniform["entropy_bits"] == pytest.approx(math.log2(20))
    assert uniform["conservation"] == pytest.approx(0)


def test_missing_and_ambiguous_do_not_fake_conservation():
    scores = list(column_scores({"a": "A-XA", "b": "--BC", "c": "--Z-"}))
    assert scores[0]["conservation"] == 1
    assert not scores[0]["sufficient_support"]
    assert scores[1]["entropy_bits"] is None
    assert scores[1]["conservation"] is None
    assert scores[1]["gap_fraction"] == 1
    assert scores[2]["n_ambiguous"] == 3
    assert scores[2]["entropy_bits"] is None
    assert scores[3]["entropy_bits"] == 1
    assert scores[3]["sufficient_support"]


def test_coordinates_roundtrip_with_terminal_internal_gaps_and_ambiguity():
    alignment = {"a": "-AX-C-", "b": "AC--DE"}
    rows = list(coordinate_rows(alignment))
    for key, aligned in alignment.items():
        subset = [row for row in rows if row["sequence_id"] == key]
        assert "".join(row["residue"] for row in subset) == aligned
        residues = [row for row in subset if row["residue_position"] is not None]
        assert [row["residue_position"] for row in residues] == list(range(1, len(residues) + 1))
        assert "".join(row["residue"] for row in residues) == aligned.replace("-", "")
        assert all(row["residue_position"] is None for row in subset if row["residue"] == "-")
    assert rows[2] == dict(sequence_id="a", msa_position=3, residue_position=2, residue="X")


@pytest.mark.parametrize("contents,match", [
    ("", "no sequences"), ("ACD", "before FASTA"), (">\nACD", "empty FASTA"),
    (">a\nAC\n>a description\nAC", "duplicate"), (">a\n>b\nAC", "empty"),
    (">a\n---", "gapped"), (">a\nAC\n>b\nA", "equal lengths"),
    (">a\nAC*", "invalid"), (">a\nA C", "invalid"), (">a\nA1C", "invalid"),
])
def test_malformed_fasta(tmp_path, contents, match):
    path = tmp_path / "input.fa"
    path.write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError, match=match):
        read_fasta(path)


def test_fasta_normalization_and_wrapping(tmp_path):
    path = tmp_path / "input.fa"
    path.write_text("\ufeff>a description\nac.\nx\n>b\nAC-X\n", encoding="utf-8")
    assert read_fasta(path) == {"a": "AC-X", "b": "AC-X"}
    with pytest.raises(ValueError, match="invalid"):
        read_fasta(path, aligned=False)


def test_rejects_a3m_as_plain_fasta_and_unicode_case_expansion(tmp_path):
    a3m = tmp_path / "input.a3m"
    a3m.write_text(">a\nAqC\n>b\nACq\n", encoding="utf-8")
    with pytest.raises(ValueError, match="A3M"):
        read_fasta(a3m)
    fasta = tmp_path / "input.fa"
    fasta.write_text(">a\nA\u00dfC\n", encoding="utf-8")
    with pytest.raises(ValueError, match="non-ASCII"):
        read_fasta(fasta)


@pytest.mark.parametrize("originals,match", [({"other": "AC"}, "identifiers"), ({"a": "CA"}, "reconstruct")])
def test_original_mismatch(originals, match):
    with pytest.raises(ValueError, match=match):
        validate_originals({"a": "A-C"}, originals)


@pytest.mark.parametrize("minimum,fraction", [(1, 0.5), (2, -0.1), (2, 1.1), (2, math.nan)])
def test_invalid_support_settings(minimum, fraction):
    with pytest.raises(ValueError):
        list(column_scores({"a": "A"}, minimum, fraction))


def command(out):
    return [sys.executable, "-m", "bioinformatics.msa_analysis",
            "--alignment", str(EXAMPLES / "synthetic_alignment.fasta"),
            "--original", str(EXAMPLES / "synthetic_original.fasta"),
            "--out-dir", str(out), "--run-id", "test-synthetic",
            "--dataset-version", "synthetic-v1", "--purpose", "synthetic",
            "--alignment-tool", "hand-constructed", "--alignment-tool-version", "1",
            "--title", "Synthetic fixture - not a biological result"]


def test_cli_artifacts_hashes_and_no_overwrite(tmp_path):
    out = tmp_path / "analysis with spaces"
    result = subprocess.run(command(out), cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["n_sequences"] == 4
    assert manifest["alignment_length"] == 20
    assert manifest["original_verified"]
    assert manifest["purpose"] == "synthetic"
    for name, digest in manifest["outputs"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == digest
    assert (out / "conservation.png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    with (out / "coordinate_map.csv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 80
    assert rows[3]["residue_position"] == ""
    with (out / "column_scores.csv").open(encoding="utf-8", newline="") as stream:
        scores = list(csv.DictReader(stream))
    assert scores[13]["entropy_bits"] == ""
    assert scores[18]["conservation"] == ""
    before = (out / "manifest.json").read_bytes()
    again = subprocess.run(command(out), cwd=ROOT, capture_output=True, text=True)
    assert again.returncode == 2
    assert (out / "manifest.json").read_bytes() == before


def test_cli_rejects_input_before_creating_outputs(tmp_path):
    path = tmp_path / "bad.fa"
    path.write_text(">a\nAC\n>b\nA", encoding="utf-8")
    out = tmp_path / "output"
    args = command(out)
    args[args.index("--alignment") + 1] = str(path)
    result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 2
    assert not out.exists()
