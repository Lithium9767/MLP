import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest
from bioinformatics.msa_analysis import read_fasta


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "bioinformatics" / "examples"


@pytest.mark.parametrize("format,alignment,original,total", [
    ("a3m", "synthetic.a3m", "synthetic_a3m_original.fasta", 18),
    ("fasta", "synthetic_alignment.fasta", "synthetic_original.fasta", 69),
])
def test_full_workflow(tmp_path, format, alignment, original, total):
    out = tmp_path / "joined result"
    command = [sys.executable, "-m", "bioinformatics.workflow", "--format", format,
               "--alignment", str(EXAMPLES / alignment), "--original", str(EXAMPLES / original),
               "--out-dir", str(out), "--run-id", "test", "--dataset-version", "synthetic-v1",
               "--purpose", "synthetic", "--alignment-tool", "hand-constructed",
               "--alignment-tool-version", "1", "--min-conservation", "0.9", "--min-length", "1",
               "--title", "Synthetic fixture - not a biological result"]
    if format == "fasta":
        tsv = EXAMPLES / "synthetic_interproscan.tsv"
    else:
        tsv = tmp_path / "a3m_hits.tsv"
        # Full query is MACDEEFG (8 residues). The hit spans inserted D/E and match E.
        tsv.write_text("query\t" + hashlib.md5(b"MACDEEFG").hexdigest()
                       + "\t8\tPfam\tPF00741\tSYNTHETIC\t4\t6\t1e-5\tT\t09-09-2026\n")
    command += ["--interproscan-tsv", str(tsv), "--pfam-accession", "PF00741",
                "--interproscan-version", "synthetic", "--pfam-version", "synthetic",
                "--source-description", "Synthetic, no database scan"]
    process = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert process.returncode == 0, process.stderr
    with (out / "residue_annotations.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == total
    for identifier, sequence in read_fasta(EXAMPLES / original, aligned=False).items():
        selected = [r for r in rows if r["sequence_id"] == identifier]
        assert [int(r["residue_position"]) for r in selected] == list(range(1, len(sequence) + 1))
        assert "".join(r["residue"] for r in selected) == sequence
    if format == "a3m":
        query = {int(row["residue_position"]): row for row in rows if row["sequence_id"] == "query"}
        assert query[6]["residue"] == "E"
        assert query[6]["alignment_position"] == "4"
        assert query[6]["in_selected_hit"] == "true"
        assert query[4]["in_selected_hit"] == "true"
        assert query[5]["in_selected_hit"] == "true"
        assert query[3]["in_selected_hit"] == "false"
        for position in [1, 4, 5, 8]:
            assert query[position]["state"] == "insert"
            assert query[position]["conservation"] == ""
            assert query[position]["alignment_position"] == ""
            assert query[position]["region_id"] == ""
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["original_verified"]
    assert manifest["pfam_joined"]
    for name, digest in manifest["outputs"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == digest
    before = (out / "manifest.json").read_bytes()
    again = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert again.returncode == 2
    assert (out / "manifest.json").read_bytes() == before
