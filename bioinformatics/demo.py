"""Reproduce the complete C tool acceptance run using synthetic fixtures only."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

from .msa_analysis import fingerprint, git_state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True, help="New directory; never overwritten")
    args = parser.parse_args()
    out = args.out_dir.resolve()
    root = Path(__file__).resolve().parents[1]
    examples = root / "bioinformatics" / "examples"
    try:
        out.mkdir(parents=True, exist_ok=False)
    except OSError as error:
        parser.exit(2, f"error: {error}\n")
    common = ["--purpose", "synthetic", "--dataset-version", "synthetic-v1",
              "--alignment-tool", "hand-constructed", "--alignment-tool-version", "1",
              "--min-conservation", "0.9", "--min-length", "2",
              "--title", "Synthetic fixture - not a biological result"]
    commands = [
        [sys.executable, "-m", "bioinformatics.data_audit",
         "--fasta", str(examples / "synthetic_original.fasta"),
         "--metadata", str(examples / "synthetic_metadata.csv"), "--out-dir", str(out / "data_audit")],
        [sys.executable, "-m", "bioinformatics.workflow", "--format", "fasta",
         "--alignment", str(examples / "synthetic_alignment.fasta"),
         "--original", str(examples / "synthetic_original.fasta"),
         "--out-dir", str(out / "fasta_pfam"), "--run-id", "C-demo-fasta-pfam", *common,
         "--interproscan-tsv", str(examples / "synthetic_interproscan.tsv"),
         "--pfam-accession", "PF00741", "--pfam-version", "synthetic-not-scanned",
         "--interproscan-version", "synthetic-not-run",
         "--source-description", "Hand-constructed synthetic fixture; not a biological result"],
        [sys.executable, "-m", "bioinformatics.workflow", "--format", "a3m",
         "--alignment", str(examples / "synthetic.a3m"),
         "--original", str(examples / "synthetic_a3m_original.fasta"),
         "--out-dir", str(out / "a3m"), "--run-id", "C-demo-a3m", *common],
    ]
    try:
        for command in commands:
            subprocess.run(command, cwd=root, check=True)
        manifest = {"purpose": "synthetic acceptance only; not a biological experiment",
                    "created_utc": datetime.now(timezone.utc).isoformat(), "git": git_state(),
                    "commands": commands,
                    "outputs": {str(p.relative_to(out)): fingerprint(p)["sha256"]
                                for p in sorted(out.rglob("*")) if p.is_file()}}
        (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    except Exception:
        (out / "FAILED.txt").write_text("Incomplete synthetic acceptance run.\n", encoding="utf-8")
        raise
    print(f"All synthetic acceptance steps completed -> {out}")


if __name__ == "__main__":
    main()
