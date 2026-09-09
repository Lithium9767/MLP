"""Analyze an existing protein FASTA alignment; never infer functional labels."""

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys


CANONICAL = frozenset("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS = frozenset("BJOUXZ")
MAP_FIELDS = ["sequence_id", "msa_position", "residue_position", "residue"]
SCORE_FIELDS = [
    "msa_position", "n_sequences", "n_canonical", "n_ambiguous", "n_gaps",
    "gap_fraction", "canonical_fraction", "consensus", "entropy_bits",
    "conservation", "sufficient_support",
]


def read_fasta(path, aligned=True):
    """First header token is the ID. Normalize case and '.' gaps; reject bad input."""
    records = {}
    identifier = None
    for number, raw in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            tokens = line[1:].split()
            if not tokens:
                raise ValueError(f"line {number}: empty FASTA identifier")
            identifier = tokens[0]
            if identifier in records:
                raise ValueError(f"line {number}: duplicate identifier {identifier}")
            records[identifier] = []
        else:
            if identifier is None:
                raise ValueError(f"line {number}: sequence before FASTA header")
            records[identifier].append(line.upper().replace(".", "-"))
    sequences = {key: "".join(parts) for key, parts in records.items()}
    if not sequences:
        raise ValueError("FASTA contains no sequences")
    allowed = CANONICAL | AMBIGUOUS | ({"-"} if aligned else set())
    for key, sequence in sequences.items():
        if not sequence or not sequence.strip("-"):
            raise ValueError(f"{key}: empty or entirely gapped sequence")
        invalid = set(sequence) - allowed
        if invalid:
            raise ValueError(f"{key}: invalid residues {sorted(invalid)}")
    if aligned and len({len(sequence) for sequence in sequences.values()}) != 1:
        raise ValueError("aligned sequences must have equal lengths")
    return sequences


def validate_originals(alignment, originals):
    if alignment.keys() != originals.keys():
        raise ValueError("original FASTA and alignment identifiers differ")
    for key, sequence in alignment.items():
        if sequence.replace("-", "") != originals[key]:
            raise ValueError(f"{key}: alignment does not reconstruct original sequence")


def coordinate_rows(alignment):
    """One row per MSA column per sequence. Gaps have no residue coordinate."""
    for key, sequence in alignment.items():
        position = 0
        for column, residue in enumerate(sequence, 1):
            if residue != "-":
                position += 1
            yield dict(sequence_id=key, msa_position=column,
                       residue_position=position if residue != "-" else None,
                       residue=residue)


def column_scores(alignment, min_canonical=2, min_fraction=0.5):
    """Unweighted Shannon entropy over canonical residues only, in bits."""
    if min_canonical < 2 or not 0 <= min_fraction <= 1:
        raise ValueError("min_canonical must be >= 2; min_fraction must be in [0, 1]")
    n = len(alignment)
    for position, column in enumerate(zip(*alignment.values()), 1):
        counts = Counter(residue for residue in column if residue in CANONICAL)
        canonical = sum(counts.values())
        gaps = column.count("-")
        entropy = -sum((count / canonical) * math.log2(count / canonical)
                       for count in counts.values()) if canonical else None
        entropy = max(0.0, entropy) if entropy is not None else None
        # All tied modes are retained, rather than arbitrarily naming one residue.
        consensus = "".join(sorted(r for r, count in counts.items()
                                   if count == max(counts.values()))) if counts else ""
        yield dict(
            msa_position=position, n_sequences=n, n_canonical=canonical,
            n_ambiguous=n - canonical - gaps, n_gaps=gaps,
            gap_fraction=gaps / n, canonical_fraction=canonical / n,
            consensus=consensus, entropy_bits=entropy,
            conservation=1 - entropy / math.log2(20) if entropy is not None else None,
            sufficient_support=canonical >= min_canonical and canonical / n >= min_fraction,
        )


def write_csv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_scores(scores, destination, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    x = [row["msa_position"] for row in scores]
    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True, layout="constrained")
    for ax, field, label, limit in zip(
        axes[:2], ["entropy_bits", "conservation"],
        ["Shannon entropy (bits)", "Conservation"], [math.log2(20), 1],
    ):
        values = [row[field] if row[field] is not None else math.nan for row in scores]
        ax.plot(x, values, color="#176b87", marker=".", markersize=4)
        low = [(row["msa_position"], row[field]) for row in scores
               if not row["sufficient_support"] and row[field] is not None]
        if low:
            ax.scatter(*zip(*low), color="#c05a27", marker="x", label="Low support", zorder=3)
            ax.legend(loc="upper right", fontsize=8)
        ax.set_ylabel(label)
        ax.set_ylim(-0.05 * limit, 1.1 * limit)
    axes[2].plot(x, [r["canonical_fraction"] for r in scores], label="Canonical", color="#176b87")
    axes[2].plot(x, [r["gap_fraction"] for r in scores], label="Gap", color="#c05a27")
    axes[2].set_ylabel("Column fraction")
    axes[2].set_ylim(-0.05, 1.1)
    axes[2].set_xlabel("MSA position (1-based)")
    axes[2].xaxis.set_major_locator(MaxNLocator(integer=True))
    axes[2].legend(loc="upper right", fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.suptitle(title)
    try:
        fig.savefig(destination / "conservation.png", dpi=160)
        fig.savefig(destination / "conservation.svg")
    finally:
        plt.close(fig)
    return matplotlib.__version__


def fingerprint(path):
    return {"path": str(Path(path).resolve()),
            "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}


def git_state():
    root = Path(__file__).resolve().parent.parent
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root,
                                capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=root,
                               capture_output=True, text=True, check=True).stdout.strip()
        return {"commit": commit, "dirty": bool(dirty)}
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "dirty": None}


def run(args):
    alignment = read_fasta(args.alignment)
    if args.original:
        validate_originals(alignment, read_fasta(args.original, aligned=False))
    scores = list(column_scores(alignment, args.min_canonical, args.min_fraction))
    # Fail before writing if plotting cannot be imported.
    import matplotlib  # noqa: F401
    destination = Path(args.out_dir)
    destination.mkdir(parents=True, exist_ok=False)
    try:
        write_csv(destination / "coordinate_map.csv", coordinate_rows(alignment), MAP_FIELDS)
        write_csv(destination / "column_scores.csv", scores, SCORE_FIELDS)
        version = plot_scores(scores, destination, args.title)
        files = sorted(destination.iterdir())
        manifest = {
            "schema_version": "1.0", "run_id": args.run_id,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "purpose": args.purpose, "dataset_version": args.dataset_version,
            "alignment": fingerprint(args.alignment),
            "original": fingerprint(args.original) if args.original else None,
            "original_verified": args.original is not None,
            "coordinate_convention": "1-based; gaps have empty residue_position",
            "n_sequences": len(alignment), "alignment_length": len(scores),
            "n_sufficient_support_columns": sum(r["sufficient_support"] for r in scores),
            "min_canonical": args.min_canonical, "min_fraction": args.min_fraction,
            "entropy_policy": "unweighted; canonical 20 only; log2; missing if no canonical residues",
            "conservation_formula": "1 - entropy_bits / log2(20)",
            "alignment_tool": args.alignment_tool,
            "alignment_tool_version": args.alignment_tool_version,
            "python": platform.python_version(), "matplotlib": version,
            "git": git_state(), "command_argv": sys.argv,
            "implementation": fingerprint(__file__),
            "outputs": {p.name: fingerprint(p)["sha256"] for p in files},
            "limitations": ["Descriptive statistics, not functional validation",
                            "No phylogenetic or sequence-redundancy weighting",
                            "Support cutoffs are tool settings, not frozen project evaluation rules"],
        }
        (destination / "manifest.json").write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    except Exception:
        (destination / "FAILED.txt").write_text(
            "Run incomplete. Do not use partial outputs. Rerun into a new directory.\n",
            encoding="utf-8")
        raise
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--original", type=Path)
    parser.add_argument("--out-dir", required=True, type=Path, help="New directory; never overwritten")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--dataset-version", required=True)
    parser.add_argument("--purpose", required=True, choices=["synthetic", "exploratory"])
    parser.add_argument("--alignment-tool", required=True)
    parser.add_argument("--alignment-tool-version", required=True)
    parser.add_argument("--min-canonical", type=int, default=2)
    parser.add_argument("--min-fraction", type=float, default=0.5)
    parser.add_argument("--title", default="Protein MSA descriptive statistics")
    args = parser.parse_args()
    try:
        manifest = run(args)
    except (ValueError, OSError, ImportError) as error:
        parser.exit(2, f"error: {error}\n")
    print(f"Completed: {manifest['n_sequences']} sequences, "
          f"{manifest['alignment_length']} columns -> {args.out_dir}")


if __name__ == "__main__":
    main()
