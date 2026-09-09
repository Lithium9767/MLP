"""Extract contiguous conserved columns under explicit descriptive thresholds."""

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path

from .msa_analysis import fingerprint, git_state, write_csv

REGION_FIELDS = ["region_id", "start", "end", "length", "min_conservation", "mean_conservation"]


def read_scores(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError("score table has duplicate column names")
        required = {"msa_position", "conservation", "sufficient_support"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f"scores require columns {sorted(required)}")
        rows = []
        for number, row in enumerate(reader, 1):
            if None in row or any(value is None for value in row.values()):
                raise ValueError("score table row has the wrong field count")
            if int(row["msa_position"]) != number:
                raise ValueError("score positions must be consecutive, unique and start at 1")
            support = row["sufficient_support"].lower()
            if support not in {"true", "false"}:
                raise ValueError(f"position {number}: invalid support flag")
            score = float(row["conservation"]) if row["conservation"] else None
            if score is not None and (not math.isfinite(score) or not -1e-12 <= score <= 1 + 1e-12):
                raise ValueError(f"position {number}: conservation must be finite in [0, 1]")
            if support == "true" and score is None:
                raise ValueError(f"position {number}: supported column has no score")
            rows.append({"position": number, "conservation": score,
                         "sufficient_support": support == "true"})
    if not rows:
        raise ValueError("score table is empty")
    return rows


def extract_regions(scores, min_conservation, min_length):
    if not math.isfinite(min_conservation) or not 0 <= min_conservation <= 1:
        raise ValueError("min_conservation must be in [0, 1]")
    if not isinstance(min_length, int) or min_length < 1:
        raise ValueError("min_length must be positive")
    regions, run = [], []

    def finish():
        if len(run) >= min_length:
            values = [r["conservation"] for r in run]
            regions.append({"region_id": f"region_{len(regions) + 1}",
                            "start": run[0]["position"], "end": run[-1]["position"],
                            "length": len(run), "min_conservation": min(values),
                            "mean_conservation": sum(values) / len(values)})
        run.clear()

    for position, row in enumerate(scores, 1):
        if row["position"] != position:
            raise ValueError("score positions must be consecutive, unique and start at 1")
        if (row["sufficient_support"] and row["conservation"] is not None
                and row["conservation"] >= min_conservation):
            run.append(row)
        else:
            finish()
    finish()
    return regions


def run(args):
    scores = read_scores(args.scores)
    regions = extract_regions(scores, args.min_conservation, args.min_length)
    memberships = {position: region["region_id"] for region in regions
                   for position in range(region["start"], region["end"] + 1)}
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=False)
    try:
        write_csv(out / "regions.csv", regions, REGION_FIELDS)
        write_csv(out / "column_membership.csv",
                  ({**row, "region_id": memberships.get(row["position"], "")} for row in scores),
                  ["position", "conservation", "sufficient_support", "region_id"])
        manifest = {"schema_version": "1.0", "created_utc": datetime.now(timezone.utc).isoformat(),
                    "run_id": args.run_id, "purpose": args.purpose,
                    "coordinate_system": args.coordinate_system,
                    "coordinate_convention": "1-based inclusive; never original sequence positions",
                    "min_conservation": args.min_conservation, "min_length": args.min_length,
                    "n_columns": len(scores), "n_regions": len(regions),
                    "input": fingerprint(args.scores), "implementation": fingerprint(__file__),
                    "outputs": {p.name: fingerprint(p)["sha256"] for p in out.iterdir()},
                    "git": git_state(),
                    "limitations": ["Descriptive contiguous runs, not functional regions or significance tests",
                                    "Low-support and missing columns split runs; no gap bridging",
                                    "Thresholds are explicit tool settings, not frozen project rules"]}
        (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    except Exception:
        (out / "FAILED.txt").write_text("Incomplete run; do not use partial outputs.\n", encoding="utf-8")
        raise
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--min-conservation", type=float, required=True)
    parser.add_argument("--min-length", type=int, required=True)
    parser.add_argument("--coordinate-system", choices=["msa", "a3m_match"], required=True)
    parser.add_argument("--purpose", choices=["synthetic", "exploratory"], required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    try:
        result = run(args)
    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")
    print(f"Completed: {result['n_regions']} descriptive regions -> {args.out_dir}")


if __name__ == "__main__":
    main()
