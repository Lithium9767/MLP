"""Run FASTA/A3M descriptive analysis and join scores to full residue coordinates."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from types import SimpleNamespace

from . import a3m, msa_analysis, pfam, regions


RESIDUE_FIELDS = ["sequence_id", "residue_position", "residue", "state",
                  "alignment_position", "insertion_anchor", "insertion_rank",
                  "entropy_bits", "conservation", "sufficient_support", "region_id"]
PFAM_FIELDS = ["pfam_accession", "pfam_annotation_status", "in_selected_hit", "n_covering_hits", "pfam_hit_ids"]


def residue_annotations(mapping, scores, intervals, pfam_rows=None):
    """Only explicit match positions receive column scores; insertions stay unscored."""
    by_column = {s["msa_position"]: s for s in scores}
    region_by_column = {p: r["region_id"] for r in intervals for p in range(r["start"], r["end"] + 1)}
    for row in mapping:
        column = row["alignment_position"]
        score = by_column[column] if column is not None else {}
        result = {**row, "entropy_bits": score.get("entropy_bits"),
               "conservation": score.get("conservation"),
               "sufficient_support": score.get("sufficient_support"),
               "region_id": region_by_column.get(column, "")}
        if pfam_rows is not None:
            annotation = pfam_rows[(row["sequence_id"], row["residue_position"])]
            if row["residue"] != annotation["residue"]:
                raise ValueError("Pfam/MSA residue identity differs at the joined coordinate")
            result.update(pfam_accession=annotation["pfam_accession"],
                          pfam_annotation_status=annotation["annotation_status"],
                          in_selected_hit=annotation["in_selected_hit"],
                          n_covering_hits=annotation["n_covering_hits"], pfam_hit_ids=annotation["hit_ids"])
        yield result


def run(args):
    # All biological sequence correspondence is checked before output creation.
    annotation_rows = None
    annotation_options = [args.interproscan_tsv, args.pfam_accession, args.interproscan_version,
                          args.pfam_version, args.source_description]
    if any(annotation_options):
        if not all(annotation_options) or args.original is None:
            raise ValueError("Pfam annotation requires --original and all five Pfam/source options")
        originals = msa_analysis.read_fasta(args.original, aligned=False)
        hits, audit = pfam.read_interproscan_tsv(args.interproscan_tsv, originals, args.pfam_accession)
        annotation_rows = {(r["sequence_id"], r["residue_position"]): r
                           for r in pfam.residue_annotations(originals, hits, audit, args.pfam_accession)}
    if args.format == "a3m":
        alignment = a3m.read_a3m(args.alignment)
        if args.original:
            a3m.validate_originals(alignment, msa_analysis.read_fasta(args.original, aligned=False))
        match_alignment = alignment.match_alignment
        mapping = ({"sequence_id": row["sequence_id"], "residue_position": row["original_position"],
                    "residue": row["residue"], "state": row["state"],
                    "alignment_position": row["match_position"],
                    "insertion_anchor": row["insertion_anchor"], "insertion_rank": row["insertion_rank"]}
                   for row in a3m.coordinate_rows(alignment))
    else:
        match_alignment = msa_analysis.read_fasta(args.alignment)
        if args.original:
            msa_analysis.validate_originals(match_alignment, msa_analysis.read_fasta(args.original, aligned=False))
        mapping = ({"sequence_id": row["sequence_id"], "residue_position": row["residue_position"],
                    "residue": row["residue"], "state": "match", "alignment_position": row["msa_position"],
                    "insertion_anchor": None, "insertion_rank": None}
                   for row in msa_analysis.coordinate_rows(match_alignment) if row["residue_position"] is not None)
    scores = list(msa_analysis.column_scores(match_alignment, args.min_canonical, args.min_fraction))
    region_scores = [{"position": row["msa_position"], "conservation": row["conservation"],
                      "sufficient_support": row["sufficient_support"]} for row in scores]
    intervals = regions.extract_regions(region_scores, args.min_conservation, args.min_length)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=False)
    try:
        if args.format == "a3m":
            a3m.run(SimpleNamespace(input=args.alignment, original=args.original,
                                   out_dir=out / "a3m", purpose=args.purpose))
            msa_input = out / "a3m" / "match_alignment.fasta"
        else:
            msa_input = args.alignment
        msa_analysis.run(SimpleNamespace(
            alignment=msa_input, original=args.original if args.format == "fasta" else None,
            out_dir=out / "msa", min_canonical=args.min_canonical, min_fraction=args.min_fraction,
            run_id=args.run_id, dataset_version=args.dataset_version, purpose=args.purpose,
            alignment_tool=args.alignment_tool, alignment_tool_version=args.alignment_tool_version,
            title=args.title))
        regions.run(SimpleNamespace(
            scores=out / "msa" / "column_scores.csv", out_dir=out / "regions",
            min_conservation=args.min_conservation, min_length=args.min_length,
            coordinate_system="a3m_match" if args.format == "a3m" else "msa",
            purpose=args.purpose, run_id=args.run_id))
        if annotation_rows is not None:
            pfam.run(SimpleNamespace(
                original=args.original, interproscan_tsv=args.interproscan_tsv,
                out_dir=out / "pfam", pfam_accession=args.pfam_accession,
                interproscan_version=args.interproscan_version, pfam_version=args.pfam_version,
                source_description=args.source_description, purpose=args.purpose,
                run_id=args.run_id, dataset_version=args.dataset_version))
        msa_analysis.write_csv(out / "residue_annotations.csv",
                               residue_annotations(mapping, scores, intervals, annotation_rows),
                               RESIDUE_FIELDS + (PFAM_FIELDS if annotation_rows is not None else []))
        manifest = {
            "schema_version": "1.0", "run_id": args.run_id, "purpose": args.purpose,
            "dataset_version": args.dataset_version, "format": args.format,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "inputs": {"alignment": msa_analysis.fingerprint(args.alignment),
                       "original": msa_analysis.fingerprint(args.original) if args.original else None},
            "original_verified": args.original is not None,
            "pfam_joined": annotation_rows is not None,
            "coordinate_authority": "residue_annotations.csv; 1-based full reconstructed residue positions",
            "n_sequences": len(match_alignment), "n_columns": len(scores), "n_regions": len(intervals),
            "git": msa_analysis.git_state(), "implementation": msa_analysis.fingerprint(__file__),
            "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
            "outputs": {str(p.relative_to(out)): msa_analysis.fingerprint(p)["sha256"]
                        for p in sorted(out.rglob("*")) if p.is_file()},
            "limitations": ["No functional labels or model attribution supplied",
                            "Full reconstructed sequence coordinates verified only if original FASTA supplied",
                            "A3M insertions retain full coordinates but have no match-column score",
                            "A3M msa/coordinate_map.csv uses match-only residue indices; use the root joined table",
                            "Conserved runs are descriptive, not statistical or functional discoveries"],
        }
        (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    except Exception:
        (out / "FAILED.txt").write_text("Incomplete run; do not use partial outputs.\n", encoding="utf-8")
        raise
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--format", choices=["fasta", "a3m"], required=True)
    parser.add_argument("--original", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--dataset-version", required=True)
    parser.add_argument("--purpose", choices=["synthetic", "exploratory"], required=True)
    parser.add_argument("--alignment-tool", required=True)
    parser.add_argument("--alignment-tool-version", required=True)
    parser.add_argument("--min-canonical", type=int, default=2)
    parser.add_argument("--min-fraction", type=float, default=0.5)
    parser.add_argument("--min-conservation", type=float, required=True)
    parser.add_argument("--min-length", type=int, required=True)
    parser.add_argument("--title", default="Protein alignment descriptive statistics")
    parser.add_argument("--interproscan-tsv", type=Path)
    parser.add_argument("--pfam-accession")
    parser.add_argument("--interproscan-version")
    parser.add_argument("--pfam-version")
    parser.add_argument("--source-description")
    args = parser.parse_args()
    try:
        result = run(args)
    except (ValueError, OSError, ImportError) as error:
        parser.exit(2, f"error: {error}\n")
    print(f"Completed: {result['n_sequences']} sequences, {result['n_regions']} regions -> {args.out_dir}")


if __name__ == "__main__":
    main()
