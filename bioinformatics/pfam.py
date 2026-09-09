"""Validate reported InterProScan 5 Pfam hits and map them to original residues.

This tool parses existing annotations; it does not run a database search or infer
functional labels. Format reference (accessed 2026-09-09):
https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html
"""

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import platform
import re
import sys

from bioinformatics.msa_analysis import fingerprint, git_state, read_fasta, write_csv


FORMAT_SOURCE = "https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html"
HIT_FIELDS = [
    "hit_id", "source_line", "sequence_id", "sequence_md5", "sequence_length",
    "analysis", "signature_accession", "pfam_accession", "signature_description",
    "start", "end", "score", "status", "run_date", "interpro_accession",
    "interpro_description", "optional_columns_json",
]
RESIDUE_FIELDS = [
    "sequence_id", "residue_position", "residue_fraction", "residue", "pfam_accession",
    "annotation_status", "in_selected_hit", "n_covering_hits", "hit_ids",
]
COVERAGE_FIELDS = [
    "sequence_id", "sequence_length", "pfam_accession", "annotation_status",
    "n_selected_hits", "covered_residue_count", "reported_interval_coverage_fraction",
]
LIMITATIONS = [
    "Reported profile hits are not experimental functional labels.",
    "TSV contains matches only; absence of a row does not establish whether a sequence was scanned.",
    "No selected-family hit is uncertain, including when another analysis has a reported match.",
    "false means outside reported intervals on a sequence with a selected-family hit; it is not evidence of absent function.",
    "Versions and source description are supplied by the caller; TSV does not prove database versions or scan completeness.",
    "Intervals use original protein coordinates, not MSA columns or model token positions.",
    "TSV start/end span is preserved; internal discontinuous fragments, if any, require richer source formats.",
]


def validate_accession(accession):
    if re.fullmatch(r"PF[0-9]{5}", accession) is None:
        raise ValueError("selected Pfam accession must have the explicit form PFxxxxx (five digits)")
    return accession


def _positive_integer(value, field, line):
    if re.fullmatch(r"[1-9][0-9]*", value) is None:
        raise ValueError(f"line {line}: {field} must be a positive integer")
    return int(value)


def read_interproscan_tsv(path, originals, pfam_accession):
    """Return selected hits and an audit after validating every nonempty row.

    Standard base rows have 11 columns; an InterPro accession/description pair
    extends them to 13. Up to two optional annotation columns extend them to 15.
    Optional annotation text is preserved without reinterpreting source options.
    """
    validate_accession(pfam_accession)
    if not originals:
        raise ValueError("original FASTA must contain sequences")
    expected_md5 = {key: hashlib.md5(seq.encode("ascii")).hexdigest()
                    for key, seq in originals.items()}
    hits = []
    reported_ids = set()
    analyses = Counter()
    total_rows = 0
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        for line_number, raw in enumerate(handle, 1):
            text = raw.rstrip("\r\n")
            if not text.strip():
                continue
            columns = text.split("\t")
            if len(columns) not in {11, 13, 14, 15}:
                raise ValueError(f"line {line_number}: expected 11 or 13-15 TSV fields; got {len(columns)}")
            if any(not value.strip() for value in columns):
                raise ValueError(f"line {line_number}: empty TSV field; use '-' for missing values")
            if any(value != value.strip() for value in columns):
                raise ValueError(f"line {line_number}: leading/trailing whitespace in TSV field")
            key, md5, length_text, analysis, signature, description = columns[:6]
            if key not in originals:
                raise ValueError(f"line {line_number}: unknown sequence_id {key!r}")
            if re.fullmatch(r"[0-9a-fA-F]{32}", md5) is None:
                raise ValueError(f"line {line_number}: invalid sequence MD5 digest")
            if md5.lower() != expected_md5[key]:
                raise ValueError(f"line {line_number}: sequence MD5 mismatch for {key}")
            length = _positive_integer(length_text, "sequence length", line_number)
            if length != len(originals[key]):
                raise ValueError(f"line {line_number}: sequence length mismatch for {key}")
            start = _positive_integer(columns[6], "start", line_number)
            end = _positive_integer(columns[7], "end", line_number)
            if not 1 <= start <= end <= length:
                raise ValueError(f"line {line_number}: require 1 <= start <= end <= sequence length")
            if analysis == "-" or signature == "-" or any(c.isspace() for c in analysis + signature):
                raise ValueError(f"line {line_number}: invalid analysis or signature accession")
            canonical_accession = ""
            if analysis.casefold() == "pfam":
                if re.fullmatch(r"PF[0-9]{5}(?:\.[0-9]+)?", signature) is None:
                    raise ValueError(f"line {line_number}: invalid Pfam signature accession")
                canonical_accession = signature.split(".", 1)[0]
            score = columns[8]
            if score != "-":
                if re.fullmatch(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?", score) is None:
                    raise ValueError(f"line {line_number}: invalid match score")
                try:
                    numeric_score = Decimal(score)
                except InvalidOperation as error:
                    raise ValueError(f"line {line_number}: invalid match score") from error
                if not numeric_score.is_finite() or (canonical_accession and numeric_score < 0):
                    raise ValueError(f"line {line_number}: invalid match score")
            if columns[9] != "T":
                raise ValueError(f"line {line_number}: expected InterProScan match status T")
            if re.fullmatch(r"[0-9]{2}-[0-9]{2}-[0-9]{4}", columns[10]) is None:
                raise ValueError(f"line {line_number}: expected run date DD-MM-YYYY")
            try:
                datetime.strptime(columns[10], "%d-%m-%Y")
            except ValueError as error:
                raise ValueError(f"line {line_number}: invalid run date") from error
            interpro_accession, interpro_description = (columns[11:13] if len(columns) >= 13 else ["-", "-"])
            if interpro_accession != "-" and re.fullmatch(r"IPR[0-9]{6}", interpro_accession) is None:
                raise ValueError(f"line {line_number}: invalid InterPro accession")
            if interpro_accession == "-" and interpro_description != "-":
                raise ValueError(f"line {line_number}: InterPro description without accession")
            total_rows += 1
            reported_ids.add(key)
            analyses[analysis] += 1
            if canonical_accession != pfam_accession:
                continue
            hits.append(dict(
                hit_id=f"tsv-line-{line_number}", source_line=line_number,
                sequence_id=key, sequence_md5=md5.lower(), sequence_length=length,
                analysis=analysis, signature_accession=signature, pfam_accession=canonical_accession,
                signature_description=description, start=start, end=end, score=score,
                status=columns[9], run_date=columns[10], interpro_accession=interpro_accession,
                interpro_description=interpro_description,
                optional_columns_json=json.dumps(columns[13:], ensure_ascii=False),
            ))
    selected_ids = {hit["sequence_id"] for hit in hits}
    audit = {
        "n_original_sequences": len(originals), "n_input_rows": total_rows,
        "n_selected_hits": len(hits), "n_filtered_rows": total_rows - len(hits),
        "all_sequences_reported_hit_fraction": len(selected_ids) / len(originals),
        "reported_hit_fraction_interpretation": "fraction of input FASTA IDs with a selected-family hit in this report; not a validated scan detection rate",
        "analysis_row_counts": dict(sorted(analyses.items())),
        "selected_hit_reported": sorted(selected_ids),
        "no_selected_hit_uncertain": sorted(reported_ids - selected_ids),
        "not_reported_uncertain": sorted(originals.keys() - reported_ids),
        "all_input_rows_validated": True, "sequence_md5_and_length_verified": True,
        "duplicate_policy": "Retain all rows including identical or overlapping hits; hit_id identifies source line.",
        "limitations": LIMITATIONS,
    }
    return hits, audit


def residue_annotations(originals, hits, audit, pfam_accession):
    """Use 1-based inclusive intervals; unknown coverage is empty, never false."""
    by_sequence = defaultdict(list)
    for hit in hits:
        by_sequence[hit["sequence_id"]].append(hit)
    reported_without_selected = set(audit["no_selected_hit_uncertain"])
    for key, sequence in originals.items():
        selected = by_sequence[key]
        status = ("selected_hit_reported" if selected else
                  "no_selected_hit_uncertain" if key in reported_without_selected else
                  "not_reported_uncertain")
        for position, residue in enumerate(sequence, 1):
            covering = [hit["hit_id"] for hit in selected if hit["start"] <= position <= hit["end"]]
            yield dict(
                sequence_id=key, residue_position=position, residue=residue,
                residue_fraction=position / len(sequence),
                pfam_accession=pfam_accession, annotation_status=status,
                in_selected_hit=("true" if covering else "false") if selected else "",
                n_covering_hits=len(covering) if selected else None,
                hit_ids="|".join(covering),
            )


def sequence_coverage(originals, hits, audit, pfam_accession):
    """Count the union of reported intervals; missing selected hits remain unknown."""
    by_sequence = defaultdict(list)
    for hit in hits:
        by_sequence[hit["sequence_id"]].append(hit)
    reported_without_selected = set(audit["no_selected_hit_uncertain"])
    for key, sequence in originals.items():
        selected = by_sequence[key]
        covered = {position for hit in selected for position in range(hit["start"], hit["end"] + 1)}
        status = ("selected_hit_reported" if selected else
                  "no_selected_hit_uncertain" if key in reported_without_selected else
                  "not_reported_uncertain")
        yield dict(sequence_id=key, sequence_length=len(sequence), pfam_accession=pfam_accession,
                   annotation_status=status, n_selected_hits=len(selected),
                   covered_residue_count=len(covered) if selected else None,
                   reported_interval_coverage_fraction=len(covered) / len(sequence) if selected else None)


def _write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def run(args):
    for field in ("run_id", "dataset_version", "interproscan_version", "pfam_version", "source_description"):
        if not getattr(args, field).strip():
            raise ValueError(f"{field} must not be empty")
    originals = read_fasta(args.original, aligned=False)
    hits, audit = read_interproscan_tsv(args.interproscan_tsv, originals, args.pfam_accession)
    destination = Path(args.out_dir)
    destination.mkdir(parents=True, exist_ok=False)
    try:
        write_csv(destination / "hits.csv", hits, HIT_FIELDS)
        write_csv(destination / "residue_annotations.csv",
                  residue_annotations(originals, hits, audit, args.pfam_accession), RESIDUE_FIELDS)
        write_csv(destination / "per_sequence_coverage.csv",
                  sequence_coverage(originals, hits, audit, args.pfam_accession), COVERAGE_FIELDS)
        _write_json(destination / "audit.json", audit)
        manifest = {
            "schema_version": "1.0", "run_id": args.run_id,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "purpose": args.purpose, "dataset_version": args.dataset_version,
            "original": fingerprint(args.original), "interproscan_tsv": fingerprint(args.interproscan_tsv),
            "pfam_accession": args.pfam_accession, "interproscan_version": args.interproscan_version,
            "pfam_version": args.pfam_version, "source_description": args.source_description,
            "version_provenance": "caller supplied; not recoverable or independently verified from TSV",
            "format_source": FORMAT_SOURCE, "n_sequences": len(originals), "n_selected_hits": len(hits),
            "coordinate_convention": "original protein; 1-based inclusive start/end and residue_position",
            "residue_fraction_formula": "residue_position / original_sequence_length; (0, 1]",
            "coverage_formula": "union of selected reported inclusive intervals / original_sequence_length; unknown for sequences without selected hits",
            "join_keys": ["sequence_id", "residue_position"],
            "python": platform.python_version(), "git": git_state(), "command_argv": sys.argv,
            "implementation": fingerprint(__file__),
            "shared_implementation": fingerprint(Path(__file__).with_name("msa_analysis.py")),
            "outputs": {p.name: fingerprint(p)["sha256"] for p in sorted(destination.iterdir())},
            "limitations": LIMITATIONS,
        }
        _write_json(destination / "manifest.json", manifest)
    except Exception:
        (destination / "FAILED.txt").write_text(
            "Run incomplete. Do not use partial outputs. Rerun into a new directory.\n", encoding="utf-8")
        raise
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", required=True, type=Path, help="Original ungapped protein FASTA")
    parser.add_argument("--interproscan-tsv", required=True, type=Path)
    parser.add_argument("--pfam-accession", required=True, help="Explicit PFxxxxx; no assumed family")
    parser.add_argument("--out-dir", required=True, type=Path, help="New directory; never overwritten")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--dataset-version", required=True)
    parser.add_argument("--purpose", required=True, choices=["synthetic", "exploratory"])
    parser.add_argument("--interproscan-version", required=True)
    parser.add_argument("--pfam-version", required=True)
    parser.add_argument("--source-description", required=True, help="Origin of TSV, scan options, or synthetic provenance")
    args = parser.parse_args()
    try:
        manifest = run(args)
    except (ValueError, OSError) as error:
        parser.exit(2, f"error: {error}\n")
    print(f"Completed: {manifest['n_sequences']} sequences, {manifest['n_selected_hits']} reported hits -> {args.out_dir}")


if __name__ == "__main__":
    main()
