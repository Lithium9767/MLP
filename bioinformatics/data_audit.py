"""Audit B's FASTA/metadata handoff without assigning or inferring labels."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys


ALLOWED_AA = frozenset("ACDEFGHIKLMNPQRSTVWYXBZJUO")
CORE_FIELDS = {"sequence_id", "sequence", "sequence_length", "sequence_hash", "exclusion_reason"}
SOURCE_FIELDS = ("accession", "database", "source_url", "retrieval_date",
                 "organism", "taxonomy_id", "lineage")
EVIDENCE_FIELDS = ("label_definition", "label_source", "evidence_level")
ALLOWED_SPLITS = frozenset({"train", "val", "test"})
PENDING_VALUES = frozenset({"", "todo", "unknown", "pending", "none", "n/a", "na",
                            "unlabelled_sequence", "unlabeled_sequence",
                            "pending_prediction_target_confirmation"})


def _issue(target, code, message, **details):
    target.append(dict(code=code, message=message, **details))


def _read_csv(path, required):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError("CSV has duplicate column names")
        missing = required - set(fields)
        if missing:
            raise ValueError(f"CSV missing required columns: {sorted(missing)}")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"CSV row ending at line {reader.line_num} has wrong field count")
            rows.append({key: value.strip() for key, value in row.items()})
        return rows, fields


def _read_fasta(path):
    """Retain duplicate IDs so they can be reported; never collapse records."""
    records = []
    identifier, parts = None, []
    for number, raw in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                records.append((identifier, "".join(parts)))
            tokens = line[1:].split()
            if not tokens:
                raise ValueError(f"line {number}: empty FASTA identifier")
            identifier, parts = tokens[0], []
        elif identifier is None:
            raise ValueError(f"line {number}: sequence before FASTA header")
        else:
            parts.append(line)
    if identifier is not None:
        records.append((identifier, "".join(parts)))
    return records


def _normalized_sequence(sequence):
    # Check ASCII before upper(): Unicode such as 'ß' would otherwise become 'SS'.
    if not sequence or not sequence.isascii():
        raise ValueError("sequence must contain nonempty ASCII amino-acid letters")
    sequence = sequence.upper()
    invalid = set(sequence) - ALLOWED_AA
    if invalid:
        raise ValueError(f"invalid residues: {sorted(invalid)}")
    return sequence


def _pending(value):
    normalized = value.strip().lower()
    return normalized in PENDING_VALUES or normalized.startswith("pending_")


def audit_data(fasta, metadata, split=None):
    """Return JSON-compatible findings. Missing labels do not make format invalid."""
    report = {
        "schema_version": "1.0", "created_utc": datetime.now(timezone.utc).isoformat(),
        "hash_policy": "SHA256 of nonempty ASCII amino-acid sequence uppercased; no gaps or whitespace",
        "inputs": {}, "errors": [], "warnings": [], "summary": {},
        "label_policy": "Labels are never inferred; a missing label is unknown, never a negative example.",
        "limitations": [
            "Presence of a source or evidence field does not verify its scientific truth or license.",
            "A passing audit is not project approval, a frozen dataset, or supervised model readiness.",
            "Cluster leakage checks validate supplied assignments; they do not compute sequence similarity.",
            "Excluded metadata records are counted but their excluded sequences need not occur in FASTA.",
        ],
    }
    errors, warnings = report["errors"], report["warnings"]
    for name, path in (("fasta", fasta), ("metadata", metadata), ("split", split)):
        if path is not None:
            try:
                report["inputs"][name] = {
                    "path": str(Path(path).resolve()),
                    "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                }
            except OSError as error:
                _issue(errors, "input_unreadable", str(error), input=name)
    records, rows, fields, split_rows = [], [], [], None
    try:
        records = _read_fasta(fasta)
    except (OSError, UnicodeError, ValueError) as error:
        _issue(errors, "fasta_parse_error", str(error))
    try:
        rows, fields = _read_csv(metadata, CORE_FIELDS)
    except (OSError, UnicodeError, ValueError, csv.Error) as error:
        _issue(errors, "metadata_parse_error", str(error))
    if split is not None:
        try:
            split_rows, _ = _read_csv(split, {"sequence_id", "similarity_cluster", "split"})
        except (OSError, UnicodeError, ValueError, csv.Error) as error:
            _issue(errors, "split_parse_error", str(error))

    included = [(number, row) for number, row in enumerate(rows, 2) if not row["exclusion_reason"]]
    excluded = [row for row in rows if row["exclusion_reason"]]
    report["summary"].update(metadata_records=len(rows), included_records=len(included),
                              excluded_records=len(excluded), fasta_records=len(records))
    if not included:
        _issue(errors, "no_included_records", "No valid metadata records are available for analysis")
    if not records:
        _issue(errors, "empty_fasta", "FASTA contains no records")

    fasta_counts = Counter(identifier for identifier, _ in records)
    metadata_counts = Counter(row["sequence_id"] for _, row in included)
    for name, counts in (("fasta", fasta_counts), ("metadata", metadata_counts)):
        duplicates = sorted(key for key, count in counts.items() if count > 1)
        if duplicates:
            _issue(errors, f"duplicate_{name}_ids", "Included IDs must be unique", sequence_ids=duplicates)
    if "" in metadata_counts:
        _issue(errors, "missing_sequence_id", "Included metadata has empty sequence_id")
    for code, identifiers in (
        ("fasta_missing_ids", set(metadata_counts) - set(fasta_counts)),
        ("fasta_extra_ids", set(fasta_counts) - set(metadata_counts)),
    ):
        if identifiers:
            _issue(errors, code, "FASTA and included metadata must cover exactly the same IDs",
                   sequence_ids=sorted(identifiers))
    excluded_collisions = sorted({row["sequence_id"] for row in excluded} & set(metadata_counts))
    if excluded_collisions:
        _issue(warnings, "excluded_id_collisions", "Excluded records reuse an included ID; retain row-level provenance",
               sequence_ids=excluded_collisions)

    fasta_sequences, metadata_sequences = {}, {}
    for identifier, raw_sequence in records:
        try:
            fasta_sequences[identifier] = _normalized_sequence(raw_sequence)
        except ValueError as error:
            _issue(errors, "invalid_fasta_sequence", str(error), sequence_id=identifier)
    lengths, hash_groups = [], defaultdict(list)
    for number, row in included:
        identifier = row["sequence_id"]
        try:
            sequence = _normalized_sequence(row["sequence"])
        except ValueError as error:
            _issue(errors, "invalid_metadata_sequence", str(error), sequence_id=identifier, row=number)
            continue
        metadata_sequences[identifier] = sequence
        lengths.append(len(sequence))
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        hash_groups[digest].append(identifier)
        if row["sequence_hash"] != digest:
            _issue(errors, "sequence_hash_mismatch", "sequence_hash must equal uppercase-ASCII SHA256",
                   sequence_id=identifier, row=number, expected=digest)
        try:
            valid_length = int(row["sequence_length"]) == len(sequence)
        except ValueError:
            valid_length = False
        if not valid_length:
            _issue(errors, "sequence_length_mismatch", "sequence_length does not match sequence",
                   sequence_id=identifier, row=number, expected=len(sequence))
        if "valid_residues" in fields and row["valid_residues"].lower() != "true":
            _issue(errors, "invalid_residue_flag", "Included record must have valid_residues=true",
                   sequence_id=identifier, row=number)
        if identifier in fasta_sequences and fasta_sequences[identifier] != sequence:
            _issue(errors, "sequence_mismatch", "Complete FASTA and metadata sequences differ",
                   sequence_id=identifier, row=number)
    duplicate_sequences = [sorted(ids) for ids in hash_groups.values() if len(set(ids)) > 1]
    if duplicate_sequences:
        _issue(warnings, "retained_exact_duplicates", "Included sequences contain exact duplicates under different IDs",
               sequence_id_groups=duplicate_sequences)
    report["summary"].update(min_length=min(lengths, default=None), max_length=max(lengths, default=None),
                              unique_included_sequences=len(hash_groups))
    organisms = Counter(row.get("organism", "") for _, row in included
                        if not _pending(row.get("organism", "")))
    report["distributions"] = {
        "population": "included metadata records only; not independent biological observations",
        "length_aa": dict(sorted(Counter(lengths).items())),
        "length_unusable_records": len(included) - len(lengths),
        "organism": dict(sorted(organisms.items())),
        "organism_missing_records": len(included) - sum(organisms.values()),
    }
    report["source_completeness"] = {
        field: {"column_present": field in fields,
                "missing_included": sum(_pending(row.get(field, "")) for _, row in included),
                "missing_all": sum(_pending(row.get(field, "")) for row in rows)}
        for field in SOURCE_FIELDS
    }
    incomplete_sources = [field for field, item in report["source_completeness"].items()
                          if item["missing_included"] or not item["column_present"]]
    if incomplete_sources:
        _issue(warnings, "source_fields_incomplete", "Source fields need provenance review", fields=incomplete_sources)
    inferred_count = sum("inferred" in row.get("database", "").lower() for _, row in included)
    report["summary"]["inferred_database_records"] = inferred_count
    if inferred_count:
        _issue(warnings, "database_is_inferred", "Database naming is inferred, not independently verified", count=inferred_count)

    labelled = [(number, row) for number, row in included if not _pending(row.get("label", ""))]
    missing_evidence = []
    definitions, labels_by_sequence = set(), defaultdict(set)
    for number, row in labelled:
        missing = [field for field in EVIDENCE_FIELDS if _pending(row.get(field, ""))]
        if missing:
            missing_evidence.append(dict(sequence_id=row["sequence_id"], row=number, fields=missing))
        if not _pending(row.get("label_definition", "")):
            definitions.add(row["label_definition"])
        sequence = metadata_sequences.get(row["sequence_id"])
        if sequence is not None:
            labels_by_sequence[sequence].add(row["label"])
    label_conflicts = [
        {"sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(), "labels": sorted(labels)}
        for sequence, labels in labels_by_sequence.items() if len(labels) > 1
    ]
    if missing_evidence:
        _issue(warnings, "label_evidence_missing", "Assigned labels lack evidence fields", records=missing_evidence)
    if len(definitions) > 1:
        _issue(warnings, "label_definition_conflict", "Multiple label definitions require review; do not assume one target",
               definitions=sorted(definitions))
    if label_conflicts:
        _issue(warnings, "sequence_label_conflict", "Identical sequences have conflicting labels; biological context must be reviewed",
               conflicts=label_conflicts)
    label_counts = Counter(row["label"] for _, row in labelled)
    if not labelled:
        label_status = "unlabelled"
    elif missing_evidence or label_conflicts or len(definitions) > 1:
        label_status = "evidence_review_required"
    elif len(labelled) < len(included):
        label_status = "partially_labelled"
    else:
        label_status = "fields_present_requires_scientific_review"
    report["labels"] = {
        "status": label_status, "labelled_records": len(labelled),
        "unknown_records": len(included) - len(labelled), "value_counts": dict(sorted(label_counts.items())),
        "assigned_labels_missing_evidence": missing_evidence, "conflicts": label_conflicts,
        "supervised_training_approved": False,
        "reason": "A/B/C must confirm target, evidence and usable subset; D must confirm modelling and evaluation readiness.",
    }
    report["distributions"]["label"] = dict(sorted(label_counts.items()))
    report["distributions"]["label_unknown_records"] = len(included) - len(labelled)

    if split is None and any(row.get("split", "") or row.get("similarity_cluster", "") for _, row in included):
        split_rows = [row for _, row in included]
    report["split"] = {"status": "not_provided", "source": "file" if split else "metadata"}
    if split_rows is not None:
        split_start = len(errors)
        counts = Counter(row["sequence_id"] for row in split_rows)
        duplicate_ids = sorted(key for key, count in counts.items() if count > 1)
        if duplicate_ids:
            _issue(errors, "duplicate_split_ids", "Each sequence must occur in exactly one split row", sequence_ids=duplicate_ids)
        for code, identifiers in (
            ("split_missing_ids", set(metadata_counts) - set(counts)),
            ("split_extra_ids", set(counts) - set(metadata_counts)),
        ):
            if identifiers:
                _issue(errors, code, "Split must cover exactly the included metadata IDs", sequence_ids=sorted(identifiers))
        cluster_splits, sequence_splits = defaultdict(set), defaultdict(set)
        metadata_by_id = {row["sequence_id"]: row for _, row in included}
        for number, row in enumerate(split_rows, 2):
            identifier, group, name = row["sequence_id"], row.get("similarity_cluster", ""), row.get("split", "")
            if name not in ALLOWED_SPLITS:
                _issue(errors, "invalid_split_name", "Allowed split names are train, val and test", row=number, value=name)
            if not group:
                _issue(errors, "missing_similarity_cluster", "Every split row needs similarity_cluster", row=number)
            else:
                cluster_splits[group].add(name)
            if identifier in metadata_sequences:
                sequence_splits[metadata_sequences[identifier]].add(name)
            if split is not None and identifier in metadata_by_id:
                for field in ("similarity_cluster", "split"):
                    stored = metadata_by_id[identifier].get(field, "")
                    if stored and stored != row[field]:
                        _issue(errors, "split_metadata_conflict", "External split disagrees with populated metadata",
                               sequence_id=identifier, field=field)
        leaking_clusters = {group: sorted(names) for group, names in cluster_splits.items() if len(names) > 1}
        if leaking_clusters:
            _issue(errors, "cluster_cross_split", "Similarity clusters span multiple splits", clusters=leaking_clusters)
        leaking_sequences = {hashlib.sha256(sequence.encode("ascii")).hexdigest(): sorted(names)
                             for sequence, names in sequence_splits.items() if len(names) > 1}
        if leaking_sequences:
            _issue(errors, "exact_sequence_cross_split", "Identical sequences span multiple splits", sequences=leaking_sequences)
        report["split"].update(status="invalid" if len(errors) > split_start else "assignment_checks_passed",
                                 records=len(split_rows), clusters=len(cluster_splits),
                                 counts=dict(sorted(Counter(row.get("split", "") for row in split_rows).items())))
    elif split is not None:
        report["split"]["status"] = "invalid"
    report["format_valid"] = not errors
    report["descriptive_sequence_analysis_ready"] = not errors
    report["status"] = "invalid" if errors else "passed_with_warnings" if warnings else "passed"
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--split", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True, help="New directory; never overwritten")
    args = parser.parse_args()
    try:
        args.out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as error:
        parser.exit(2, f"error: {error}\n")
    report = audit_data(args.fasta, args.metadata, args.split)
    report["execution"] = {"python": platform.python_version(), "argv": sys.argv,
                           "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    try:
        (args.out_dir / "audit.json").write_text(
            json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    except OSError as error:
        parser.exit(2, f"error writing audit: {error}\n")
    print(f"{report['status']}: {len(report['errors'])} errors; labels={report['labels']['status']} -> {args.out_dir}")
    return 0 if report["format_valid"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
