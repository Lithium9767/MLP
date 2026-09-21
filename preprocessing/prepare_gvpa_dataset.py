#!/usr/bin/env python3
"""Build an auditable GvpA candidate dataset from the recognition JSON in gv.zip."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


DEFAULT_MEMBER = "gv/recognition/data/gvp/gvpa/GvpA_sequences.json"
CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
PARTIAL_RE = re.compile(r"\bpartial\b", re.IGNORECASE)

METADATA_FIELDS = [
    "sequence_id", "accession", "sequence", "sequence_length", "sequence_sha256",
    "organism", "description", "source_file", "representative_gene",
    "representative_gvp_types", "member_count", "redundant_member_count", "member_gvp_types",
    "is_partial_representative", "is_partial_any_member", "has_type_conflict",
    "has_gvpj_member", "source_missing", "organism_missing", "nonstandard_residues",
    "is_exact_duplicate", "exact_duplicate_of", "is_length_outlier",
    "sequence_qc_eligible", "primary_analysis_eligible", "analysis_cohort", "exclusion_reason",
]

MEMBER_FIELDS = [
    "representative_id", "member_id", "organism", "description", "source_file",
    "gene", "gvp_types", "is_partial",
]


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_json_bytes(value: Any) -> bytes:
    text = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n"
    return text.encode("utf-8")


def annotation_text(annotation: dict[str, Any]) -> str:
    return " ".join(str(annotation.get(field, "")) for field in ("description", "full_header"))


def is_partial(annotation: dict[str, Any]) -> bool:
    return bool(PARTIAL_RE.search(annotation_text(annotation)))


def normalise_types(value: Any) -> list[str]:
    if value is None:
        return []
    values: Iterable[Any]
    if isinstance(value, (list, tuple, set)):
        values = value
    else:
        values = re.split(r"[,;|]", str(value))
    return sorted({str(item).strip() for item in values if str(item).strip()})


def require_dict(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be an object")
    return value


def load_source(zip_path: Path, member: str = DEFAULT_MEMBER) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not zip_path.is_file():
        raise FileNotFoundError(zip_path)
    with zipfile.ZipFile(zip_path) as archive:
        try:
            raw = archive.read(member)
        except KeyError as exc:
            matches = [name for name in archive.namelist() if name.lower().endswith("/gvpa_sequences.json")]
            detail = f"; possible members: {matches}" if matches else ""
            raise ValueError(f"Missing expected ZIP member {member}{detail}") from exc
    try:
        payload = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid JSON in ZIP member {member}: {exc}") from exc
    if not isinstance(payload, list) or not payload:
        raise ValueError("Recognition JSON must be a non-empty list")
    records = [require_dict(record, f"record {index}") for index, record in enumerate(payload, start=1)]
    provenance = {
        "source_zip": zip_path.name,
        "source_zip_sha256": sha256_file(zip_path),
        "source_member": member,
        "source_member_size": len(raw),
        "source_member_sha256": sha256_bytes(raw),
    }
    return records, provenance


def _bool(value: bool) -> str:
    return "true" if value else "false"


def build_rows(records: list[dict[str, Any]]) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, Any]]:
    rows: list[dict[str, str]] = []
    member_rows: list[dict[str, str]] = []
    first_by_hash: dict[str, str] = {}
    seen_ids: Counter[str] = Counter()

    for index, record in enumerate(records, start=1):
        required = (
            "unique_sequence_id", "sequence", "representative_annotation",
            "cluster_statistics", "redundant_members_details",
        )
        for field in required:
            if field not in record:
                raise ValueError(f"record {index} is missing required field {field}")
        sequence_id = str(record["unique_sequence_id"]).strip()
        sequence = re.sub(r"\s+", "", str(record["sequence"])).upper()
        representative = require_dict(record["representative_annotation"], f"record {index} representative_annotation")
        statistics_record = require_dict(record["cluster_statistics"], f"record {index} cluster_statistics")
        members_value = record["redundant_members_details"]
        if not isinstance(members_value, list):
            raise ValueError(f"record {index} redundant_members_details must be a list")
        members = [require_dict(member, f"record {index} member") for member in members_value]
        try:
            total_original = int(statistics_record["total_original_sequences"])
            redundant_count = int(statistics_record["redundant_count"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"record {index} has invalid cluster statistics") from exc
        if total_original != len(members) + 1 or redundant_count != len(members):
            raise ValueError(
                f"record {index} cluster statistics disagree with redundant_members_details"
            )

        digest = sha256_bytes(sequence.encode("utf-8"))
        duplicate_of = first_by_hash.get(digest, "")
        if not duplicate_of:
            first_by_hash[digest] = sequence_id or f"record_{index}"
        seen_ids[sequence_id] += 1

        rep_types = normalise_types(representative.get("gvp_types"))
        all_annotations = [representative, *members]
        all_types = sorted({item for member in all_annotations for item in normalise_types(member.get("gvp_types"))})
        non_gvpa_types = [item for item in all_types if item.casefold() != "gvpa"]
        partial_representative = is_partial(representative)
        partial_any = any(is_partial(member) for member in all_annotations)
        source_missing = not representative.get("source_file") or any(not member.get("source_file") for member in members)
        organism_missing = not representative.get("organism")
        invalid = "".join(sorted(set(sequence) - CANONICAL_AA))

        reasons: list[str] = []
        if not sequence_id:
            reasons.append("missing_sequence_id")
        if not sequence:
            reasons.append("empty_sequence")
        if invalid:
            reasons.append(f"nonstandard_residues:{invalid}")
        if duplicate_of:
            reasons.append(f"exact_duplicate_of:{duplicate_of}")
        if source_missing:
            reasons.append("source_missing")

        if invalid or not sequence_id or not sequence or duplicate_of:
            cohort = "excluded_sequence_qc"
        elif partial_any and non_gvpa_types:
            cohort = "sensitivity_partial_and_type_conflict"
        elif partial_any:
            cohort = "sensitivity_partial"
        elif non_gvpa_types:
            cohort = "sensitivity_type_conflict"
        else:
            cohort = "primary"
        eligible = not invalid and bool(sequence_id) and bool(sequence) and not duplicate_of

        rows.append({
            "sequence_id": sequence_id,
            "accession": str(representative.get("id", "")).strip(),
            "sequence": sequence,
            "sequence_length": str(len(sequence)),
            "sequence_sha256": digest,
            "organism": str(representative.get("organism", "")).strip(),
            "description": str(representative.get("description", "")).strip(),
            "source_file": str(representative.get("source_file", "")).strip(),
            "representative_gene": str(representative.get("gene", "")).strip(),
            "representative_gvp_types": "|".join(rep_types),
            "member_count": str(total_original),
            "redundant_member_count": str(redundant_count),
            "member_gvp_types": "|".join(all_types),
            "is_partial_representative": _bool(partial_representative),
            "is_partial_any_member": _bool(partial_any),
            "has_type_conflict": _bool(bool(non_gvpa_types)),
            "has_gvpj_member": _bool(any(item.casefold() == "gvpj" for item in all_types)),
            "source_missing": _bool(source_missing),
            "organism_missing": _bool(organism_missing),
            "nonstandard_residues": invalid,
            "is_exact_duplicate": _bool(bool(duplicate_of)),
            "exact_duplicate_of": duplicate_of,
            "is_length_outlier": "false",
            "sequence_qc_eligible": _bool(eligible),
            "primary_analysis_eligible": _bool(eligible and cohort == "primary"),
            "analysis_cohort": cohort,
            "exclusion_reason": ";".join(reasons),
        })

        for member in members:
            member_rows.append({
                "representative_id": sequence_id,
                "member_id": str(member.get("id", "")).strip(),
                "organism": str(member.get("organism", "")).strip(),
                "description": str(member.get("description", "")).strip(),
                "source_file": str(member.get("source_file", "")).strip(),
                "gene": str(member.get("gene", "")).strip(),
                "gvp_types": "|".join(normalise_types(member.get("gvp_types"))),
                "is_partial": _bool(is_partial(member)),
            })

    duplicate_ids = {item for item, count in seen_ids.items() if item and count > 1}
    for row in rows:
        if row["sequence_id"] in duplicate_ids:
            row["sequence_qc_eligible"] = "false"
            row["primary_analysis_eligible"] = "false"
            row["analysis_cohort"] = "excluded_sequence_qc"
            row["exclusion_reason"] = ";".join(filter(None, (row["exclusion_reason"], "duplicate_sequence_id")))

    eligible_lengths = [int(row["sequence_length"]) for row in rows if row["sequence_qc_eligible"] == "true"]
    mean_length = statistics.fmean(eligible_lengths) if eligible_lengths else math.nan
    std_length = statistics.pstdev(eligible_lengths) if len(eligible_lengths) > 1 else 0.0
    lower = mean_length - 3 * std_length
    upper = mean_length + 3 * std_length
    for row in rows:
        length = int(row["sequence_length"])
        row["is_length_outlier"] = _bool(bool(eligible_lengths) and (length < lower or length > upper))

    cohort_counts = Counter(row["analysis_cohort"] for row in rows)
    summary = {
        "candidate_records": len(rows),
        "total_original_sequence_records": sum(int(row["member_count"]) for row in rows),
        "redundant_member_records": len(member_rows),
        "unique_sequence_ids": len({row["sequence_id"] for row in rows}),
        "unique_sequence_hashes": len({row["sequence_sha256"] for row in rows}),
        "sequence_qc_eligible": sum(row["sequence_qc_eligible"] == "true" for row in rows),
        "primary_cohort": cohort_counts["primary"],
        "cohort_counts": dict(sorted(cohort_counts.items())),
        "partial_representative": sum(row["is_partial_representative"] == "true" for row in rows),
        "partial_any_member": sum(row["is_partial_any_member"] == "true" for row in rows),
        "type_conflict": sum(row["has_type_conflict"] == "true" for row in rows),
        "gvpj_member": sum(row["has_gvpj_member"] == "true" for row in rows),
        "source_missing": sum(row["source_missing"] == "true" for row in rows),
        "organism_missing": sum(row["organism_missing"] == "true" for row in rows),
        "nonstandard_sequence": sum(bool(row["nonstandard_residues"]) for row in rows),
        "exact_duplicate_sequence": sum(row["is_exact_duplicate"] == "true" for row in rows),
        "duplicate_sequence_id_records": sum(row["sequence_id"] in duplicate_ids for row in rows),
        "length": {
            "minimum": min(eligible_lengths) if eligible_lengths else None,
            "maximum": max(eligible_lengths) if eligible_lengths else None,
            "mean": mean_length if eligible_lengths else None,
            "population_std": std_length if eligible_lengths else None,
            "outlier_lower_exclusive": lower if eligible_lengths else None,
            "outlier_upper_exclusive": upper if eligible_lengths else None,
            "outlier_count": sum(row["is_length_outlier"] == "true" for row in rows),
        },
        "label_status": "unlabelled_no_function_classes",
        "homology_split_status": "pending_mmseqs2",
    }
    return rows, member_rows, summary


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(rows: list[dict[str, str]], member_rows: list[dict[str, str]], summary: dict[str, Any], provenance: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = output_dir / "metadata.csv"
    members_path = output_dir / "members.csv"
    fasta_path = output_dir / "sequences_for_clustering.fasta"
    write_csv(metadata_path, METADATA_FIELDS, rows)
    write_csv(members_path, MEMBER_FIELDS, member_rows)
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for row in rows:
            if row["sequence_qc_eligible"] != "true":
                continue
            handle.write(f">{row['sequence_id']}\n")
            sequence = row["sequence"]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    dataset_version = f"gvpa-recognition-{provenance['source_member_sha256'][:12]}"
    completed_summary = {
        "schema_version": "1.0", "dataset_version": dataset_version,
        "source": provenance, **summary,
        "output_sha256": {
            "metadata.csv": sha256_file(metadata_path),
            "members.csv": sha256_file(members_path),
            "sequences_for_clustering.fasta": sha256_file(fasta_path),
        },
    }
    summary_path = output_dir / "audit_summary.json"
    summary_path.write_bytes(stable_json_bytes(completed_summary))
    manifest = {
        "schema_version": "1.0",
        "dataset_version": dataset_version,
        "prediction_target": "candidate_functional_regions_and_residues",
        "function_labels_present": False,
        "sequence_qc_rule": "canonical amino acids; non-empty unique sequence ID; exact-sequence representative",
        "primary_analysis_rule": "sequence QC eligible; no partial annotation; no non-GvpA member type",
        "sensitivity_flags": ["partial", "non_GvpA member type", "length outlier", "missing source"],
        "homology_split_status": "pending_mmseqs2",
        "audit_summary_sha256": sha256_file(summary_path),
    }
    (output_dir / "dataset_manifest.json").write_bytes(stable_json_bytes(manifest))
    return completed_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-zip", type=Path, required=True)
    parser.add_argument("--member", default=DEFAULT_MEMBER)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    records, provenance = load_source(args.input_zip, args.member)
    rows, members, summary = build_rows(records)
    completed = write_outputs(rows, members, summary, provenance, args.output_dir)
    print(json.dumps(completed, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
