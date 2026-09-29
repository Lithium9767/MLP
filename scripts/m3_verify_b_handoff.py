"""Verify B's M3 data handoff ZIP against the frozen M2 release, offline.

The archive is treated as untrusted data: nothing is extracted or executed.
This check is an automated input receipt, not approval of a biological result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import stat
import sys
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT_SUMMARY = ROOT / "results" / "data_audit" / "gvpa_v1_audit_summary.json"
SPLIT_SUMMARY = ROOT / "results" / "data_audit" / "gvpa_v1_split_summary.json"
PREFIX = "m3_b_handoff/"
MANIFESTS = (
    "discovery_primary_manifest.csv",
    "validation_primary_manifest.csv",
    "sensitivity_manifest.csv",
)
INNER_FILES = frozenset((*MANIFESTS, "metadata.csv", "sequences_for_clustering.fasta", "split_manifest.csv"))
OUTER_FILES = frozenset((PREFIX, *(PREFIX + name for name in MANIFESTS), PREFIX + "m3_b_handoff_bundle.zip"))
MANIFEST_COLUMNS = (
    "internal_id", "sequence_id", "sequence_sha256", "sequence_length",
    "analysis_cohort", "split", "homology_cluster",
)
MAX_ARCHIVE_BYTES = 32 * 1024 * 1024
MAX_MEMBER_BYTES = 16 * 1024 * 1024
EXPECTED_COUNTS = {
    "discovery_primary_manifest.csv": 1202,
    "validation_primary_manifest.csv": 519,
    "sensitivity_manifest.csv": 355,
}
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")


class VerificationError(ValueError):
    """The handoff fails a byte-level or semantic integrity check."""


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _check_hash(actual: bytes, expected: str, label: str) -> None:
    if not isinstance(expected, str) or not SHA256_RE.fullmatch(expected):
        raise VerificationError(f"invalid frozen SHA-256 for {label}")
    if _sha(actual) != expected:
        raise VerificationError(f"frozen SHA-256 mismatch: {label}")


def _check_names(archive: zipfile.ZipFile, expected: frozenset[str], label: str) -> None:
    infos = archive.infolist()
    names = [info.filename for info in infos]
    if len(names) != len(set(names)):
        raise VerificationError(f"duplicate {label} ZIP member")
    if len(infos) > 32:
        raise VerificationError(f"too many {label} ZIP members")
    for info in infos:
        name = info.filename
        if (name.startswith("/") or "\\" in name or ":" in name or
                any(part in {"", ".", ".."} for part in name.rstrip("/").split("/"))):
            raise VerificationError(f"unsafe {label} ZIP path")
        if stat.S_ISLNK(info.external_attr >> 16) or info.flag_bits & 1:
            raise VerificationError(f"symlink or encrypted {label} ZIP member")
        if info.file_size > MAX_MEMBER_BYTES or info.file_size < 0:
            raise VerificationError(f"oversized {label} ZIP member")
    if set(names) != expected:
        raise VerificationError(
            f"{label} ZIP members differ: missing={sorted(expected - set(names))}, "
            f"extra={sorted(set(names) - expected)}"
        )
    if sum(info.file_size for info in infos) > MAX_MEMBER_BYTES:
        raise VerificationError(f"{label} ZIP uncompressed total exceeds limit")


def _read(archive: zipfile.ZipFile, name: str) -> bytes:
    info = archive.getinfo(name)
    if info.file_size > MAX_MEMBER_BYTES:
        raise VerificationError(f"oversized ZIP member: {name}")
    with archive.open(info) as source:
        data = source.read(MAX_MEMBER_BYTES + 1)
    if len(data) != info.file_size:
        raise VerificationError(f"ZIP member size mismatch: {name}")
    return data


def _csv_rows(data: bytes, name: str, *, columns: tuple[str, ...] | None = None,
              required: frozenset[str] = frozenset()) -> list[dict[str, str]]:
    try:
        text = data.decode("utf-8-sig")
        reader = csv.DictReader(io.StringIO(text, newline=""))
        header = reader.fieldnames
        if header is None or len(header) != len(set(header)):
            raise VerificationError(f"missing or duplicate CSV header: {name}")
        if columns is not None and tuple(header) != columns:
            raise VerificationError(f"unexpected CSV columns: {name}")
        if not required.issubset(header):
            raise VerificationError(f"missing required CSV columns: {name}")
        rows = list(reader)
    except (UnicodeError, csv.Error) as exc:
        raise VerificationError(f"unreadable CSV: {name}") from exc
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise VerificationError(f"malformed CSV row: {name}")
    return rows


def _unique_rows(rows: list[dict[str, str]], name: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        internal_id = row.get("internal_id", "")
        if not internal_id or internal_id in result:
            raise VerificationError(f"blank or duplicate internal_id: {name}")
        result[internal_id] = row
    return result


def _fasta_records(data: bytes) -> dict[str, str]:
    try:
        text = data.decode("ascii")
    except UnicodeError as exc:
        raise VerificationError("non-ASCII FASTA") from exc
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if current is not None:
                if not chunks:
                    raise VerificationError("empty FASTA sequence")
                records[current] = "".join(chunks)
            current = line[1:]
            if not current or any(char.isspace() for char in current) or current in records:
                raise VerificationError("blank or duplicate FASTA ID")
            chunks = []
        elif current is None or not line or not re.fullmatch(r"[A-Z]+", line):
            raise VerificationError("malformed FASTA sequence")
        else:
            chunks.append(line)
    if current is not None:
        if not chunks:
            raise VerificationError("empty FASTA sequence")
        records[current] = "".join(chunks)
    if not records:
        raise VerificationError("empty FASTA")
    return records


def verify_handoff(
    archive_bytes: bytes,
    *,
    audit_summary: dict | None = None,
    split_summary: dict | None = None,
    expected_counts: dict[str, int] | None = None,
) -> dict:
    """Return a sequence-free receipt only after all frozen and semantic checks pass."""
    if len(archive_bytes) > MAX_ARCHIVE_BYTES:
        raise VerificationError("outer ZIP exceeds size limit")
    if audit_summary is None:
        audit_summary = json.loads(AUDIT_SUMMARY.read_text(encoding="utf-8"))
    if split_summary is None:
        split_summary = json.loads(SPLIT_SUMMARY.read_text(encoding="utf-8"))
    if expected_counts is None:
        expected_counts = EXPECTED_COUNTS
    if set(expected_counts) != set(MANIFESTS) or any(
        not isinstance(value, int) or isinstance(value, bool) or value < 0
        for value in expected_counts.values()
    ):
        raise VerificationError("invalid expected manifest counts")
    try:
        with zipfile.ZipFile(io.BytesIO(archive_bytes)) as outer:
            _check_names(outer, OUTER_FILES, "outer")
            if not outer.getinfo(PREFIX).is_dir():
                raise VerificationError("outer directory marker is not a directory")
            outer_manifests = {name: _read(outer, PREFIX + name) for name in MANIFESTS}
            nested_bytes = _read(outer, PREFIX + "m3_b_handoff_bundle.zip")
        with zipfile.ZipFile(io.BytesIO(nested_bytes)) as inner:
            _check_names(inner, INNER_FILES, "inner")
            inner_files = {name: _read(inner, name) for name in INNER_FILES}
    except (zipfile.BadZipFile, RuntimeError, OSError, EOFError, NotImplementedError) as exc:
        raise VerificationError("unreadable or unsupported handoff ZIP") from exc

    for name in MANIFESTS:
        if outer_manifests[name] != inner_files[name]:
            raise VerificationError(f"outer/inner manifest differs: {name}")
    audit_hashes = audit_summary.get("output_sha256", {})
    frozen_hashes = {
        "metadata.csv": audit_hashes.get("metadata.csv"),
        "sequences_for_clustering.fasta": audit_hashes.get("sequences_for_clustering.fasta"),
        "split_manifest.csv": split_summary.get("split_manifest_sha256"),
    }
    if split_summary.get("metadata_sha256") != frozen_hashes["metadata.csv"]:
        raise VerificationError("frozen M2 summaries disagree on metadata SHA-256")
    for name, digest in frozen_hashes.items():
        _check_hash(inner_files[name], digest, name)

    metadata_rows = _csv_rows(
        inner_files["metadata.csv"], "metadata.csv",
        required=frozenset((*MANIFEST_COLUMNS[:4], "analysis_cohort", "primary_analysis_eligible", "sequence")),
    )
    metadata = _unique_rows(metadata_rows, "metadata.csv")
    sequence_ids = [row["sequence_id"] for row in metadata_rows]
    if any(not value for value in sequence_ids) or len(sequence_ids) != len(set(sequence_ids)):
        raise VerificationError("blank or duplicate metadata sequence_id")
    split_rows = _csv_rows(
        inner_files["split_manifest.csv"], "split_manifest.csv",
        columns=("internal_id", "sequence_id", "analysis_cohort", "primary_analysis_eligible",
                 "homology_cluster", "split"),
    )
    split = _unique_rows(split_rows, "split_manifest.csv")
    fasta = _fasta_records(inner_files["sequences_for_clustering.fasta"])
    if set(split) != set(fasta):
        raise VerificationError("FASTA IDs differ from frozen split IDs")
    if not set(split).issubset(metadata):
        raise VerificationError("frozen split contains unknown metadata ID")
    if len(split) != sum(expected_counts.values()):
        raise VerificationError("frozen split total differs from handoff counts")

    cluster_split: dict[str, str] = {}
    cohort_counts: dict[str, Counter[str]] = {"discovery": Counter(), "validation": Counter()}
    for internal_id, row in split.items():
        meta = metadata[internal_id]
        if row["split"] not in cohort_counts or not row["homology_cluster"]:
            raise VerificationError("invalid frozen split or homology cluster")
        cluster = row["homology_cluster"]
        if cluster in cluster_split and cluster_split[cluster] != row["split"]:
            raise VerificationError("homology cluster crosses discovery/validation")
        cluster_split[cluster] = row["split"]
        for field in ("sequence_id", "analysis_cohort", "primary_analysis_eligible"):
            if row[field] != meta[field]:
                raise VerificationError(f"frozen split/metadata mismatch: {field}")
        if row["primary_analysis_eligible"] != ("true" if row["analysis_cohort"] == "primary" else "false"):
            raise VerificationError("primary eligibility/cohort mismatch")
        sequence = fasta[internal_id]
        if sequence != meta["sequence"] or len(sequence) != int(meta["sequence_length"]):
            raise VerificationError("FASTA/metadata sequence or length mismatch")
        _check_hash(sequence.encode("ascii"), meta["sequence_sha256"], "metadata sequence")
        cohort_counts[row["split"]][row["analysis_cohort"]] += 1
    if split_summary.get("cohort_counts_by_split") != {
        name: dict(counts) for name, counts in cohort_counts.items()
    }:
        raise VerificationError("frozen M2 cohort counts differ from split data")

    seen: set[str] = set()
    manifest_hashes: dict[str, dict] = {}
    for name in MANIFESTS:
        rows = _csv_rows(inner_files[name], name, columns=MANIFEST_COLUMNS)
        index = _unique_rows(rows, name)
        if len(rows) != expected_counts[name]:
            raise VerificationError(f"wrong manifest count: {name}")
        for internal_id, row in index.items():
            if internal_id in seen or internal_id not in split:
                raise VerificationError("duplicate or unknown manifest internal_id")
            seen.add(internal_id)
            frozen_split = split[internal_id]
            meta = metadata[internal_id]
            expected_name = (
                f"{frozen_split['split']}_primary_manifest.csv"
                if frozen_split["analysis_cohort"] == "primary"
                else "sensitivity_manifest.csv"
            )
            if name != expected_name:
                raise VerificationError("manifest record is in wrong cohort or split file")
            expected_row = {
                "internal_id": internal_id,
                "sequence_id": meta["sequence_id"],
                "sequence_sha256": meta["sequence_sha256"],
                "sequence_length": meta["sequence_length"],
                "analysis_cohort": frozen_split["analysis_cohort"],
                "split": frozen_split["split"],
                "homology_cluster": frozen_split["homology_cluster"],
            }
            if row != expected_row:
                raise VerificationError("manifest record differs from frozen M2 data")
        manifest_hashes[name] = {"sha256": _sha(inner_files[name]), "rows": len(rows)}
    if seen != set(split):
        raise VerificationError("handoff manifests omit frozen M2 split IDs")
    if split_summary.get("primary_sequence_counts") != {
        name: cohort_counts[name]["primary"] for name in ("discovery", "validation")
    }:
        raise VerificationError("frozen M2 primary counts differ from split data")

    return {
        "schema_version": "1.0",
        "status": "verified_automated_input_integrity",
        "verified_utc": datetime.now(timezone.utc).isoformat(),
        "archive": {"sha256": _sha(archive_bytes), "bytes": len(archive_bytes)},
        "nested_archive": {"sha256": _sha(nested_bytes), "bytes": len(nested_bytes)},
        "frozen_dataset_version": audit_summary.get("dataset_version"),
        "frozen_split_version": split_summary.get("split_version"),
        "frozen_input_sha256": frozen_hashes,
        "manifests": manifest_hashes,
        "split_records": len(split),
        "homology_clusters": len(cluster_split),
        "cluster_leakage": False,
        "scope_note": (
            "Automated local ZIP check only; it does not establish a fixed shared location, "
            "independent human handoff acceptance, or an M3 experiment result."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zip", dest="zip_path", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.zip_path.stat().st_size > MAX_ARCHIVE_BYTES:
            raise VerificationError("outer ZIP exceeds size limit")
        audit_bytes = AUDIT_SUMMARY.read_bytes()
        split_bytes = SPLIT_SUMMARY.read_bytes()
        receipt = verify_handoff(
            args.zip_path.read_bytes(),
            audit_summary=json.loads(audit_bytes),
            split_summary=json.loads(split_bytes),
        )
        receipt["verifier_sha256"] = _sha(Path(__file__).read_bytes())
        receipt["frozen_summary_sha256"] = {
            "gvpa_v1_audit_summary.json": _sha(audit_bytes),
            "gvpa_v1_split_summary.json": _sha(split_bytes),
        }
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"M3 B handoff verification failed: {exc}", file=sys.stderr)
        return 1
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("M3 B handoff input integrity verified; receipt written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
