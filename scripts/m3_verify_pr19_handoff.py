"""Verify the public PR #19 M2-C handoff without extracting its ZIP files.

This is an automated byte-level check, not a human review or M2 acceptance.
Only the pinned Release URL and the three frozen B inputs are accepted.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
import sys
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "reports" / "M2" / "c_scan_002_handoff.json"
DEFAULT_RECEIPT = ROOT / "reports" / "M3" / "pr19_handoff_verification.json"
RELEASE_URL = (
    "https://github.com/Lithium9767/MLP/releases/download/"
    "c-m2-scan-002-handoff/MLP_C_M2_scan_002_handoff.zip"
)
B_ZIP_NAME = "B_input/gvpa_v1_reproduction.zip"
MAX_DOWNLOAD_BYTES = 64 * 1024 * 1024
PINNED_ZIP_MANIFEST = {
    "bytes": 2218,
    "sha256": "483ed7a43f9d9f5244bc409ffe10173928cb4687808d4242c0d6498f3116d8b9",
}
FROZEN_B_INPUTS = {
    "gvpa_v1_reproduction/metadata.csv": {
        "bytes": 1041627,
        "sha256": "62893e31a2dbc245f29c0e9af0a282dcb5099d0e95b96db1a9cead2ac933004e",
    },
    "gvpa_v1_reproduction/sequences_for_clustering.fasta": {
        "bytes": 300688,
        "sha256": "4ff33c04be1d28f42ef5afa7d67261229d19c8e301784ea59750030196e8b122",
    },
    "gvpa_v1_reproduction/split/split_manifest.csv": {
        "bytes": 137073,
        "sha256": "0a6bd3b44d284a407e1cc9d5be758f0f9439c28395deb6d62d2c787f95341e9c",
    },
}


class VerificationError(ValueError):
    """The handoff cannot be verified against the pinned evidence."""


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def download_release(url: str, *, max_bytes: int = MAX_DOWNLOAD_BYTES) -> bytes:
    """Download the fixed public asset to memory with a strict size limit."""
    if url != RELEASE_URL:
        raise VerificationError(f"unexpected Release URL: {url}")
    request = urllib.request.Request(url, headers={"User-Agent": "MLP-M3-handoff-verifier/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        if response.status != 200:
            raise VerificationError(f"Release HTTP status is {response.status}, not 200")
        result = io.BytesIO()
        while chunk := response.read(1024 * 1024):
            if result.tell() + len(chunk) > max_bytes:
                raise VerificationError(f"Release exceeds {max_bytes} bytes")
            result.write(chunk)
    return result.getvalue()


def _expected_files(manifest: dict) -> dict[str, dict]:
    if not isinstance(manifest, dict):
        raise VerificationError("handoff manifest is not a JSON object")
    if manifest.get("shared_uri") != RELEASE_URL:
        raise VerificationError("manifest does not name the pinned Release URL")
    if manifest.get("archive_name") != RELEASE_URL.rsplit("/", 1)[-1]:
        raise VerificationError("manifest archive name does not match Release URL")
    if not isinstance(manifest.get("files"), list):
        raise VerificationError("manifest files is not a list")
    expected = {}
    for item in manifest["files"]:
        if not isinstance(item, dict) or not isinstance(item.get("name"), str):
            raise VerificationError("invalid manifest file entry")
        name = item["name"]
        if name in expected:
            raise VerificationError(f"duplicate manifest file: {name}")
        expected[name] = item
    if B_ZIP_NAME not in expected:
        raise VerificationError(f"manifest omits {B_ZIP_NAME}")
    return expected


def _verify_members(
    archive: zipfile.ZipFile,
    expected: dict[str, dict],
    *,
    exact_names: bool,
) -> list[dict]:
    entries = archive.infolist()
    names = [entry.filename for entry in entries]
    if len(names) != len(set(names)):
        raise VerificationError("ZIP has duplicate member names")
    missing = set(expected) - set(names)
    extra = set(names) - set(expected) if exact_names else set()
    if missing or extra:
        raise VerificationError(f"ZIP members differ: missing={sorted(missing)}, extra={sorted(extra)}")

    results = []
    for name, item in expected.items():
        info = archive.getinfo(name)
        expected_bytes = item.get("bytes")
        expected_hash = item.get("sha256")
        if not isinstance(expected_bytes, int) or isinstance(expected_bytes, bool) or expected_bytes < 0:
            raise VerificationError(f"invalid expected byte size for {name}")
        if not isinstance(expected_hash, str) or len(expected_hash) != 64:
            raise VerificationError(f"invalid expected SHA-256 for {name}")
        if info.file_size != expected_bytes:
            raise VerificationError(f"size mismatch for {name}: {info.file_size} != {expected_bytes}")
        digest = hashlib.sha256()
        count = 0
        with archive.open(info) as member:
            while chunk := member.read(1024 * 1024):
                count += len(chunk)
                if count > expected_bytes:
                    raise VerificationError(f"member exceeds expected size: {name}")
                digest.update(chunk)
        actual_hash = digest.hexdigest()
        if count != expected_bytes or actual_hash != expected_hash:
            raise VerificationError(f"content mismatch for {name}: {actual_hash} ({count} bytes)")
        results.append({"name": name, "bytes": count, "sha256": actual_hash, "status": "match"})
    return results


def verify_handoff(
    archive_bytes: bytes,
    manifest: dict,
    *,
    frozen_b_inputs: dict[str, dict] = FROZEN_B_INPUTS,
    internal_manifest_pin: dict = PINNED_ZIP_MANIFEST,
) -> dict:
    """Verify every outer member and the three frozen files within B's ZIP."""
    expected = _expected_files(manifest)
    expected["manifest.json"] = {"name": "manifest.json", **internal_manifest_pin}
    expected_bytes = manifest.get("archive_bytes")
    expected_hash = manifest.get("archive_sha256")
    if not isinstance(expected_bytes, int) or isinstance(expected_bytes, bool):
        raise VerificationError("invalid archive byte count in manifest")
    if not isinstance(expected_hash, str) or len(expected_hash) != 64:
        raise VerificationError("invalid archive SHA-256 in manifest")
    actual_hash = sha256(archive_bytes)
    if len(archive_bytes) != expected_bytes or actual_hash != expected_hash:
        raise VerificationError(
            f"Release archive mismatch: {actual_hash} ({len(archive_bytes)} bytes)"
        )
    try:
        with zipfile.ZipFile(io.BytesIO(archive_bytes)) as outer:
            outer_results = _verify_members(outer, expected, exact_names=True)
            internal_manifest = json.loads(outer.read("manifest.json"))
            if internal_manifest.get("files") != manifest["files"]:
                raise VerificationError("embedded manifest file list differs from external handoff manifest")
            if (internal_manifest.get("run_commit") != manifest.get("run_commit") or
                    internal_manifest.get("run_dirty") != manifest.get("run_dirty")):
                raise VerificationError("embedded run provenance differs from external handoff manifest")
            b_bytes = outer.read(B_ZIP_NAME)
        with zipfile.ZipFile(io.BytesIO(b_bytes)) as b_archive:
            b_results = _verify_members(b_archive, frozen_b_inputs, exact_names=False)
    except (zipfile.BadZipFile, OSError) as exc:
        raise VerificationError(f"unreadable handoff ZIP: {exc}") from exc
    return {
        "archive": {"bytes": len(archive_bytes), "sha256": actual_hash, "status": "match"},
        "outer_members": outer_results,
        "frozen_b_inputs": b_results,
    }


def make_receipt(
    manifest_bytes: bytes, archive_bytes: bytes, *, manifest_path: Path = DEFAULT_MANIFEST
) -> dict:
    manifest = json.loads(manifest_bytes)
    verified = verify_handoff(archive_bytes, manifest)
    resolved_manifest = manifest_path.resolve()
    try:
        manifest_label = resolved_manifest.relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        manifest_label = str(resolved_manifest)
    return {
        "schema_version": "1.0",
        "status": "verified",
        "verified_utc": datetime.now(timezone.utc).isoformat(),
        "reviewer_type": "assistant_automated",
        "release_url": RELEASE_URL,
        "manifest_path": manifest_label,
        "manifest_sha256": sha256(manifest_bytes),
        "implementation_sha256": sha256(Path(__file__).read_bytes()),
        "python_version": platform.python_version(),
        **verified,
        "scope_note": (
            "Automated public download and byte/hash verification only; "
            "not a non-author human review, PR approval, or M2/M3 acceptance."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--out", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()
    try:
        manifest_bytes = args.manifest.read_bytes()
        manifest = json.loads(manifest_bytes)
        if not isinstance(manifest, dict):
            raise VerificationError("handoff manifest is not a JSON object")
        archive_bytes = download_release(manifest.get("shared_uri"))
        receipt = make_receipt(manifest_bytes, archive_bytes, manifest_path=args.manifest)
    except (OSError, ValueError, KeyError, urllib.error.URLError) as exc:
        print(f"PR #19 handoff verification failed: {exc}", file=sys.stderr)
        return 1
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"verified {len(receipt['outer_members'])} Release members and 3 frozen B inputs; {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
