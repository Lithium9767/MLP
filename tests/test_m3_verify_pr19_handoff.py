"""Offline tests for the pinned PR #19 byte-level handoff verifier."""

import hashlib
import io
import json
import sys
import unittest
import warnings
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import m3_verify_pr19_handoff as verifier


def zip_bytes(members):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, content in members.items():
            archive.writestr(name, content)
    return buffer.getvalue()


class HandoffVerificationTest(unittest.TestCase):
    def setUp(self):
        self.b_files = {
            "gvpa_v1_reproduction/metadata.csv": b"metadata\n",
            "gvpa_v1_reproduction/sequences_for_clustering.fasta": b">id\nAA\n",
            "gvpa_v1_reproduction/split/split_manifest.csv": b"split\n",
        }
        self.expected_b = {
            name: {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
            for name, content in self.b_files.items()
        }
        self.outer_files = {
            verifier.B_ZIP_NAME: zip_bytes(self.b_files),
            "C_scan_002/raw/failures.json": b'{"failures": []}\n',
        }

    def fixture(self, *, outer_files=None):
        listed_files = self.outer_files if outer_files is None else outer_files
        file_list = [
            {"name": name, "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
            for name, content in listed_files.items()
        ]
        internal = json.dumps(
            {"files": file_list, "run_commit": "fixture", "run_dirty": False},
            sort_keys=True,
        ).encode()
        self.internal_manifest_bytes = internal
        files = {**listed_files, "manifest.json": internal}
        archive = zip_bytes(files)
        manifest = {
            "shared_uri": verifier.RELEASE_URL,
            "archive_name": verifier.RELEASE_URL.rsplit("/", 1)[-1],
            "archive_bytes": len(archive),
            "archive_sha256": hashlib.sha256(archive).hexdigest(),
            "files": file_list,
            "run_commit": "fixture",
            "run_dirty": False,
        }
        pin = {"bytes": len(internal), "sha256": hashlib.sha256(internal).hexdigest()}
        return archive, manifest, pin

    def test_valid_outer_and_nested_inputs(self):
        archive, manifest, pin = self.fixture()
        result = verifier.verify_handoff(
            archive, manifest, frozen_b_inputs=self.expected_b, internal_manifest_pin=pin
        )
        self.assertEqual(result["archive"]["status"], "match")
        self.assertEqual(len(result["outer_members"]), 3)
        self.assertEqual(len(result["frozen_b_inputs"]), 3)
        self.assertTrue(all(row["status"] == "match" for row in result["frozen_b_inputs"]))

    def test_wrong_release_bytes_fail_before_zip_inspection(self):
        archive, manifest, pin = self.fixture()
        with self.assertRaisesRegex(verifier.VerificationError, "Release archive mismatch"):
            verifier.verify_handoff(
                archive + b"x", manifest, frozen_b_inputs=self.expected_b, internal_manifest_pin=pin
            )

    def test_unlisted_outer_member_is_rejected(self):
        _, manifest, pin = self.fixture()
        expanded = {**self.outer_files, "manifest.json": self.internal_manifest_bytes, "unlisted.txt": b"unexpected"}
        expanded_archive = zip_bytes(expanded)
        manifest["archive_bytes"] = len(expanded_archive)
        manifest["archive_sha256"] = hashlib.sha256(expanded_archive).hexdigest()
        with self.assertRaisesRegex(verifier.VerificationError, "extra=\\['unlisted.txt'\\]"):
            verifier.verify_handoff(
                expanded_archive, manifest, frozen_b_inputs=self.expected_b, internal_manifest_pin=pin
            )

    def test_missing_listed_outer_member_is_rejected(self):
        _, manifest, pin = self.fixture()
        missing_archive = zip_bytes({
            verifier.B_ZIP_NAME: self.outer_files[verifier.B_ZIP_NAME],
            "manifest.json": self.internal_manifest_bytes,
        })
        manifest["archive_bytes"] = len(missing_archive)
        manifest["archive_sha256"] = hashlib.sha256(missing_archive).hexdigest()
        with self.assertRaisesRegex(verifier.VerificationError, "missing="):
            verifier.verify_handoff(
                missing_archive, manifest, frozen_b_inputs=self.expected_b, internal_manifest_pin=pin
            )

    def test_outer_member_hash_is_checked(self):
        archive, manifest, pin = self.fixture()
        manifest["files"][1]["sha256"] = "0" * 64
        with self.assertRaisesRegex(verifier.VerificationError, "content mismatch"):
            verifier.verify_handoff(
                archive, manifest, frozen_b_inputs=self.expected_b, internal_manifest_pin=pin
            )

    def test_nested_b_hash_is_frozen_not_trusted_from_outer_manifest(self):
        changed_b = dict(self.b_files)
        changed_b["gvpa_v1_reproduction/metadata.csv"] = b"changed metadata\n"
        changed_outer = dict(self.outer_files)
        changed_outer[verifier.B_ZIP_NAME] = zip_bytes(changed_b)
        archive, manifest, pin = self.fixture(outer_files=changed_outer)
        with self.assertRaisesRegex(verifier.VerificationError, "metadata.csv"):
            verifier.verify_handoff(
                archive, manifest, frozen_b_inputs=self.expected_b, internal_manifest_pin=pin
            )

    def test_duplicate_zip_member_is_rejected(self):
        _, manifest, pin = self.fixture()
        stream = io.BytesIO()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            with zipfile.ZipFile(stream, "w") as target:
                target.writestr(verifier.B_ZIP_NAME, self.outer_files[verifier.B_ZIP_NAME])
                target.writestr(verifier.B_ZIP_NAME, self.outer_files[verifier.B_ZIP_NAME])
                target.writestr("C_scan_002/raw/failures.json", self.outer_files["C_scan_002/raw/failures.json"])
                target.writestr("manifest.json", self.internal_manifest_bytes)
        duplicate = stream.getvalue()
        manifest["archive_bytes"] = len(duplicate)
        manifest["archive_sha256"] = hashlib.sha256(duplicate).hexdigest()
        with self.assertRaisesRegex(verifier.VerificationError, "duplicate member"):
            verifier.verify_handoff(
                duplicate, manifest, frozen_b_inputs=self.expected_b, internal_manifest_pin=pin
            )

    def test_release_url_cannot_be_silently_changed(self):
        archive, manifest, pin = self.fixture()
        manifest["shared_uri"] = "https://example.org/handoff.zip"
        with self.assertRaisesRegex(verifier.VerificationError, "pinned Release URL"):
            verifier.verify_handoff(
                archive, manifest, frozen_b_inputs=self.expected_b, internal_manifest_pin=pin
            )
        with self.assertRaises(verifier.VerificationError):
            verifier.download_release("https://example.org/handoff.zip")

    def test_repository_manifest_names_all_twelve_entries(self):
        manifest = json.loads(verifier.DEFAULT_MANIFEST.read_text(encoding="utf-8"))
        self.assertEqual(len(verifier._expected_files(manifest)), 12)


if __name__ == "__main__":
    unittest.main()
