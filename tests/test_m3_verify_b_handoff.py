"""Offline adversarial checks for B's nested M3 input ZIP verifier."""

from __future__ import annotations

import csv
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
import m3_verify_b_handoff as verifier


def csv_bytes(columns, rows):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode()


def zip_bytes(items):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_STORED) as archive:
        for name, content in items:
            archive.writestr(name, content)
    return buffer.getvalue()


class BZipVerificationTest(unittest.TestCase):
    def setUp(self):
        self.records = [
            ("GVPA_A", "A.1", "AAAA", "primary", "discovery", "CL_A"),
            ("GVPA_B", "B.1", "BBBB", "primary", "validation", "CL_B"),
            ("GVPA_C", "C.1", "CCCC", "sensitivity_partial", "discovery", "CL_C"),
        ]
        metadata = []
        split = []
        manifests = {name: [] for name in verifier.MANIFESTS}
        for internal_id, sequence_id, sequence, cohort, bucket, cluster in self.records:
            digest = hashlib.sha256(sequence.encode()).hexdigest()
            base = {
                "internal_id": internal_id, "sequence_id": sequence_id,
                "sequence_sha256": digest, "sequence_length": str(len(sequence)),
                "analysis_cohort": cohort,
            }
            metadata.append({**base, "primary_analysis_eligible": str(cohort == "primary").lower(),
                             "sequence": sequence})
            split.append({"internal_id": internal_id, "sequence_id": sequence_id,
                          "analysis_cohort": cohort,
                          "primary_analysis_eligible": str(cohort == "primary").lower(),
                          "homology_cluster": cluster, "split": bucket})
            name = f"{bucket}_primary_manifest.csv" if cohort == "primary" else "sensitivity_manifest.csv"
            manifests[name].append({key: {**base, "split": bucket, "homology_cluster": cluster}[key]
                                    for key in verifier.MANIFEST_COLUMNS})
        self.files = {
            "metadata.csv": csv_bytes(tuple(metadata[0]), metadata),
            "split_manifest.csv": csv_bytes(
                ("internal_id", "sequence_id", "analysis_cohort", "primary_analysis_eligible",
                 "homology_cluster", "split"), split),
            "sequences_for_clustering.fasta": b">GVPA_A\nAAAA\n>GVPA_B\nBBBB\n>GVPA_C\nCCCC\n",
            **{name: csv_bytes(verifier.MANIFEST_COLUMNS, rows) for name, rows in manifests.items()},
        }
        self.expected_counts = dict.fromkeys(verifier.MANIFESTS, 1)

    def summaries(self, files=None):
        files = self.files if files is None else files
        digest = lambda data: hashlib.sha256(data).hexdigest()
        return (
            {"dataset_version": "test-dataset", "output_sha256": {
                "metadata.csv": digest(files["metadata.csv"]),
                "sequences_for_clustering.fasta": digest(files["sequences_for_clustering.fasta"]),
            }},
            {"split_version": "test-split", "metadata_sha256": digest(files["metadata.csv"]),
             "split_manifest_sha256": digest(files["split_manifest.csv"]),
             "cohort_counts_by_split": {
                 "discovery": {"primary": 1, "sensitivity_partial": 1},
                 "validation": {"primary": 1},
             }, "primary_sequence_counts": {"discovery": 1, "validation": 1}},
        )

    def archive(self, files=None, *, outer_manifest=None, outer_extra=(), inner_extra=()):
        files = self.files if files is None else files
        nested = zip_bytes([*files.items(), *inner_extra])
        return zip_bytes([
            (verifier.PREFIX, b""),
            *((verifier.PREFIX + name, files[name] if outer_manifest is None else outer_manifest[name])
              for name in verifier.MANIFESTS),
            (verifier.PREFIX + "m3_b_handoff_bundle.zip", nested),
            *outer_extra,
        ])

    def verify(self, files=None, *, archive=None, summaries=None):
        files = self.files if files is None else files
        audit, split = self.summaries(files) if summaries is None else summaries
        return verifier.verify_handoff(
            self.archive(files) if archive is None else archive,
            audit_summary=audit, split_summary=split,
            expected_counts=self.expected_counts,
        )

    def test_valid_nested_handoff_receipt_has_no_sequences_or_paths(self):
        receipt = self.verify()
        self.assertEqual(receipt["status"], "verified_automated_input_integrity")
        self.assertEqual(receipt["split_records"], 3)
        self.assertEqual(receipt["homology_clusters"], 3)
        self.assertEqual([receipt["manifests"][name]["rows"] for name in verifier.MANIFESTS], [1, 1, 1])
        serialized = json.dumps(receipt)
        self.assertNotIn("AAAA", serialized)
        self.assertNotIn("GVPA_A", serialized)
        self.assertNotIn("C:\\", serialized)

    def test_outer_and_inner_manifest_must_match_byte_for_byte(self):
        changed = dict(self.files)
        changed["discovery_primary_manifest.csv"] += b"\n"
        outer = {name: self.files[name] for name in verifier.MANIFESTS}
        with self.assertRaisesRegex(verifier.VerificationError, "outer/inner manifest differs"):
            self.verify(changed, archive=self.archive(changed, outer_manifest=outer))

    def test_frozen_metadata_hash_cannot_be_replaced_by_zip(self):
        changed = dict(self.files)
        changed["metadata.csv"] += b"\n"
        with self.assertRaisesRegex(verifier.VerificationError, "frozen SHA-256 mismatch"):
            self.verify(changed, summaries=self.summaries())

    def test_unknown_manifest_id_is_rejected(self):
        changed = dict(self.files)
        changed["discovery_primary_manifest.csv"] = changed["discovery_primary_manifest.csv"].replace(
            b"GVPA_A", b"GVPA_X")
        with self.assertRaisesRegex(verifier.VerificationError, "unknown manifest internal_id"):
            self.verify(changed)

    def test_wrong_sequence_digest_is_rejected_even_when_repinned(self):
        changed = dict(self.files)
        changed["metadata.csv"] = changed["metadata.csv"].replace(
            hashlib.sha256(b"AAAA").hexdigest().encode(), b"0" * 64)
        with self.assertRaisesRegex(verifier.VerificationError, "metadata sequence"):
            self.verify(changed)

    def test_cross_split_homology_cluster_is_rejected(self):
        changed = dict(self.files)
        changed["split_manifest.csv"] = changed["split_manifest.csv"].replace(b"CL_B", b"CL_A")
        with self.assertRaisesRegex(verifier.VerificationError, "homology cluster crosses"):
            self.verify(changed)

    def test_duplicate_manifest_id_is_rejected(self):
        changed = dict(self.files)
        changed["discovery_primary_manifest.csv"] += changed["discovery_primary_manifest.csv"].splitlines(keepends=True)[1]
        with self.assertRaisesRegex(verifier.VerificationError, "duplicate internal_id"):
            self.verify(changed)

    def test_unexpected_or_traversal_member_is_rejected(self):
        for extra in (("extra.txt", b"x"), ("../escape", b"x")):
            with self.subTest(extra=extra[0]):
                with self.assertRaises(verifier.VerificationError):
                    self.verify(archive=self.archive(outer_extra=(extra,)))

    def test_nested_extra_member_is_rejected(self):
        with self.assertRaisesRegex(verifier.VerificationError, "inner ZIP members differ"):
            self.verify(archive=self.archive(inner_extra=(("run.py", b"print('no')"),)))

    def test_duplicate_zip_member_is_rejected(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            duplicate = self.archive(outer_extra=((verifier.PREFIX + "discovery_primary_manifest.csv", b"x"),))
        with self.assertRaisesRegex(verifier.VerificationError, "duplicate outer ZIP member"):
            self.verify(archive=duplicate)

    def test_corrupted_zip_payload_is_rejected_by_crc(self):
        archive = self.archive()
        at = archive.index(b">GVPA_A\nAAAA")
        corrupted = archive[:at] + b"?" + archive[at + 1:]
        with self.assertRaisesRegex(verifier.VerificationError, "unreadable or unsupported handoff ZIP"):
            self.verify(archive=corrupted)


if __name__ == "__main__":
    unittest.main()
