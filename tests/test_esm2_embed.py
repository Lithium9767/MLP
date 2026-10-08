import csv
import hashlib
import tempfile
import unittest
import subprocess
import sys
from pathlib import Path

from features.esm2_embed import iter_windows, normalise_sequence, read_embedding_records


class ESM2InputTest(unittest.TestCase):
    def test_split_and_primary_filter_are_enforced(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata.csv"
            sequence_a = "ACDEFG"
            sequence_b = "ACDEFA"
            fields = [
                "internal_id", "sequence_id", "sequence", "sequence_sha256",
                "sequence_qc_eligible", "primary_analysis_eligible", "analysis_cohort",
            ]
            with metadata.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows([
                    {
                        "internal_id": "GVPA_000001", "sequence_id": "raw|one", "sequence": sequence_a,
                        "sequence_sha256": hashlib.sha256(sequence_a.encode()).hexdigest(),
                        "sequence_qc_eligible": "true", "primary_analysis_eligible": "true",
                        "analysis_cohort": "primary",
                    },
                    {
                        "internal_id": "GVPA_000002", "sequence_id": "raw|two", "sequence": sequence_b,
                        "sequence_sha256": hashlib.sha256(sequence_b.encode()).hexdigest(),
                        "sequence_qc_eligible": "true", "primary_analysis_eligible": "false",
                        "analysis_cohort": "sensitivity_partial",
                    },
                ])
            split = root / "split.csv"
            with split.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["internal_id", "sequence_id", "split"])
                writer.writeheader()
                writer.writerows([
                    {"internal_id": "GVPA_000001", "sequence_id": "raw|one", "split": "discovery"},
                    {"internal_id": "GVPA_000002", "sequence_id": "raw|two", "split": "discovery"},
                ])
            primary = read_embedding_records(metadata, split)
            all_records = read_embedding_records(metadata, split, primary_only=False)
            self.assertEqual([record.internal_id for record in primary], ["GVPA_000001"])
            self.assertEqual(len(all_records), 2)

    def test_sequence_and_window_validation(self):
        self.assertEqual(normalise_sequence("ac de\nfg"), "ACDEFG")
        with self.assertRaisesRegex(ValueError, "Non-canonical"):
            normalise_sequence("ACDEX")
        self.assertEqual(iter_windows(40, 30, 5), [(1, 30), (6, 35), (11, 40)])
        self.assertEqual(iter_windows(20, 30, 5), [])

    def test_validation_cli_is_locked_before_reading_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run(
                [sys.executable, str(Path(__file__).resolve().parents[1] / "features/esm2_embed.py"),
                 "--metadata", "missing.csv", "--split-manifest", "missing.csv",
                 "--output-dir", directory, "--split", "validation"],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("Validation is locked", result.stderr)
            self.assertEqual(list(Path(directory).iterdir()), [])
