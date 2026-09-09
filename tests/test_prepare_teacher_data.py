import json
import tempfile
import unittest
from pathlib import Path

from preprocessing.prepare_teacher_data import audit, load_records, write_outputs


class PrepareTeacherDataTest(unittest.TestCase):
    def test_load_and_audit_teacher_json(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "gvpa"
            source.mkdir()
            payload = [
                {
                    "unique_sequence_id": "WP_000001.1",
                    "sequence": "ACD EFG",
                    "representative_annotation": {"organism": "Example species"},
                }
            ]
            (source / "GvpA_sequences.json").write_text(json.dumps(payload))

            rows = load_records(root)
            summary = audit(rows)

            self.assertEqual(rows[0]["label"], "GvpA")
            self.assertEqual(rows[0]["sequence"], "ACDEFG")
            self.assertEqual(rows[0]["valid_residues"], "true")
            self.assertEqual(summary["total_records"], 1)
            self.assertEqual(summary["cross_label_exact_sequence_conflicts"], 0)

    def test_write_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input" / "gvpc"
            source.mkdir(parents=True)
            (source / "records.json").write_text(
                json.dumps([{"unique_sequence_id": "P12345", "sequence": "ACDX"}])
            )
            rows = load_records(root / "input")

            summary = write_outputs(rows, root / "output")

            self.assertEqual(summary["total_records"], 1)
            self.assertTrue((root / "output" / "metadata.csv").is_file())
            self.assertTrue((root / "output" / "sequences.fasta").is_file())
            self.assertTrue((root / "output" / "audit_summary.json").is_file())


if __name__ == "__main__":
    unittest.main()
