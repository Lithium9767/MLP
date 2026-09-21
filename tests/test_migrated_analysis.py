import csv
import hashlib
import tempfile
import unittest
from pathlib import Path

from bioinformatics.hmm_coordinates import (
    conservation_by_state,
    map_window_to_hmm_states,
    parse_alignment_columns,
)
from bioinformatics.structure_mapping import (
    build_structure_hmm_rows,
    parse_pdb_chain,
    structure_sequence,
)
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


class HMMCoordinateTest(unittest.TestCase):
    def test_insertions_and_deletions_preserve_coordinates(self):
        rows = parse_alignment_columns("AB.CD", "AXY-D", hmm_from=5, target_from=10)
        self.assertEqual(
            [(row["raw_position"], row["hmm_match_state"], row["is_insertion"]) for row in rows],
            [(10, 5, False), (11, 6, False), (12, None, True), (13, 8, False)],
        )
        self.assertEqual(map_window_to_hmm_states(rows, 11, 13), [6, 8])

    def test_conservation_uses_total_sequence_occupancy(self):
        rows = [
            {"hmm_match_state": 1, "residue": "A"},
            {"hmm_match_state": 1, "residue": "A"},
            {"hmm_match_state": 2, "residue": "C"},
        ]
        table = conservation_by_state(rows, total_sequences=4)
        self.assertEqual(table[0]["consensus"], "A")
        self.assertEqual(table[0]["occupancy"], 0.5)
        self.assertEqual(table[1]["occupancy"], 0.25)


class StructureMappingTest(unittest.TestCase):
    @staticmethod
    def atom(serial, residue, chain, number):
        return (
            f"ATOM  {serial:5d}  CA  {residue:>3s} {chain}{number:4d}    "
            "   0.000   0.000   0.000  1.00 20.00           C"
        )

    def test_pdb_chain_and_hmm_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.pdb"
            path.write_text(
                self.atom(1, "ALA", "N", 10) + "\n" +
                self.atom(2, "GLY", "N", 11) + "\n" +
                self.atom(3, "VAL", "N", 12) + "\n",
                encoding="ascii",
            )
            residues = parse_pdb_chain(path, "N")
            self.assertEqual(structure_sequence(residues), "AGV")
            rows = build_structure_hmm_rows(
                residues,
                [
                    {"raw_position": 1, "hmm_match_state": 3, "is_insertion": False},
                    {"raw_position": 2, "hmm_match_state": 4, "is_insertion": False},
                ],
            )
            self.assertEqual(rows[0]["pdb_residue_number"], 10)
            self.assertEqual(rows[1]["hmm_match_state"], 4)


if __name__ == "__main__":
    unittest.main()
