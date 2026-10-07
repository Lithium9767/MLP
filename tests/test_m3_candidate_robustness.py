import unittest

from scripts.m3_candidate_robustness import (
    match_same_sequence_position,
    sequence_level_delta,
    shorter_length_auc,
    validate_discovery_rows,
)


def window(identifier, start, position, cluster, has_state, length=120):
    return {
        "internal_id": identifier,
        "raw_start": str(start),
        "raw_end": str(start + 29),
        "sequence_length": str(length),
        "window_length": "30",
        "split": "discovery",
        "cluster": str(cluster),
        "normalized_start": str(position),
        "hmm_states": "25" if has_state else "24",
        "mapped_fraction": "0.8",
    }


class CandidateRobustnessTest(unittest.TestCase):
    def test_position_matching_never_borrows_another_sequence(self):
        rows = [window("a", 1, .20, 3, True),
                window("a", 6, .23, 2, False),
                window("a", 11, .90, 2, False),
                window("b", 1, .20, 2, False)]
        pairs, unmatched = match_same_sequence_position(rows, 3, tolerance=.05)
        self.assertEqual(len(pairs), 1)
        self.assertEqual(unmatched, 0)
        self.assertEqual(pairs[0][1]["internal_id"], "a")
        self.assertEqual(pairs[0][1]["raw_start"], "6")

    def test_unmatched_target_is_reported_not_replaced_by_distant_window(self):
        rows = [window("a", 1, .20, 3, True), window("a", 6, .90, 2, False)]
        pairs, unmatched = match_same_sequence_position(rows, 3, tolerance=.05)
        self.assertEqual(pairs, [])
        self.assertEqual(unmatched, 1)

    def test_sequence_aggregation_prevents_window_pseudoreplication(self):
        pairs = [(window("a", 1, .20, 3, True), window("a", 6, .22, 2, False)),
                 (window("a", 11, .30, 3, True), window("a", 16, .32, 2, False)),
                 (window("b", 1, .20, 3, False), window("b", 6, .22, 2, True))]
        self.assertEqual(sequence_level_delta(pairs, 25), {"a": 1.0, "b": -1.0})

    def test_validation_input_is_rejected(self):
        rows = [window("a", 1, .2, 3, True)]
        rows[0]["split"] = "validation"
        with self.assertRaisesRegex(ValueError, "discovery"):
            validate_discovery_rows(rows)

    def test_length_only_auc_handles_ties_and_direction(self):
        self.assertEqual(shorter_length_auc([80, 90], [110, 120]), 1.0)
        self.assertEqual(shorter_length_auc([100], [100]), 0.5)
        self.assertEqual(shorter_length_auc([120], [80]), 0.0)


if __name__ == "__main__":
    unittest.main()
