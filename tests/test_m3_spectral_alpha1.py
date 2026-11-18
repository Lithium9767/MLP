import unittest

import numpy as np

from scripts.m3_spectral_alpha1 import matched_rates


class MatchedRatesTests(unittest.TestCase):
    def test_unshared_cells_do_not_affect_comparison(self):
        rows = [{'internal_id': str(i), 'homology_cluster': str(i)} for i in range(3)]
        result, support = matched_rates(rows, [True, False, True], [True, False, False],
                                       ['same', 'same', 'candidate only'], np.arange(3))
        self.assertEqual(result['difference'], 1.)
        self.assertEqual(result['comparable_mass'], .5)
        self.assertEqual(support, [{'0'}, {'1'}])

    def test_no_estimate_when_site_is_missing_in_sibling(self):
        rows = [{'internal_id': str(i), 'homology_cluster': str(i)} for i in range(2)]
        result, _ = matched_rates(rows, [True, False], [True, None], ['a', 'a'], np.arange(2))
        self.assertIsNone(result['difference'])
        self.assertEqual(result['groups_candidate'], 0)
