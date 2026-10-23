"""Contract tests for audited bundles and discovery-only transfer evaluation."""
from __future__ import annotations

import copy
import csv
import hashlib
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.model_selection import GroupShuffleSplit

from scripts import m3_spectral


class BundleContractTests(unittest.TestCase):
    """Tiny sources exercise identity boundaries without model fitting."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.features_path = self.root / 'features.npz'
        self.windows_path = self.root / 'windows.csv'
        self.features = np.arange(18, dtype=np.float32).reshape(6, 3)
        self.rows = []
        for sequence in ('a', 'b', 'c'):
            for start in (1, 6):
                self.rows.append({
                    'internal_id': sequence,
                    'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
                    'homology_cluster': f'group_{sequence}',
                    'split': 'discovery', 'sequence_length': '20',
                    'raw_start': str(start), 'raw_end': str(start + 9),
                    'window_length': '10', 'accepted_mapping_count': '8',
                })
        self.config = {
            'feature_schema': 'synthetic_test_schema:v1',
            'feature_dimension': 3, 'window_length': 10, 'window_step': 5,
            'minimum_accepted_mapping': 8,
            'discovery_expectations': {'windows': 6, 'sequences': 3, 'groups': 3},
        }
        self.write_fixture()

    def write_fixture(self, *, rows=None, arrays=None, fieldnames=None):
        rows = copy.deepcopy(self.rows if rows is None else rows)
        fields = list(self.rows[0]) if fieldnames is None else fieldnames
        with self.windows_path.open('w', encoding='utf-8-sig', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
        payload = {
            'features': self.features,
            'internal_id': np.array([row['internal_id'] for row in rows]),
            'raw_start': np.array([int(row['raw_start']) for row in rows]),
            'raw_end': np.array([int(row['raw_end']) for row in rows]),
        }
        if arrays:
            payload.update(arrays)
        np.savez_compressed(self.features_path, **payload)
        self.source = {
            'feature_schema': self.config['feature_schema'],
            'features_sha256': m3_spectral.sha256_file(self.features_path),
            'windows_sha256': m3_spectral.sha256_file(self.windows_path),
        }

    def load(self, *, fitting=True, config=None):
        return m3_spectral.load_bundle(
            self.features_path, self.windows_path, self.source,
            self.config if config is None else config, fitting=fitting,
        )

    def test_verified_bundle_preserves_feature_and_window_order(self):
        features, rows = self.load()
        np.testing.assert_array_equal(features, self.features)
        self.assertEqual(rows, self.rows)

    def test_byte_changes_fail_even_when_contents_remain_readable(self):
        for target in (self.features_path, self.windows_path):
            with self.subTest(file=target.name):
                self.write_fixture()
                with target.open('ab') as handle:
                    handle.write(b'\n')
                with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                    self.load()

    def test_schema_mismatch_is_rejected(self):
        self.source['feature_schema'] = 'a_different_preprocessing_definition'
        with self.assertRaisesRegex(ValueError, 'schema mismatch'):
            self.load()

    def test_reordered_npz_keys_cannot_relabel_the_csv_rows(self):
        self.write_fixture(arrays={
            'raw_start': np.array([6, 1, 1, 6, 1, 6]),
            'raw_end': np.array([15, 10, 10, 15, 10, 15]),
        })
        with self.assertRaisesRegex(ValueError, 'reordered'):
            self.load()

    def test_duplicate_csv_window_keys_are_rejected(self):
        rows = copy.deepcopy(self.rows)
        rows[1] = copy.deepcopy(rows[0])
        self.write_fixture(rows=rows)
        with self.assertRaisesRegex(ValueError, 'Duplicate window key'):
            self.load()

    def test_required_column_cannot_be_omitted(self):
        self.write_fixture(fieldnames=[name for name in self.rows[0]
                                      if name != 'accepted_mapping_count'])
        with self.assertRaisesRegex(ValueError, 'Window table requires'):
            self.load()

    def test_empty_window_population_is_rejected(self):
        self.write_fixture(rows=[])
        with self.assertRaisesRegex(ValueError, 'No windows'):
            self.load()

    def test_validation_is_rejected_for_fitting_and_prediction(self):
        rows = copy.deepcopy(self.rows)
        rows[-1]['split'] = 'validation'
        self.write_fixture(rows=rows)
        for fitting in (True, False):
            with self.subTest(fitting=fitting):
                with self.assertRaisesRegex(ValueError, 'Validation remains locked'):
                    self.load(fitting=fitting)

    def test_new_split_is_predictable_but_cannot_train(self):
        rows = copy.deepcopy(self.rows)
        for row in rows:
            row['split'] = 'new'
        self.write_fixture(rows=rows)
        with self.assertRaisesRegex(ValueError, 'Validation remains locked'):
            self.load(fitting=True)
        config = copy.deepcopy(self.config)
        config['discovery_expectations'] = {'windows': 999, 'sequences': 999, 'groups': 999}
        _, loaded = self.load(fitting=False, config=config)
        self.assertEqual({row['split'] for row in loaded}, {'new'})

    def test_each_discovery_population_count_is_checked(self):
        for field in ('windows', 'sequences', 'groups'):
            with self.subTest(field=field):
                config = copy.deepcopy(self.config)
                config['discovery_expectations'][field] += 1
                with self.assertRaisesRegex(ValueError, 'Discovery population mismatch'):
                    self.load(config=config)

    def test_invalid_identity_fields_are_rejected(self):
        for field, value in (
            ('internal_id', ''), ('homology_cluster', ''),
            ('sequence_sha256', 'z' * 64), ('sequence_sha256', 'a' * 63),
        ):
            with self.subTest(field=field, value=value):
                rows = copy.deepcopy(self.rows)
                rows[0][field] = value
                self.write_fixture(rows=rows)
                with self.assertRaisesRegex(ValueError, 'sequence identity'):
                    self.load()

    def test_one_sequence_cannot_have_conflicting_metadata(self):
        for field, value in (
            ('sequence_sha256', 'f' * 64), ('homology_cluster', 'different_group'),
            ('sequence_length', '25'), ('split', 'new'),
        ):
            with self.subTest(field=field):
                rows = copy.deepcopy(self.rows)
                rows[1][field] = value
                self.write_fixture(rows=rows)
                with self.assertRaisesRegex(ValueError, 'Conflicting sequence metadata'):
                    self.load(fitting=False)

    def test_identical_sequences_cannot_cross_homology_groups(self):
        rows = copy.deepcopy(self.rows)
        for row in rows[2:4]:
            row['sequence_sha256'] = rows[0]['sequence_sha256']
        self.write_fixture(rows=rows)
        with self.assertRaisesRegex(ValueError, 'crosses homology groups'):
            self.load()

    def test_window_bounds_step_and_mapping_count_are_checked(self):
        changes = (
            {'raw_start': '0', 'raw_end': '9'},
            {'raw_start': '2', 'raw_end': '11'},
            {'raw_start': '16', 'raw_end': '25'},
            {'raw_end': '11'}, {'sequence_length': '9'},
            {'window_length': '20'}, {'accepted_mapping_count': '7'},
            {'accepted_mapping_count': '11'},
        )
        for change in changes:
            with self.subTest(change=change):
                rows = copy.deepcopy(self.rows)
                rows[0].update(change)
                self.write_fixture(rows=rows)
                with self.assertRaisesRegex(ValueError, 'bounds/representation'):
                    self.load()

    def test_csv_fractional_position_and_mapping_count_are_rejected(self):
        # Rewrite the CSV after saving valid integer NPZ keys.
        for field in ('raw_start', 'raw_end', 'accepted_mapping_count'):
            with self.subTest(field=field):
                self.write_fixture()
                rows = copy.deepcopy(self.rows)
                rows[0][field] = '8.5'
                with self.windows_path.open('w', encoding='utf-8', newline='') as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                    writer.writeheader()
                    writer.writerows(rows)
                self.source['windows_sha256'] = m3_spectral.sha256_file(self.windows_path)
                with self.assertRaises(ValueError):
                    self.load()

    def test_npz_key_arrays_must_have_exactly_one_row_per_window(self):
        valid = {
            'internal_id': np.array([row['internal_id'] for row in self.rows]),
            'raw_start': np.array([int(row['raw_start']) for row in self.rows]),
            'raw_end': np.array([int(row['raw_end']) for row in self.rows]),
        }
        for name, values in valid.items():
            for invalid in (values[:-1], np.append(values, values[-1]), values[:, None]):
                with self.subTest(name=name, shape=invalid.shape):
                    self.write_fixture(arrays={name: invalid})
                    with self.assertRaisesRegex(ValueError, 'equal-length 1D'):
                        self.load()

    def test_numeric_window_key_arrays_cannot_be_float_or_bool(self):
        for name in ('raw_start', 'raw_end'):
            values = np.array([int(row[name]) for row in self.rows])
            for invalid in (values.astype(float), values.astype(bool), values.astype(str)):
                with self.subTest(name=name, dtype=str(invalid.dtype)):
                    self.write_fixture(arrays={name: invalid})
                    with self.assertRaisesRegex(ValueError, 'strings/integers'):
                        self.load()

    def test_feature_matrix_shape_numeric_type_and_finiteness_are_checked(self):
        nan_features = self.features.copy()
        nan_features[0, 0] = np.nan
        inf_features = self.features.copy()
        inf_features[0, 0] = np.inf
        for invalid in (
            self.features[:-1], self.features[:, :2], self.features.reshape(-1),
            self.features.astype(int), self.features.astype(str), nan_features, inf_features,
        ):
            with self.subTest(shape=invalid.shape, dtype=str(invalid.dtype)):
                self.write_fixture(arrays={'features': invalid})
                with self.assertRaisesRegex(ValueError, 'Invalid feature shape/type/values'):
                    self.load()


class GroupHoldoutContractTests(unittest.TestCase):
    """Mock model fitting; retain real splits and label-matching calculations."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.output = Path(temporary.name)
        self.features = np.column_stack((np.arange(12), np.arange(12) + .1,
                                         np.arange(12) + .2)).astype(float)
        self.rows = [
            {'internal_id': f'sequence_{index // 2}', 'homology_cluster': f'group_{index // 2}',
             'sequence_sha256': hashlib.sha256(f'sequence_{index // 2}'.encode()).hexdigest(),
             'split': 'discovery',
             'raw_start': str(1 + 5 * (index % 2)), 'raw_end': str(10 + 5 * (index % 2))}
            for index in range(12)
        ]
        self.config = {
            'feature_schema': 'synthetic_test_schema:v1',
            'model': {'n_clusters': 2, 'pca_components': 3, 'graph_neighbors': 3,
                      'transfer_neighbors': 3, 'seed': 42},
            'evaluation': {'heldout_group_fraction': 1 / 3, 'seeds': [3, 4]},
        }
        self.reference = SimpleNamespace(
            config=m3_spectral.SpectralConfig(**self.config['model']),
            pca_mean=np.zeros(3), components=np.eye(3),
            training_scores=self.features.copy(), labels=np.arange(12) % 2,
            diagnostics={'training_cohort_sha256': m3_spectral.cohort_identity(self.rows)},
        )
        self.fitted_indices = []
        self.predicted_indices = []

    def fake_fit(self, features, config):
        indices = features[:, 0].astype(int)
        self.fitted_indices.append(indices.copy())

        def predict(query):
            indices = query[:, 0].astype(int)
            self.predicted_indices.append(indices.copy())
            return {'cluster': 1 - indices % 2, 'vote_fraction': np.ones(len(indices)),
                    'nearest_distance': np.full(len(indices), .1)}

        return SimpleNamespace(labels=1 - indices % 2, predict=predict, diagnostics={})

    def evaluate(self):
        return m3_spectral.evaluate_command(
            self.features, self.rows, self.config, self.output / 'reference.npz', self.output,
        )

    def test_group_holdout_never_fits_heldout_windows_and_matches_using_train_only(self):
        with patch.object(m3_spectral, 'load_model', return_value=self.reference), \
             patch.object(m3_spectral, 'fit_spectral', side_effect=self.fake_fit), \
             patch.object(m3_spectral, 'save_model'), \
             patch.object(m3_spectral, 'linear_sum_assignment', wraps=linear_sum_assignment) as matching:
            summary = self.evaluate()
        groups = np.array([row['homology_cluster'] for row in self.rows])
        for run, seed in enumerate(self.config['evaluation']['seeds']):
            train, heldout = next(GroupShuffleSplit(
                n_splits=1, test_size=1 / 3, random_state=seed,
            ).split(self.features, groups=groups))
            np.testing.assert_array_equal(self.fitted_indices[run], train)
            np.testing.assert_array_equal(self.predicted_indices[run], heldout)
            self.assertFalse(set(groups[train]) & set(groups[heldout]))
            counts = -matching.call_args_list[run].args[0]
            self.assertEqual(int(counts.sum()), len(train))
            expected = np.zeros((2, 2), dtype=int)
            np.add.at(expected, (1 - train % 2, self.reference.labels[train]), 1)
            np.testing.assert_array_equal(counts, expected)
            self.assertEqual(summary['runs'][run]['heldout_window_agreement'], 1.)
        self.assertTrue(summary['reference_includes_heldout_groups'])
        self.assertEqual(summary['label_matching_uses'], 'training rows only')
        self.assertFalse(summary['independent_function_accuracy'])
        self.assertFalse(summary['validation_used'])
        self.assertTrue((self.output / 'group_holdout_predictions.csv').is_file())

    def test_heldout_reference_labels_cannot_influence_hungarian_mapping(self):
        train, heldout = np.arange(8), np.arange(8, 12)
        self.reference.labels[heldout] = 1 - self.reference.labels[heldout]
        self.config['evaluation']['seeds'] = [42]
        splitter = Mock()
        splitter.split.return_value = iter([(train, heldout)])
        with patch.object(m3_spectral, 'load_model', return_value=self.reference), \
             patch.object(m3_spectral, 'GroupShuffleSplit', return_value=splitter), \
             patch.object(m3_spectral, 'fit_spectral', side_effect=self.fake_fit), \
             patch.object(m3_spectral, 'save_model'), \
             patch.object(m3_spectral, 'linear_sum_assignment', wraps=linear_sum_assignment) as matching:
            summary = self.evaluate()
        np.testing.assert_array_equal(-matching.call_args.args[0], [[0, 4], [4, 0]])
        self.assertEqual(summary['runs'][0]['heldout_window_agreement'], 0.)

    def test_reference_from_another_order_fails_before_fitting(self):
        self.reference.training_scores = self.reference.training_scores[::-1].copy()
        with patch.object(m3_spectral, 'load_model', return_value=self.reference), \
             patch.object(m3_spectral, 'fit_spectral') as fitting:
            with self.assertRaisesRegex(ValueError, 'ordered feature cohort'):
                self.evaluate()
        fitting.assert_not_called()

    def test_model_parameters_must_match_for_prediction_and_evaluation(self):
        self.config['model']['transfer_neighbors'] = 4
        for operation in ('predict', 'evaluate'):
            with self.subTest(operation=operation):
                with patch.object(m3_spectral, 'load_model', return_value=self.reference), \
                     patch.object(m3_spectral, 'fit_spectral') as fitting:
                    with self.assertRaisesRegex(ValueError, 'Model parameters do not match'):
                        if operation == 'predict':
                            m3_spectral.predict_command(
                                self.features, self.rows, self.config,
                                self.output / 'reference.npz', self.output,
                            )
                        else:
                            self.evaluate()
                fitting.assert_not_called()

    def test_reference_window_identity_and_group_changes_fail_before_fitting(self):
        originals = copy.deepcopy(self.rows)
        for field, value in (
            ('internal_id', 'replacement_sequence'), ('sequence_sha256', 'f' * 64),
            ('homology_cluster', 'replacement_group'), ('split', 'new'),
            ('raw_start', '6'), ('raw_end', '15'),
        ):
            with self.subTest(field=field):
                self.rows = copy.deepcopy(originals)
                self.rows[0][field] = value
                with patch.object(m3_spectral, 'load_model', return_value=self.reference), \
                     patch.object(m3_spectral, 'fit_spectral') as fitting:
                    with self.assertRaisesRegex(ValueError, 'window identities/groups do not match'):
                        self.evaluate()
                fitting.assert_not_called()

    def test_fit_binds_ordered_window_identity_into_saved_model(self):
        with patch.object(m3_spectral, 'fit_spectral', side_effect=self.fake_fit), \
             patch.object(m3_spectral, 'save_model') as saving:
            summary = m3_spectral.fit_command(
                self.features, self.rows, self.config, self.output,
            )
        saved_model = saving.call_args.args[0]
        self.assertEqual(saved_model.diagnostics['training_cohort_sha256'],
                         m3_spectral.cohort_identity(self.rows))
        self.assertEqual(summary['model_diagnostics']['training_cohort_sha256'],
                         saved_model.diagnostics['training_cohort_sha256'])
        self.assertFalse(summary['validation_used'])

    def test_leaking_group_split_fails_before_any_fitting(self):
        splitter = Mock()
        splitter.split.return_value = iter([(np.arange(0, 12, 2), np.arange(1, 12, 2))])
        with patch.object(m3_spectral, 'load_model', return_value=self.reference), \
             patch.object(m3_spectral, 'GroupShuffleSplit', return_value=splitter), \
             patch.object(m3_spectral, 'fit_spectral') as fitting:
            with self.assertRaisesRegex(ValueError, 'Homology group leakage'):
                self.evaluate()
        fitting.assert_not_called()


if __name__ == '__main__':
    unittest.main()
