"""Small contracts for the explicit, frozen validation entry point."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np

from scripts import m3_spectral, m3_spectral_validate as validation


class FrozenValidationContractTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        root_patch = patch.object(validation, 'ROOT', self.root)
        root_patch.start()
        self.addCleanup(root_patch.stop)
        self.algorithm = self.root / 'feature_code.py'
        self.algorithm.write_text('# frozen mapping implementation\n', encoding='utf-8')
        self.config = {'feature_schema': 'synthetic:v1'}
        files = {}
        for name in ('config', 'model', 'metadata', 'split_manifest', 'coordinates',
                     'reference', 'discovery_features', 'discovery_windows'):
            files[name] = self.root / name
            files[name].write_text(json.dumps(self.config) if name == 'config' else name,
                                   encoding='utf-8')
        self.args = SimpleNamespace(**files, freeze=self.root / 'freeze.json',
                                    prepared=self.root / 'prepared', output=self.root / 'output')
        self.freeze = {
            'authorization': {'quote': 'Use the frozen validation population'},
            'evaluation': {'refit_PCA_or_clustering': False,
                           'permutation_seed': 42, 'within_sequence_label_permutations': 2},
            'config_sha256': m3_spectral.sha256_file(self.args.config),
            'spectral_model_sha256': m3_spectral.sha256_file(self.args.model),
            'implementation_sha256': {'feature_code.py': m3_spectral.sha256_file(self.algorithm)},
            'transfer_software': {'numpy': 'synthetic-test-version'},
            'inputs_sha256': {name: m3_spectral.sha256_file(files[name]) for name in
                             ('metadata', 'split_manifest', 'coordinates', 'reference',
                              'discovery_features', 'discovery_windows')},
            'expected_validation': {'sequences': 2, 'homology_groups': 2},
            'expected_discovery': {'sequences': 2, 'homology_groups': 2},
        }
        self.save_freeze()
        self.populations = {
            split: [{'internal_id': f'{split}_{i}',
                     'sequence_sha256': hashlib.sha256(f'{split}_{i}'.encode()).hexdigest(),
                     'homology_cluster': f'{split}_group_{i}'} for i in range(2)]
            for split in ('discovery', 'validation')
        }

    def save_freeze(self):
        self.args.freeze.write_text(json.dumps(self.freeze), encoding='utf-8')

    def test_protocol_accepts_exact_artifacts_without_loading_validation(self):
        with patch.object(validation.subprocess, 'check_output', return_value=''), \
                patch.object(validation, 'load_primary_records') as load_records:
            freeze, config = validation.checked_protocol(self.args)
        self.assertEqual(freeze, self.freeze)
        self.assertEqual(config, self.config)
        load_records.assert_not_called()

    def test_protocol_rejects_changed_model_config_or_biological_inputs(self):
        for name in ('config', 'model', 'metadata', 'split_manifest', 'coordinates', 'reference'):
            with self.subTest(changed=name):
                path = getattr(self.args, name)
                original = path.read_bytes()
                # Whitespace preserves valid JSON while changing its frozen bytes.
                path.write_bytes(original + b'\n')
                try:
                    with patch.object(validation.subprocess, 'check_output', return_value=''):
                        with self.assertRaisesRegex(ValueError, 'changed|mismatch'):
                            validation.checked_protocol(self.args)
                finally:
                    path.write_bytes(original)

    def test_protocol_requires_clean_source_explicit_scope_and_no_refit(self):
        with patch.object(validation.subprocess, 'check_output', return_value=' M script.py\n'):
            with self.assertRaisesRegex(ValueError, 'clean committed source'):
                validation.checked_protocol(self.args)
        original = copy.deepcopy(self.freeze)
        for missing_authorization in (True, False):
            with self.subTest(missing_authorization=missing_authorization):
                self.freeze = copy.deepcopy(original)
                if missing_authorization:
                    self.freeze['authorization']['quote'] = ''
                else:
                    self.freeze['evaluation']['refit_PCA_or_clustering'] = True
                self.save_freeze()
                with self.assertRaisesRegex(ValueError, 'Explicit frozen validation protocol'):
                    validation.checked_protocol(self.args)

    def test_protocol_rejects_committed_implementation_drift_or_missing_identity(self):
        self.algorithm.write_text('# changed mapping implementation\n', encoding='utf-8')
        with patch.object(validation.subprocess, 'check_output', return_value=''):
            with self.assertRaisesRegex(ValueError, 'Frozen implementation changed'):
                validation.checked_protocol(self.args)
        self.freeze.pop('implementation_sha256')
        self.save_freeze()
        with self.assertRaisesRegex(ValueError, 'Frozen implementation identities required'):
            validation.checked_protocol(self.args)

    def check_population(self, populations=None, freeze=None):
        populations = self.populations if populations is None else populations
        with patch.object(validation, 'load_primary_records',
                          side_effect=lambda metadata, manifest, split: populations[split]):
            return validation.checked_population(self.args, self.freeze if freeze is None else freeze)

    def test_population_accepts_separate_cohorts_and_checks_both_count_fields(self):
        records, counts = self.check_population()
        self.assertEqual(records, self.populations['validation'])
        self.assertEqual(counts['discovery'], {'sequences': 2, 'homology_groups': 2})
        for split in ('validation', 'discovery'):
            for field in ('sequences', 'homology_groups'):
                with self.subTest(split=split, field=field):
                    freeze = copy.deepcopy(self.freeze)
                    freeze['expected_' + split][field] += 1
                    with self.assertRaisesRegex(ValueError, f'Unexpected {split} population'):
                        self.check_population(freeze=freeze)

    def test_population_rejects_id_sequence_and_homology_leakage(self):
        for key in ('internal_id', 'sequence_sha256', 'homology_cluster'):
            with self.subTest(leaked=key):
                populations = copy.deepcopy(self.populations)
                populations['validation'][0][key] = populations['discovery'][0][key]
                with self.assertRaisesRegex(ValueError, f'leakage in {key}'):
                    self.check_population(populations=populations)

    def test_upstream_stage_requires_complete_receipt_same_freeze_and_unchanged_bytes(self):
        stage_dir = self.root / 'stage'
        stage_dir.mkdir()
        output = stage_dir / 'features.npz'
        output.write_bytes(b'synthetic bytes, never loaded as features')
        receipt = {'status': 'complete', 'stage': 'prepare', 'freeze_sha256': 'frozen-id',
                   'outputs_sha256': {'features.npz': m3_spectral.sha256_file(output)}}
        receipt_path = stage_dir / 'run_receipt.json'
        receipt_path.write_text(json.dumps(receipt), encoding='utf-8')
        self.assertEqual(validation.checked_stage(stage_dir, 'frozen-id', stage='prepare'), receipt)
        for field, value in [('status', 'started'), ('stage', 'embed'),
                             ('freeze_sha256', 'different-freeze')]:
            with self.subTest(field=field):
                altered = {**receipt, field: value}
                receipt_path.write_text(json.dumps(altered), encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'not completed under this freeze'):
                    validation.checked_stage(stage_dir, 'frozen-id', stage='prepare')
        receipt_path.write_text(json.dumps(receipt), encoding='utf-8')
        output.write_bytes(b'tampered output')
        with self.assertRaisesRegex(ValueError, 'Upstream output changed'):
            validation.checked_stage(stage_dir, 'frozen-id', stage='prepare')

    def test_ordinary_bundle_loader_keeps_validation_locked_until_explicit_scope(self):
        feature_path, window_path = self.root / 'features.npz', self.root / 'windows.csv'
        row = {'internal_id': 'v1', 'sequence_sha256': 'a' * 64, 'homology_cluster': 'vg1',
               'split': 'validation', 'sequence_length': '10', 'raw_start': '1', 'raw_end': '10',
               'window_length': '10', 'accepted_mapping_count': '8'}
        m3_spectral.write_csv(window_path, [row])
        np.savez(feature_path, features=np.ones((1, 2)), internal_id=np.array(['v1']),
                 raw_start=np.array([1]), raw_end=np.array([10]))
        config = {**self.config, 'window_length': 10, 'window_step': 5,
                  'feature_dimension': 2, 'minimum_accepted_mapping': 8}
        source = {'feature_schema': self.config['feature_schema'],
                  'features_sha256': m3_spectral.sha256_file(feature_path),
                  'windows_sha256': m3_spectral.sha256_file(window_path)}
        for fitting in (True, False):
            with self.subTest(fitting=fitting):
                with self.assertRaisesRegex(ValueError, 'Validation remains locked'):
                    m3_spectral.load_bundle(feature_path, window_path, source, config, fitting=fitting)
        _, rows = m3_spectral.load_bundle(feature_path, window_path, source, config,
                                         fitting=False, prediction_splits=frozenset({'validation'}))
        self.assertEqual(rows, [row])

    def test_transfer_grants_only_validation_scope_and_does_not_refit_model(self):
        features = np.ones((2, 2))
        rows = self.populations['validation']
        model = SimpleNamespace(diagnostics={'training_cohort_sha256': 'ordered-discovery'},
                                labels=np.array([0, 1]), predict=Mock(return_value={'cluster': np.array([0, 1])}))
        with patch.object(validation, 'checked_stage'), \
                patch.object(validation.importlib.metadata, 'version', return_value='synthetic-test-version'), \
                patch.object(validation, 'read_json', return_value={}), \
                patch.object(validation, 'load_bundle', side_effect=[(features, rows),
                                   (features, self.populations['discovery'])]) as load_bundle, \
                patch.object(validation, 'checked_model', return_value=model), \
                patch.object(validation, 'cohort_identity', return_value='ordered-discovery'), \
                patch.object(validation, 'predict_command'), \
                patch('models.spectral_validation.summarize_validation', return_value=({}, [], [])), \
                patch('models.spectral_regions.fit_spectral') as fit, \
                patch.object(validation, 'write_json'), patch.object(validation, 'write_csv'):
            result = validation.transfer_stage(self.args, self.freeze, self.config, rows)
        self.assertEqual(load_bundle.call_args_list[0].kwargs,
                         {'fitting': False, 'prediction_splits': frozenset({'validation'})})
        self.assertEqual(load_bundle.call_args_list[1].kwargs, {'fitting': True})
        fit.assert_not_called()
        self.assertFalse(result['fit_performed'])
        self.assertTrue(result['frozen_model_unchanged'])


if __name__ == '__main__':
    unittest.main()
