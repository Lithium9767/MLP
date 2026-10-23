#!/usr/bin/env python3
"""Fit spectral regions and transfer frozen labels on prepared window features.

This boundary consumes the audited upstream ESM window bundle. It does not
extract embeddings, infer structural coordinates, or unlock validation data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score
from sklearn.model_selection import GroupShuffleSplit
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from models.spectral_regions import SpectralConfig, fit_spectral, save_model, load_model


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8'))


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError('Cannot export an empty table')
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_bundle(features_path: Path, windows_path: Path, source: dict, config: dict,
                *, fitting: bool) -> tuple[np.ndarray, list[dict]]:
    """Verify bytes, per-row keys, cohort and representation before model use.

    A source manifest attests file identity, not biological validity. The
    producer remains responsible for the stated upstream feature definition.
    """
    if source['feature_schema'] != config['feature_schema']:
        raise ValueError('Feature schema mismatch')
    for path, key in [(features_path, 'features_sha256'), (windows_path, 'windows_sha256')]:
        if sha256_file(path) != source[key]:
            raise ValueError(f'Input hash mismatch: {path.name}')
    with windows_path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        required = {'internal_id', 'sequence_sha256', 'homology_cluster', 'split',
                    'raw_start', 'raw_end', 'sequence_length', 'window_length', 'accepted_mapping_count'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'Window table requires {sorted(required)}')
        rows = list(reader)
    if not rows:
        raise ValueError('No windows')
    keys = [(r['internal_id'], int(r['raw_start']), int(r['raw_end'])) for r in rows]
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicate window key')
    sequence_info, digest_group = {}, {}
    for row, (_, start, end) in zip(rows, keys):
        if row['split'] not in ({'discovery'} if fitting else {'discovery', 'new'}):
            raise ValueError('Validation remains locked; allowed splits are discovery/new')
        digest = row['sequence_sha256']
        if (not row['internal_id'] or not row['homology_cluster'] or len(digest) != 64
                or any(c not in '0123456789abcdef' for c in digest)):
            raise ValueError('Missing sequence identity or homology group')
        if (start < 1 or (start-1) % config['window_step'] != 0
                or end > int(row['sequence_length']) or end-start+1 != config['window_length']
                or int(row['window_length']) != config['window_length']
                or not config['minimum_accepted_mapping'] <= int(row['accepted_mapping_count']) <= config['window_length']):
            raise ValueError('Window bounds/representation selection mismatch')
        info = (digest, row['homology_cluster'], row['split'], int(row['sequence_length']))
        if sequence_info.setdefault(row['internal_id'], info) != info:
            raise ValueError('Conflicting sequence metadata')
        if digest_group.setdefault(digest, row['homology_cluster']) != row['homology_cluster']:
            raise ValueError('Identical sequence hash crosses homology groups')
    with np.load(features_path, allow_pickle=False) as bundle:
        features = bundle['features'].copy()
        ids, starts, ends = (bundle[key] for key in ('internal_id', 'raw_start', 'raw_end'))
        if (any(value.shape != (len(rows),) for value in (ids, starts, ends))
                or ids.dtype.kind not in 'US' or starts.dtype.kind not in 'iu' or ends.dtype.kind not in 'iu'):
            raise ValueError('Window key arrays must be equal-length 1D strings/integers')
        array_keys = list(zip(ids.tolist(), starts.tolist(), ends.tolist()))
    if array_keys != keys:
        raise ValueError('Feature rows and window rows differ or were reordered')
    if (features.shape != (len(rows), config['feature_dimension'])
            or not np.issubdtype(features.dtype, np.floating) or not np.isfinite(features).all()):
        raise ValueError('Invalid feature shape/type/values')
    if fitting:
        actual = {'windows': len(rows), 'sequences': len(sequence_info),
                  'groups': len({r['homology_cluster'] for r in rows})}
        if actual != config['discovery_expectations']:
            raise ValueError(f'Discovery population mismatch: {actual}')
    return features, rows


def assignment_rows(rows: list[dict], labels: np.ndarray, diagnostics: dict | None = None) -> list[dict]:
    result = []
    for i, (row, label) in enumerate(zip(rows, labels)):
        record = {key: row[key] for key in ('internal_id', 'sequence_sha256', 'homology_cluster',
                                           'split', 'raw_start', 'raw_end')}
        record['cluster'] = int(label)
        if diagnostics is not None:
            record.update({key: float(value[i]) for key, value in diagnostics.items() if key != 'cluster'})
        result.append(record)
    return result


def cohort_identity(rows: list[dict]) -> str:
    """Bind ordered window identities and homology groups to a fitted reference."""
    keys = [(r['internal_id'], r['sequence_sha256'], r['homology_cluster'], r['split'],
             int(r['raw_start']), int(r['raw_end'])) for r in rows]
    return hashlib.sha256(json.dumps(keys, separators=(',', ':')).encode('utf-8')).hexdigest()


def checked_model(path: Path, config: dict):
    model = load_model(path, expected_schema=config['feature_schema'])
    if model.config != SpectralConfig(**config['model']):
        raise ValueError('Model parameters do not match the requested configuration')
    return model


def fit_command(features: np.ndarray, rows: list[dict], config: dict, output: Path) -> dict:
    model = fit_spectral(features, SpectralConfig(**config['model']))
    model.diagnostics['training_cohort_sha256'] = cohort_identity(rows)
    save_model(model, output / 'model.npz', config['feature_schema'])
    write_csv(output / 'assignments.csv', assignment_rows(rows, model.labels))
    summary = {'operation': 'spectral_fit', 'windows': len(rows),
               'sequences': len({r['internal_id'] for r in rows}),
               'groups': len({r['homology_cluster'] for r in rows}),
               'cluster_window_counts': np.bincount(model.labels).tolist(),
               'model_diagnostics': model.diagnostics,
               'feature_schema': config['feature_schema'],
               'measured_function_accuracy': None, 'validation_used': False}
    write_json(output / 'summary.json', summary)
    return summary


def predict_command(features: np.ndarray, rows: list[dict], config: dict, model_path: Path,
                    output: Path) -> dict:
    model = checked_model(model_path, config)
    prediction = model.predict(features)
    exported = assignment_rows(rows, prediction['cluster'], prediction)
    for row in exported:
        row['candidate_cluster'] = row.pop('cluster')
    write_csv(output / 'predictions.csv', exported)
    summary = {'operation': 'frozen_knn_label_transfer', 'windows': len(rows),
               'cluster_window_counts': np.bincount(prediction['cluster'], minlength=model.config.n_clusters).tolist(),
               'vote_fraction_median': float(np.median(prediction['vote_fraction'])),
               'vote_fraction_is_calibrated_probability': False,
               'out_of_distribution_rejection_calibrated': False,
               'assignment_is_candidate_only': True,
               'validation_used': False, 'model_sha256': sha256_file(model_path)}
    write_json(output / 'summary.json', summary)
    return summary


def evaluate_command(features: np.ndarray, rows: list[dict], config: dict, model_path: Path,
                     output: Path) -> dict:
    """Measure transfer agreement with a transductive reference, never function accuracy.

    Reference labels include the held-out groups. That is deliberate for a
    reproducibility diagnostic, so its result is not an independent test score.
    Each predictor's PCA, graph and labels use training groups only. Label-ID
    matching uses training rows only, before applying it to held-out predictions.
    """
    reference = checked_model(model_path, config)
    if reference.diagnostics.get('training_cohort_sha256') != cohort_identity(rows):
        raise ValueError('Reference training window identities/groups do not match')
    # A reference from another cohort/order would invalidate the agreement metric.
    projected = (features.astype(np.float64)-reference.pca_mean) @ reference.components.T
    if projected.shape != reference.training_scores.shape or not np.allclose(projected, reference.training_scores, atol=1e-9, rtol=1e-9):
        raise ValueError('Reference model is not fitted to this ordered feature cohort')
    groups = np.array([r['homology_cluster'] for r in rows])
    k = reference.config.n_clusters
    run_rows, predicted_rows, alignment_rows = [], [], []
    settings = config['evaluation']
    for seed in settings['seeds']:
        splitter = GroupShuffleSplit(n_splits=1, test_size=settings['heldout_group_fraction'], random_state=seed)
        train, heldout = next(splitter.split(features, groups=groups))
        train_groups, heldout_groups = set(groups[train]), set(groups[heldout])
        if train_groups & heldout_groups:
            raise ValueError('Homology group leakage')
        predictor = fit_spectral(features[train], reference.config)
        predictor.diagnostics['training_cohort_sha256'] = cohort_identity([rows[i] for i in train])
        contingency = np.zeros((k, k), dtype=int)
        np.add.at(contingency, (predictor.labels, reference.labels[train]), 1)
        source_ids, reference_ids = linear_sum_assignment(-contingency)
        relabel = np.empty(k, dtype=int)
        relabel[source_ids] = reference_ids
        alignment_rows.extend({'seed': seed, 'local_cluster': int(local), 'reference_cluster': int(target),
                               'training_overlap_windows': int(contingency[local, target])}
                              for local, target in zip(source_ids, reference_ids))
        prediction = predictor.predict(features[heldout])
        matched = relabel[prediction['cluster']]
        correct = matched == reference.labels[heldout]
        group_agreement = [float(correct[groups[heldout] == group].mean()) for group in sorted(heldout_groups)]
        run_rows.append({'seed': seed, 'train_windows': len(train), 'heldout_windows': len(heldout),
                         'train_groups': len(train_groups), 'heldout_groups': len(heldout_groups),
                         'group_overlap': 0, 'train_refit_ARI': float(adjusted_rand_score(reference.labels[train], predictor.labels)),
                         'heldout_transfer_ARI': float(adjusted_rand_score(reference.labels[heldout], matched)),
                         'heldout_window_agreement': float(correct.mean()),
                         'heldout_equal_group_agreement': float(np.mean(group_agreement)),
                         'median_vote_fraction': float(np.median(prediction['vote_fraction']))})
        for local_i, original_i in enumerate(heldout):
            row = rows[original_i]
            predicted_rows.append({'seed': seed, 'internal_id': row['internal_id'],
                                   'homology_cluster': row['homology_cluster'],
                                   'raw_start': row['raw_start'], 'raw_end': row['raw_end'],
                                   'reference_cluster': int(reference.labels[original_i]),
                                   'predicted_cluster_matched_to_reference': int(matched[local_i]),
                                   'vote_fraction': float(prediction['vote_fraction'][local_i]),
                                   'nearest_distance': float(prediction['nearest_distance'][local_i])})
        save_model(predictor, output / f'heldout_seed{seed}_model.npz', config['feature_schema'])
    write_csv(output / 'group_holdout_metrics.csv', run_rows)
    write_csv(output / 'group_holdout_predictions.csv', predicted_rows)
    write_csv(output / 'label_alignment.csv', alignment_rows)
    summary = {'operation': 'discovery_group_holdout_transfer_consistency', 'runs': run_rows,
               'reference_includes_heldout_groups': True, 'label_matching_uses': 'training rows only',
               'independent_function_accuracy': False, 'validation_used': False,
               'note': 'Agreement with a full-discovery transductive partition, not independent biological truth.'}
    write_json(output / 'summary.json', summary)
    return summary


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['fit', 'predict', 'evaluate'])
    parser.add_argument('--features', type=Path, required=True)
    parser.add_argument('--windows', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True, help='Hash and feature-schema manifest')
    parser.add_argument('--config', type=Path, default=ROOT / 'configs/m3_spectral_k12.json')
    parser.add_argument('--model', type=Path, help='Required for predict/evaluate')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=3)
    args = parser.parse_args(argv)
    if args.command != 'fit' and args.model is None:
        parser.error('--model is required for predict/evaluate')
    if args.threads < 1:
        parser.error('--threads must be positive')
    config, source = read_json(args.config), read_json(args.source)
    features, rows = load_bundle(args.features, args.windows, source, config, fitting=args.command != 'predict')
    # Formal model runs start at a clean source commit. Ignored output folders are allowed.
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True).strip()
    if dirty:
        raise ValueError('Run requires clean committed source')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    args.output.mkdir(parents=True, exist_ok=False)
    receipt = {'status': 'started', 'operation': args.command, 'commit': commit, 'dirty': False,
               'command': [sys.executable, *sys.argv],
               'utc': datetime.now(timezone.utc).isoformat(),
               'inputs_sha256': {'features': sha256_file(args.features), 'windows': sha256_file(args.windows),
                                 'source': sha256_file(args.source), 'config': sha256_file(args.config)},
               'software': {name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'scikit-learn')},
               'threads': args.threads, 'validation_used': False}
    if args.model:
        receipt['inputs_sha256']['model'] = sha256_file(args.model)
    write_json(args.output / 'run_receipt.json', receipt)
    try:
        with threadpool_limits(limits=args.threads):
            if args.command == 'fit':
                fit_command(features, rows, config, args.output)
            elif args.command == 'predict':
                predict_command(features, rows, config, args.model, args.output)
            else:
                evaluate_command(features, rows, config, args.model, args.output)
    except Exception as exc:
        receipt.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        write_json(args.output / 'run_receipt.json', receipt)
        raise
    receipt.update(status='complete', outputs_sha256={p.name: sha256_file(p) for p in sorted(args.output.iterdir()) if p.name != 'run_receipt.json'})
    write_json(args.output / 'run_receipt.json', receipt)
    print(json.dumps({'status': receipt['status'], 'operation': args.command, 'output': str(args.output)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
