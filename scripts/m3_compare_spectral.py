"""Paired comparison of flat and hierarchical spectral partitions on the same windows."""
import argparse
import hashlib
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from models.spectral_regions import SpectralConfig, fit_spectral
from models.spectral_validation import _balanced_weights
from scripts import m3_spectral as io


def partitions(train, evaluation):
    results = {}
    for k in (6, 12):
        model = fit_spectral(train, SpectralConfig(n_clusters=k))
        results[f'flat{k}'] = model.labels, model.predict(evaluation)['cluster']
    root = fit_spectral(train, SpectralConfig(n_clusters=3))
    branches = root.predict(evaluation)['cluster']
    train_labels = np.empty(len(train), dtype=int)
    eval_labels = np.empty(len(evaluation), dtype=int)
    for label in range(3):
        group = root.labels == label
        target = branches == label
        child = fit_spectral(train[group], SpectralConfig(n_clusters=2))
        train_labels[group] = 2 * label + child.labels
        eval_labels[target] = 2 * label + child.predict(evaluation[target])['cluster']
    results['hierarchical6'] = train_labels, eval_labels
    return results


def archive_numbering(pair, reference):
    """Recover presentation IDs by discovery size, then verify both full label vectors."""
    train, _ = pair
    sizes = reference['label_counts']
    if len(set(sizes)) != len(sizes) or sorted(np.bincount(train)) != sorted(sizes):
        raise ValueError('Discovery partition differs from the archived run')
    mapping = np.array([sizes.index(n) for n in np.bincount(train)])
    result = tuple(mapping[v] for v in pair)
    for split, labels in zip(('discovery', 'validation'), result):
        digest = hashlib.sha256(labels.astype('uint8').tobytes()).hexdigest()
        if digest != reference[split + '_sha256']:
            raise ValueError(f'{split} assignments do not reproduce the archived run')
    return result


def midpoints(rows):
    return np.array([np.median([int(n) for n in r['reference_positions'].split(';')]) for r in rows])


def position_error(train_rows, eval_rows, trained, predicted):
    """Use only the cluster ID to predict a mapped midpoint. Fit medians on discovery."""
    train_position = midpoints(train_rows)
    eval_position = midpoints(eval_rows)
    estimate = np.empty(len(eval_rows))
    for cluster in np.unique(trained):
        idx = np.flatnonzero(trained == cluster)
        idx = idx[np.argsort(train_position[idx], kind='stable')]
        cumulative = _balanced_weights(train_rows, idx).cumsum()
        median = train_position[idx[np.searchsorted(cumulative, .5)]]
        estimate[predicted == cluster] = median
    weights = _balanced_weights(eval_rows, np.arange(len(eval_rows)))
    return float(weights @ np.abs(estimate - eval_position))


def main(args):
    settings = io.read_json(ROOT / 'configs/m3_spectral_comparison.json')
    source = io.read_json(ROOT / 'configs/m3_spectral_alpha1.json')
    config = io.read_json(ROOT / 'configs/m3_spectral_k12.json')
    inputs = []
    for split in ('discovery', 'validation'):
        spec = source[split]
        folder = getattr(args, split)
        inputs.append(io.load_bundle(folder / 'features.npz', folder / spec['window_table'],
            dict(spec, feature_schema=config['feature_schema']), config, fitting=split == 'discovery',
            prediction_splits=frozenset({'validation'})))
    (x, rows), (v, vr) = inputs
    for field in ('sequence_sha256', 'internal_id', 'homology_cluster'):
        if {r[field] for r in rows} & {r[field] for r in vr}:
            raise ValueError(f'Cohorts overlap on {field}')
    args.output.mkdir(parents=True, exist_ok=False)
    full = partitions(x, v)
    for name, reference in settings['archived'].items():
        full[name] = archive_numbering(full[name], reference)

    crosswalk = []
    for split_index, split in enumerate(('discovery', 'validation')):
        labels_six = full['hierarchical6'][split_index]
        labels_twelve = full['flat12'][split_index]
        for k in range(12):
            counts = np.bincount(labels_six[labels_twelve == k], minlength=6)
            crosswalk.append(dict(split=split, cluster12=k, windows=int(counts.sum()),
                                  **{f'region_{j//2}_{j%2}': int(n) for j, n in enumerate(counts)}))
    io.write_csv(args.output / 'crosswalk.csv', crosswalk)
    groups = np.array([r['homology_cluster'] for r in rows])
    unique_groups = np.unique(groups)
    predictions = {}
    for seed in settings['seeds']:
        keep = np.random.default_rng(seed).choice(unique_groups, settings['training_groups'], replace=False)
        fitted = partitions(x[np.isin(groups, keep)], v)
        predictions[seed] = {name: pair[1] for name, pair in fitted.items()}
        print(f'Completed paired refit {seed}', flush=True)
    scores = []
    for first, second in combinations(settings['seeds'], 2):
        scores.append(dict(seed_a=first, seed_b=second, **{name: adjusted_rand_score(
            predictions[first][name], predictions[second][name]) for name in full}))
    io.write_csv(args.output / 'paired_refits.csv', scores)
    summary = {}
    for name, (train, predicted) in full.items():
        values = [r[name] for r in scores]
        summary[name] = dict(median_ari=float(np.median(values)), min_ari=float(min(values)),
                             max_ari=float(max(values)),
                             midpoint_mae=position_error(rows, vr, train, predicted))
    io.write_json(args.output / 'summary.json', dict(discovery_windows=len(rows),
        validation_windows=len(vr), paired_runs=len(predictions), pair_comparisons=len(scores),
        methods=summary, archived_labels_verified=list(settings['archived'])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('discovery', 'validation', 'output'):
        parser.add_argument('--' + option, type=Path, required=True)
    with threadpool_limits(limits=3):
        main(parser.parse_args())
