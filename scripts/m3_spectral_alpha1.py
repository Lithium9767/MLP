"""Audit the nominated alpha1/N-terminal branch with group omission and window shifts."""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path, PureWindowsPath

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from models.spectral_validation import _balanced_weights
from scripts import m3_spectral as spectral


def checked_csv(path, expected):
    if spectral.sha256_file(path) != expected:
        raise ValueError(f'Input digest mismatch: {path.name}')
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def fit_path(features, rows, levels):
    """Bind branch identity to its discovery members; cluster numbers can permute."""
    selected = np.arange(len(rows))
    path = []
    for level in levels:
        model = spectral.fit_spectral(features[selected], spectral.SpectralConfig(n_clusters=level['k']))
        for label in range(level['k']):
            members = selected[model.labels == label]
            if len(members) == level['windows']:
                if spectral.cohort_identity([rows[i] for i in members]) == level['cohort_sha256']:
                    break
        else:
            raise ValueError(f"Cannot reproduce the archived branch at K={level['k']}")
        path.append((model, label, selected, members))
        selected = members
    return path


def route(path, features):
    active = np.arange(len(features))
    for model, label, _, _ in path:
        previous = active
        active = previous[model.predict(features[previous])['cluster'] == label]
    return previous, active


def matched_rates(rows, mask, classes, strata, indices):
    """Weight all parent rows first, then standardize the two arms to common candidate strata."""
    weights = _balanced_weights(rows, indices)
    # Each cell accumulates total weight and S/T weight per arm.
    totals = defaultdict(lambda: np.zeros((2, 2)))
    groups = defaultdict(lambda: [set(), set()])
    for i, weight in zip(indices, weights):
        if classes[i] is None:
            continue
        arm = 0 if mask[i] else 1
        cell = strata[i]
        totals[cell][arm] += (weight, weight * classes[i])
        groups[cell][arm].add(rows[i]['homology_cluster'])
    shared = [cell for cell, value in totals.items() if np.all(value[:, 0] > 0)]
    support = [set().union(*(groups[cell][arm] for cell in shared)) for arm in (0, 1)]
    summary = {'groups_candidate': len(support[0]), 'groups_siblings': len(support[1]),
               'candidate_rate': None, 'sibling_rate': None, 'difference': None, 'comparable_mass': 0.0}
    if shared:
        common_mass = sum(totals[cell][0, 0] for cell in shared)
        rates = sum(totals[cell][0, 0] * totals[cell][:, 1] / totals[cell][:, 0]
                    for cell in shared) / common_mass
        summary.update(candidate_rate=float(rates[0]), sibling_rate=float(rates[1]),
                       difference=float(rates[0] - rates[1]),
                       comparable_mass=float(common_mass / sum(value[0, 0] for value in totals.values())))
    return summary, support


def chemistry(rows, mask, mappings, site, residues):
    classes, footprints = [], []
    for row in rows:
        start, end = int(row['raw_start']), int(row['raw_end'])
        raw = next((pos for pos in range(start, end + 1)
                    if mappings[row['internal_id']].get(pos) == site), None)
        classes.append(None if raw is None else row['fragment'][raw - start] in residues)
        length = int(row['sequence_length'])
        footprints.append((row['reference_positions'], int(length > 100) + int(length > 130)))
    all_rows = np.arange(len(rows))
    output = []
    for within in (False, True):
        strata = [(fp, row['homology_cluster']) for fp, row in zip(footprints, rows)] if within else footprints
        full, support = matched_rates(rows, mask, classes, strata, all_rows)
        differences = []
        for group in sorted(support[0] | support[1]):
            kept = np.array([i for i, r in enumerate(rows) if r['homology_cluster'] != group])
            reduced, _ = matched_rates(rows, mask, classes, strata, kept)
            if reduced['difference'] is not None:
                differences.append(reduced['difference'])
        output.append(dict(within_homology_group=within, **full, omissions_with_estimate=len(differences),
                           omission_min=min(differences, default=None), omission_max=max(differences, default=None)))
    return output


def shifted_features(rows, original, mappings, manifest, embeddings, offsets):
    """Only accepted reference positions are available to the shift diagnostic."""
    data = {offset: [] for offset in offsets}
    kept = {offset: [] for offset in offsets}
    by_sequence = defaultdict(list)
    for i, row in enumerate(rows):
        by_sequence[row['internal_id']].append(i)
    for sid, indices in by_sequence.items():
        record = manifest[sid]
        # Resolve the archived cache suffix against the supplied directory,
        # including when a Windows producer's manifest is opened on Linux.
        parts = PureWindowsPath(record['embedding_path']).parts
        if 'cache' not in parts or '..' in parts:
            raise ValueError('Expected a path within the supplied embedding cache')
        location = embeddings.joinpath(*parts[parts.index('cache'):])
        if spectral.sha256_file(location) != record['embedding_sha256']:
            raise ValueError('Embedding hash mismatch')
        with np.load(location, allow_pickle=False) as archive:
            embedding = archive['residue_embedding']
            if str(archive['sequence_sha256']) != rows[indices[0]]['sequence_sha256']:
                raise ValueError('Embedding sequence mismatch')
        mapping = mappings[sid]
        anchors = sorted(raw for raw, ref in mapping.items() if 11 <= ref <= 49)
        background = embedding[np.array(anchors) - 1].mean(axis=0)
        for i in indices:
            start = int(rows[i]['raw_start'])
            for offset in offsets:
                begin = start + offset
                if begin < 1 or begin + 9 > len(embedding):
                    continue
                if sum(pos in mapping for pos in range(begin, begin + 10)) < 8:
                    continue
                kept[offset].append(i)
                data[offset].append(embedding[begin - 1:begin + 9].mean(axis=0) - background)
    if sorted(kept[0]) != list(range(len(rows))) or not np.array_equal(original[kept[0]], data[0]):
        raise ValueError('Zero-shift features do not reproduce the original candidate windows')
    return kept, data


def main(args):
    lock = spectral.read_json(ROOT / 'configs/m3_spectral_alpha1.json')
    config = spectral.read_json(ROOT / 'configs/m3_spectral_k12.json')
    args.output.mkdir(parents=True, exist_ok=False)
    stats, shifts, populations = [], [], {}
    for split in ('discovery', 'validation'):
        spec = lock[split]
        bundle = getattr(args, split)
        x, rows = spectral.load_bundle(bundle / 'features.npz', bundle / spec['window_table'],
            dict(spec, feature_schema=config['feature_schema']), config, fitting=split == 'discovery',
            prediction_splits=frozenset({'validation'}))
        populations[split] = [{r[key] for r in rows} for key in ('internal_id', 'sequence_sha256', 'homology_cluster')]
        if split == 'discovery':
            path = fit_path(x, rows, lock['path'])
            parent, candidate = path[-1][2:]
        else:
            if any(a & b for a, b in zip(populations['discovery'], populations['validation'])):
                raise ValueError('Cohort overlap')
            parent, candidate = route(path, x)
        mapping_rows = checked_csv(getattr(args, split + '_mapping'), spec['mapping_sha256'])
        mappings = defaultdict(dict)
        for record in mapping_rows:
            if record['reference_position']:
                mappings[record['internal_id']][int(record['raw_position'])] = int(record['reference_position'])
        parent_rows = [rows[i] for i in parent]
        stats.extend(dict(split=split, **item) for item in chemistry(parent_rows, np.isin(parent, candidate),
                     mappings, lock['site'], lock['residues']))
        embeddings = getattr(args, split + '_embeddings')
        manifest = checked_csv(embeddings / 'embedding_manifest.csv', spec['embedding_manifest_sha256'])
        candidate_rows = [rows[i] for i in candidate]
        offsets = (0, -1, 1, -2, 2)
        kept, data = shifted_features(candidate_rows, x[candidate], mappings,
                         {r['internal_id']: r for r in manifest}, embeddings, offsets)
        for offset in offsets:
            indices = kept[offset]
            retained = []
            if indices:
                _, retained = route(path, np.asarray(data[offset]))
            rate = float(_balanced_weights(candidate_rows, indices)[retained].sum()) if indices else None
            shifts.append(dict(split=split, offset=offset, original=len(candidate), eligible=len(indices),
                               retained=len(retained), group_balanced_retention=rate))
    spectral.write_csv(args.output / 'chemistry.csv', stats)
    spectral.write_csv(args.output / 'window_shift.csv', shifts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for split in ('discovery', 'validation'):
        for suffix in ('', '-mapping', '-embeddings'):
            parser.add_argument('--' + split + suffix, type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    with threadpool_limits(limits=3):
        main(parser.parse_args())
