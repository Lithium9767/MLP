#!/usr/bin/env python3
"""One frozen K12 validation: embed, prepare identical windows, then transfer.

Each stage uses a fresh output directory and records its immutable inputs.
Only the embedding stage needs torch/transformers. This entry does not fit PCA,
build a new spectral graph, change K, or use validation to set a threshold.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from bioinformatics.m2_pf00741 import csv_rows
from features.spectral_windows import (load_primary_records, load_hmm_anchors,
                                      accepted_reference_mapping, build_window_features)
from scripts.m3_spectral import (sha256_file, read_json, write_json, write_csv,
                                 load_bundle, checked_model, predict_command, cohort_identity)


def checked_protocol(args):
    """Bind user-authorized validation to specific discovery artifacts first."""
    freeze, config = read_json(args.freeze), read_json(args.config)
    if not freeze.get('authorization', {}).get('quote') or freeze['evaluation']['refit_PCA_or_clustering']:
        raise ValueError('Explicit frozen validation protocol required')
    if sha256_file(args.config) != freeze['config_sha256'] or sha256_file(args.model) != freeze['spectral_model_sha256']:
        raise ValueError('Frozen config or learned model changed')
    if not freeze.get('implementation_sha256'):
        raise ValueError('Frozen implementation identities required')
    for name, digest in freeze['implementation_sha256'].items():
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or sha256_file(path) != digest:
            raise ValueError(f'Frozen implementation changed: {name}')
    for name in ('metadata', 'split_manifest', 'coordinates', 'reference'):
        if sha256_file(getattr(args, name)) != freeze['inputs_sha256'][name]:
            raise ValueError(f'Frozen input mismatch: {name}')
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True).strip()
    if dirty:
        raise ValueError('Validation requires clean committed source')
    return freeze, config


def checked_population(args, freeze):
    validation = load_primary_records(args.metadata, args.split_manifest, 'validation')
    discovery = load_primary_records(args.metadata, args.split_manifest, 'discovery')
    counts = {}
    for name, records in [('validation', validation), ('discovery', discovery)]:
        counts[name] = {'sequences': len(records), 'homology_groups': len({r['homology_cluster'] for r in records})}
        if counts[name] != freeze[f'expected_{name}']:
            raise ValueError(f'Unexpected {name} population')
    for key in ('internal_id', 'sequence_sha256', 'homology_cluster'):
        if {r[key] for r in validation} & {r[key] for r in discovery}:
            raise ValueError(f'Discovery/validation leakage in {key}')
    return validation, counts


def checked_stage(path, freeze, *, stage):
    receipt = read_json(path / 'run_receipt.json')
    if (receipt['status'] != 'complete' or receipt['stage'] != stage
            or receipt['freeze_sha256'] != freeze):
        raise ValueError('Upstream stage was not completed under this freeze')
    for name, digest in receipt['outputs_sha256'].items():
        if sha256_file(path / name) != digest:
            raise ValueError(f'Upstream output changed: {name}')
    return receipt


def embed_stage(args, freeze, records):
    import torch
    from transformers import AutoTokenizer, AutoModel
    settings = freeze['esm']
    versions = {name: importlib.metadata.version(name) for name in settings['software']}
    if versions != settings['software']:
        raise ValueError(f'Embedding runtime differs from discovery: {versions}')
    for name, digest in settings['files'].items():
        if sha256_file(args.esm_dir / name) != digest:
            raise ValueError(f'ESM model/tokenizer changed: {name}')
    torch.set_num_threads(settings['threads'])
    torch.manual_seed(42)
    torch.use_deterministic_algorithms(True)
    tokenizer = AutoTokenizer.from_pretrained(str(args.esm_dir), local_files_only=True)
    model = AutoModel.from_pretrained(str(args.esm_dir), local_files_only=True, use_safetensors=True)
    model.eval()
    model.to('cpu')
    if model.config.hidden_size != settings['dimension'] or model.config.num_hidden_layers != settings['layer']:
        raise ValueError('Unexpected ESM architecture')
    cache = args.output / 'cache'
    cache.mkdir()
    manifest = []
    repeat_difference = None
    for offset in range(0, len(records), settings['batch_size']):
        batch = records[offset:offset+settings['batch_size']]
        encoded = tokenizer([r['sequence'] for r in batch], return_tensors='pt', padding=True, truncation=False)
        for index, row in enumerate(batch):
            length = len(row['sequence'])
            tokens = tokenizer.convert_ids_to_tokens(encoded['input_ids'][index].tolist())
            if ''.join(tokens[1:length+1]) != row['sequence'] or int(encoded['attention_mask'][index].sum()) != length+2:
                raise ValueError('Token/residue mismatch')
        with torch.inference_mode():
            output = model(**encoded).last_hidden_state
            if repeat_difference is None:
                repeated = model(**encoded).last_hidden_state
                repeat_difference = float((output-repeated).abs().max())
                if repeat_difference > 1e-6:
                    raise ValueError('Repeat inference drift')
        for index, row in enumerate(batch):
            matrix = output[index, 1:len(row['sequence'])+1].cpu().numpy().astype(np.float32)
            if matrix.shape != (len(row['sequence']), settings['dimension']) or not np.isfinite(matrix).all():
                raise ValueError('Invalid ESM output')
            name = f"cache/{row['sequence_sha256']}.npz"
            np.savez_compressed(args.output / name, residue_embedding=matrix,
                                full_length_mean=matrix.mean(axis=0),
                                sequence_sha256=np.array(row['sequence_sha256']),
                                model_revision=np.array(settings['revision']))
            manifest.append({'internal_id': row['internal_id'], 'sequence_sha256': row['sequence_sha256'],
                             'split': 'validation', 'analysis_cohort': 'primary',
                             'sequence_length': len(row['sequence']), 'embedding_path': name,
                             'embedding_sha256': sha256_file(args.output / name), 'dimension': settings['dimension']})
        print(f"EMBED {min(offset+len(batch),len(records))}/{len(records)}", flush=True)
    write_csv(args.output / 'embedding_manifest.csv', manifest)
    return {'software': versions, 'sequences': len(records), 'max_repeat_difference': repeat_difference,
            'ESM_weights_updated': False}


def prepare_stage(args, freeze, config, records):
    if np.__version__ != freeze['esm']['software']['numpy'] or importlib.metadata.version('biopython') != freeze['mapping']['biopython_version']:
        raise ValueError('Feature/mapping runtime differs from discovery')
    checked_stage(args.embeddings, sha256_file(args.freeze), stage='embed')
    manifest = csv_rows(args.embeddings / 'embedding_manifest.csv', {'internal_id', 'sequence_sha256', 'embedding_path', 'embedding_sha256', 'split'})
    if len(manifest) != len(records) or {r['internal_id'] for r in manifest} != {r['internal_id'] for r in records}:
        raise ValueError('Embedding manifest does not cover all validation records')
    by_id = {r['internal_id']: r for r in records}
    embeddings = {}
    for row in manifest:
        expected = by_id[row['internal_id']]
        path = (args.embeddings / row['embedding_path']).resolve()
        if not path.is_relative_to(args.embeddings.resolve()) or sha256_file(path) != row['embedding_sha256']:
            raise ValueError('Invalid embedding path/hash')
        if row['split'] != 'validation' or row['sequence_sha256'] != expected['sequence_sha256']:
            raise ValueError('Invalid embedding sequence/split')
        with np.load(path, allow_pickle=False) as archive:
            if str(archive['sequence_sha256']) != expected['sequence_sha256'] or str(archive['model_revision']) != freeze['esm']['revision']:
                raise ValueError('Embedding archive identity mismatch')
            embeddings[row['internal_id']] = archive['residue_embedding'].copy()
    reference_rows = csv_rows(args.reference, {'deposited_position', 'residue', 'hmm_match_state', 'modeled'})
    if [int(r['deposited_position']) for r in reference_rows] != list(range(1, 89)):
        raise ValueError('7R1C reference positions changed')
    if {int(r['hmm_match_state']): int(r['deposited_position']) for r in reference_rows if r['hmm_match_state']} != {h: h+10 for h in range(1, 40)}:
        raise ValueError('7R1C HMM reference anchors changed')
    reference = ''.join(r['residue'] for r in reference_rows)
    anchors = load_hmm_anchors(args.coordinates, records)
    mappings, quality_rows, mapping_rows = {}, [], []
    for row in records:
        sid = row['internal_id']
        mappings[sid], quality = accepted_reference_mapping(row['sequence'], anchors[sid], reference)
        quality_rows.append({'internal_id': sid, 'present_HMM_anchors': len(anchors[sid]),
                             'accepted_reference_residues': len(mappings[sid]),
                             **{f'{name}_{key}': value for name in ('N', 'alpha2', 'C_context')
                                for key, value in quality[name].items()}})
        mapping_rows.extend({'internal_id': sid, 'raw_position': raw, 'reference_position': ref, 'split': 'validation'}
                            for raw, ref in sorted(mappings[sid].items()))
    rows, features, all_rows = build_window_features(records, embeddings, anchors, mappings,
                                                   width=config['window_length'], step=config['window_step'],
                                                   min_accepted=config['minimum_accepted_mapping'])
    write_csv(args.output / 'windows.csv', rows)
    write_csv(args.output / 'all_windows.csv', all_rows)
    write_csv(args.output / 'sequence_mapping_quality.csv', quality_rows)
    write_csv(args.output / 'accepted_reference_mapping.csv', mapping_rows)
    np.savez_compressed(args.output / 'features.npz', features=features,
                        internal_id=np.array([r['internal_id'] for r in rows]),
                        raw_start=np.array([int(r['raw_start']) for r in rows]),
                        raw_end=np.array([int(r['raw_end']) for r in rows]))
    write_json(args.output / 'source.json', {'feature_schema': config['feature_schema'],
               'features_sha256': sha256_file(args.output / 'features.npz'),
               'windows_sha256': sha256_file(args.output / 'windows.csv'),
               'freeze_sha256': sha256_file(args.freeze), 'validation_used': True,
               'purpose': 'Frozen validation transfer; no parameter selection'})
    return {'all_windows': len(all_rows), 'selected_windows': len(rows),
            'sequences_with_selected_windows': len({r['internal_id'] for r in rows}), 'fit_performed': False}


def transfer_stage(args, freeze, config, records):
    from models.spectral_validation import summarize_validation
    from threadpoolctl import threadpool_limits
    versions = {name: importlib.metadata.version(name) for name in freeze['transfer_software']}
    if versions != freeze['transfer_software']:
        raise ValueError('Transfer runtime differs from frozen discovery model')
    checked_stage(args.prepared, sha256_file(args.freeze), stage='prepare')
    for path, key in [(args.discovery_features, 'discovery_features'), (args.discovery_windows, 'discovery_windows')]:
        if sha256_file(path) != freeze['inputs_sha256'][key]:
            raise ValueError(f'Frozen discovery source changed: {key}')
    source = read_json(args.prepared / 'source.json')
    features, rows = load_bundle(args.prepared / 'features.npz', args.prepared / 'windows.csv', source, config,
                                fitting=False, prediction_splits=frozenset({'validation'}))
    if not {r['internal_id'] for r in rows} <= {r['internal_id'] for r in records}:
        raise ValueError('Prepared windows include unknown validation sequences')
    discovery_source = {'feature_schema': config['feature_schema'],
                        'features_sha256': freeze['inputs_sha256']['discovery_features'],
                        'windows_sha256': freeze['inputs_sha256']['discovery_windows']}
    train_features, train_rows = load_bundle(args.discovery_features, args.discovery_windows, discovery_source, config, fitting=True)
    model = checked_model(args.model, config)
    if cohort_identity(train_rows) != model.diagnostics['training_cohort_sha256']:
        raise ValueError('Discovery training identities changed')
    with threadpool_limits(limits=3):
        predict_command(features, rows, config, args.model, args.output)
        prediction = model.predict(features)
        summary, clusters, diagnostics = summarize_validation(train_rows, train_features, model.labels,
              rows, features, prediction, model, seed=freeze['evaluation']['permutation_seed'],
              permutations=freeze['evaluation']['within_sequence_label_permutations'])
    write_json(args.output / 'validation_summary.json', summary)
    for row in clusters:
        for key in ('discovery_reference_profile_1_66_and_unmapped', 'validation_reference_profile_1_66_and_unmapped'):
            row[key] = json.dumps(row[key], allow_nan=False)
    write_csv(args.output / 'cluster_transfer_summary.csv', clusters)
    write_csv(args.output / 'window_diagnostics.csv', diagnostics)
    if sha256_file(args.model) != freeze['spectral_model_sha256']:
        raise ValueError('Frozen model changed during inference')
    return {'selected_windows': len(rows), 'selected_sequences': len({r['internal_id'] for r in rows}),
            'fit_performed': False, 'frozen_model_unchanged': True, 'validation_used': True,
            'software': versions}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('embed', 'prepare', 'transfer'))
    for name in ('metadata', 'split-manifest', 'coordinates', 'reference', 'model', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--freeze', type=Path, default=ROOT / 'configs/m3_spectral_validation_freeze.json')
    parser.add_argument('--config', type=Path, default=ROOT / 'configs/m3_spectral_k12.json')
    for name in ('esm-dir', 'embeddings', 'prepared', 'discovery-features', 'discovery-windows'):
        parser.add_argument('--'+name, type=Path)
    args = parser.parse_args(argv)
    requirements = {'embed': ('esm_dir',), 'prepare': ('embeddings',),
                    'transfer': ('prepared', 'discovery_features', 'discovery_windows')}
    if any(getattr(args, name) is None for name in requirements[args.stage]):
        parser.error(f'{args.stage} requires {requirements[args.stage]}')
    freeze, config = checked_protocol(args)
    records, population = checked_population(args, freeze)
    args.output.mkdir(parents=True, exist_ok=False)
    receipt = {'status': 'started', 'stage': args.stage, 'validation_used': True,
               'freeze_sha256': sha256_file(args.freeze), 'protocol': freeze['protocol'],
               'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
               'dirty': False, 'command': [sys.executable, *sys.argv],
               'utc': datetime.now(timezone.utc).isoformat(), 'population': population,
               'software': {'python': sys.version, 'numpy': np.__version__},
               'frozen_inputs_sha256': freeze['inputs_sha256'], 'fit_performed': False,
               'spectral_model_sha256': freeze['spectral_model_sha256']}
    write_json(args.output / 'run_receipt.json', receipt)
    try:
        if args.stage == 'embed':
            result = embed_stage(args, freeze, records)
        elif args.stage == 'prepare':
            result = prepare_stage(args, freeze, config, records)
        else:
            result = transfer_stage(args, freeze, config, records)
    except Exception as exc:
        receipt.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        write_json(args.output / 'run_receipt.json', receipt)
        raise
    receipt.update(status='complete', result=result,
                   outputs_sha256={p.relative_to(args.output).as_posix(): sha256_file(p)
                       for p in sorted(args.output.rglob('*')) if p.is_file() and p.name != 'run_receipt.json'})
    write_json(args.output / 'run_receipt.json', receipt)
    print({'status': 'complete', 'stage': args.stage, **result}, flush=True)


if __name__ == '__main__':
    main()
