#!/usr/bin/env python3
"""Discovery-only ESM-2 extraction and preregistered region analysis.

Run code from a clean commit. Large outputs stay in ignored data/processed/.
No validation inference, parameter search or attribution is performed here.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from features.esm2_embed import read_embedding_records, iter_windows, sha256_file


def read_csv(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text('', encoding='utf-8'); return
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()


def clean_provenance(config):
    if git('status', '--porcelain'):
        raise ValueError('Formal run requires a clean committed working tree')
    return {'commit': git('rev-parse', 'HEAD'), 'dirty': False,
            'config_sha256': sha256_file(config), 'command': sys.argv,
            'utc': datetime.now(timezone.utc).isoformat(), 'python': sys.version,
            'platform': platform.platform()}


def checked_inputs(args, config):
    for path, key in [(args.metadata, 'metadata_sha256'), (args.split_manifest, 'split_manifest_sha256'),
                      (args.coordinates, 'coordinates_sha256')]:
        if sha256_file(path) != config[key]:
            raise ValueError(f'Frozen input hash mismatch: {path}')
    records = read_embedding_records(args.metadata, args.split_manifest, primary_only=not args.include_sensitivity)
    expected = 1453 if args.include_sensitivity else config['expected_primary_discovery']
    if len(records) != expected:
        raise ValueError('Unexpected discovery primary count')
    if any(r.split != 'discovery' or (not args.include_sensitivity and r.analysis_cohort != 'primary') for r in records):
        raise ValueError('Discovery primary only')
    return records


def checked_coordinates(path, records):
    by_id = {r.internal_id: r for r in records}; mappings = defaultdict(dict)
    for row in read_csv(path):
        key = row['internal_id']
        if key not in by_id: continue
        r = by_id[key]; position = int(row['raw_position'])
        if row['split'] != 'discovery' or row['sequence_sha256'] != r.sequence_sha256:
            raise ValueError('Coordinate identity/split mismatch')
        if position < 1 or position > len(r.sequence) or r.sequence[position-1] != row['residue']:
            raise ValueError('Coordinate residue mismatch')
        if position in mappings[key]: raise ValueError('Ambiguous multi-domain coordinate')
        state = int(row['hmm_match_state']) if row['hmm_match_state'] else None
        if state is not None and not 1 <= state <= 39: raise ValueError('HMM state out of range')
        mappings[key][position] = state
    if set(mappings) != set(by_id): raise ValueError('Missing discovery coordinates')
    return mappings


def window_rows(records, embeddings, mappings, width, step):
    import numpy as np
    vectors, rows, composition = [], [], []
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    for record in records:
        matrix = embeddings[record.internal_id]
        if matrix.shape != (len(record.sequence), 480) or not np.isfinite(matrix).all():
            raise ValueError('Embedding shape/finite check failed')
        for start, end in iter_windows(len(record.sequence), width, step):
            states = sorted({mappings[record.internal_id].get(p) for p in range(start, end+1)} - {None})
            mapped = sum(mappings[record.internal_id].get(p) is not None for p in range(start, end+1))
            rows.append({'internal_id': record.internal_id, 'sequence_sha256': record.sequence_sha256,
                         'analysis_cohort': record.analysis_cohort,
                         'split': 'discovery', 'raw_start': start, 'raw_end': end, 'window_length': width,
                         'sequence_length': len(record.sequence), 'normalized_start': (start-1)/len(record.sequence),
                         'hmm_states': ';'.join(map(str, states)), 'mapped_fraction': mapped/width})
            vectors.append(matrix[start-1:end].mean(axis=0))
            fragment = record.sequence[start-1:end]
            composition.append([fragment.count(aa)/width for aa in alphabet])
    if not vectors: raise ValueError('No windows')
    return rows, np.stack(vectors), np.asarray(composition, dtype=np.float32)


def embed(args, config, records):
    import numpy as np
    import torch
    from transformers import AutoTokenizer, AutoModel
    torch.set_num_threads(args.threads); torch.manual_seed(config['seed'])
    torch.use_deterministic_algorithms(True)
    source = args.model_dir
    receipt = json.loads((source/'download_receipt.json').read_text())
    if receipt['revision'] != config['model_revision'] or receipt['model'] != config['model_id']:
        raise ValueError('Model revision receipt mismatch')
    for f in receipt['files']:
        if sha256_file(source/f['name']) != f['sha256']: raise ValueError('Model file hash mismatch')
    tokenizer = AutoTokenizer.from_pretrained(str(source), local_files_only=True)
    model = AutoModel.from_pretrained(str(source), local_files_only=True, use_safetensors=True)
    model.eval(); model.to('cpu')
    if model.config.hidden_size != 480: raise ValueError('Unexpected model dimension')
    versions = {n: importlib.metadata.version(n) for n in ['numpy','torch','transformers','safetensors']}
    cache_fingerprint = hashlib.sha256(json.dumps(versions, sort_keys=True).encode()).hexdigest()[:12]
    cache = args.output_dir/'cache'/config['model_revision']/cache_fingerprint
    cache.mkdir(parents=True, exist_ok=True)
    selected = records[:args.limit] if args.limit else records
    manifests=[]; repeated=False; max_repeat_difference=None
    for offset in range(0, len(selected), args.batch_size):
        batch = selected[offset:offset+args.batch_size]; pending=[]
        for r in batch:
            p=cache/(r.sequence_sha256+'.npz')
            if p.exists():
                with np.load(p, allow_pickle=False) as z:
                    a=z['residue_embedding']
                    if str(z['sequence_sha256'])!=r.sequence_sha256 or str(z['model_revision'])!=config['model_revision']:
                        raise ValueError('Cached provenance mismatch')
                    if a.shape!=(len(r.sequence),480) or not np.isfinite(a).all():raise ValueError('Invalid cache')
            else: pending.append(r)
        if pending:
            encoded=tokenizer([r.sequence for r in pending],return_tensors='pt',padding=True,truncation=False)
            for i,r in enumerate(pending):
                tokens=tokenizer.convert_ids_to_tokens(encoded['input_ids'][i].tolist())
                if ''.join(tokens[1:1+len(r.sequence)])!=r.sequence:raise ValueError('Token/residue alignment mismatch')
                if int(encoded['attention_mask'][i].sum())!=len(r.sequence)+2:raise ValueError('Unexpected special tokens')
            with torch.inference_mode(): output=model(**encoded).last_hidden_state
            if not repeated:
                with torch.inference_mode(): again=model(**encoded).last_hidden_state
                max_repeat_difference=float((output-again).abs().max())
                if max_repeat_difference>1e-6:raise ValueError('Repeat inference drift')
                repeated=True
            for i,r in enumerate(pending):
                a=output[i,1:1+len(r.sequence)].cpu().numpy().astype(np.float32)
                if a.shape!=(len(r.sequence),480) or not np.isfinite(a).all():raise ValueError('Invalid model output')
                np.savez_compressed(cache/(r.sequence_sha256+'.npz'),residue_embedding=a,
                                    full_length_mean=a.mean(axis=0),sequence_sha256=r.sequence_sha256,
                                    model_revision=config['model_revision'])
        for r in batch:
            p=cache/(r.sequence_sha256+'.npz')
            manifests.append({'internal_id':r.internal_id,'sequence_sha256':r.sequence_sha256,
                              'split':'discovery','analysis_cohort':r.analysis_cohort,'sequence_length':len(r.sequence),
                              'embedding_path':str(p.resolve()),'embedding_sha256':sha256_file(p),'dimension':480})
        print(f'EMBED {min(offset+len(batch),len(selected))}/{len(selected)}',flush=True)
    write_csv(args.output_dir/'embedding_manifest.csv',manifests)
    return {'n_sequences':len(selected),'eligible_parent_count':len(records),'software':versions,
            'device':'cpu','threads':args.threads,'layer':model.config.num_hidden_layers,
            'model_id':config['model_id'],'revision':config['model_revision'],'model_files':receipt['files'],
            'max_repeat_difference':max_repeat_difference,'fresh_repeat_test':repeated,
            'cache_entries_validated':len(manifests),'pooling':'mean of residues; special tokens excluded',
            'embedding_manifest_sha256':sha256_file(args.output_dir/'embedding_manifest.csv')}


def cluster_fit(vectors, dim, min_size, config):
    import numpy as np
    from sklearn.decomposition import PCA
    import hdbscan
    pca=PCA(n_components=min(dim,vectors.shape[0]-1,vectors.shape[1]),svd_solver='randomized',random_state=config['seed'])
    transformed=pca.fit_transform(vectors)
    model=hdbscan.HDBSCAN(min_cluster_size=min_size,min_samples=config['min_samples'],
                         core_dist_n_jobs=1,prediction_data=True)
    labels=model.fit_predict(transformed)
    return pca,model,transformed,labels


def candidate_table(rows, labels, config):
    import numpy as np
    rng=np.random.default_rng(config['seed']);incidence=np.zeros((len(rows),39),dtype=np.int8)
    bins=defaultdict(list)
    for i,r in enumerate(rows):
        for state in r['hmm_states'].split(';'):
            if state:incidence[i,int(state)-1]=1
        bins[(r['sequence_length']//20,int(r['normalized_start']*5))].append(i)
    output=[]
    for label in sorted(set(labels)-{-1}):
        ids=np.flatnonzero(labels==label); values=incidence[ids].mean(axis=0)
        best=int(values.argmax())+1;groups=Counter((rows[i]['sequence_length']//20,int(rows[i]['normalized_start']*5)) for i in ids)
        null=[]
        for _ in range(config['random_baseline_repeats']):
            sampled=[]
            for key,count in groups.items():
                pool=[i for i in bins[key] if labels[i]!=label]
                if not pool:pool=bins[key]
                sampled.extend(rng.choice(pool,size=count,replace=True).tolist())
            null.append(float(incidence[sampled].mean(axis=0).max()))
        support=Counter()
        for i in ids:
            for state in rows[i]['hmm_states'].split(';'):
                if state:support[int(state)]+=1
        active=[state for state,count in support.items() if count/len(ids)>=0.5]
        output.append({'cluster':int(label),'n_windows':len(ids),'n_sequences':len({rows[i]['internal_id'] for i in ids}),
                       'dominant_hmm_state':best,'dominant_state_window_fraction':float(values.max()),
                       'majority_hmm_start':min(active) if active else '', 'majority_hmm_end':max(active) if active else '',
                       'mean_mapped_fraction':float(np.mean([rows[i]['mapped_fraction'] for i in ids])),
                       'matched_random_mean':float(np.mean(null)),
                       'matched_random_p95':float(np.quantile(null,.95)),
                       'random_comparison':'exploratory; overlapping windows are not independent; no functional p-value'})
    return output


def analyze(args,config,records):
    import numpy as np
    import joblib
    from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
    from sklearn.decomposition import PCA
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    manifests=read_csv(args.embeddings/'embedding_manifest.csv'); selected={m['internal_id'] for m in manifests}
    records=[r for r in records if r.internal_id in selected]
    if len(records)!=len(manifests):raise ValueError('Embedding manifest identity mismatch')
    by_id={r.internal_id:r for r in records};embeddings={}
    for m in manifests:
        r=by_id[m['internal_id']]
        if m['split']!='discovery' or m['sequence_sha256']!=r.sequence_sha256:raise ValueError('Embedding split/hash mismatch')
        p=Path(m['embedding_path'])
        if sha256_file(p)!=m['embedding_sha256']:raise ValueError('Embedding file mismatch')
        with np.load(p,allow_pickle=False) as z: embeddings[r.internal_id]=z['residue_embedding']
    mappings=checked_coordinates(args.coordinates,records)
    full=np.stack([embeddings[r.internal_id].mean(axis=0) for r in records]);np.save(args.output_dir/'full_length_mean.npy',full)
    meta={r['internal_id']:r for r in read_csv(args.metadata)}
    full_pca=PCA(n_components=min(20,len(records)-1),random_state=config['seed']).fit_transform(full)
    write_csv(args.output_dir/'full_length_pca.csv',[{'internal_id':r.internal_id,'pc1':float(full_pca[i,0]),'pc2':float(full_pca[i,1])} for i,r in enumerate(records)])
    main_labels=None;main_rows=None;main_transformed=None;sensitivity=[];explained=[];confounds=[];runs={}
    for width in config['window_lengths']:
        rows,vectors,composition=window_rows(records,embeddings,mappings,width,config['window_step'])
        for dim in config['pca_dimensions']:
            for min_size in config['min_cluster_sizes']:
                pca,model,x,labels=cluster_fit(vectors,dim,min_size,config)
                key=f'w{width}-p{dim}-m{min_size}'
                info={'run':key,'window_length':width,'pca_requested':dim,'pca_actual':int(x.shape[1]),
                      'min_cluster_size':min_size,'n_windows':len(rows),'n_clusters':len(set(labels)-{-1}),
                      'noise_fraction':float(np.mean(labels==-1)), 'explained_variance':float(pca.explained_variance_ratio_.sum())}
                sensitivity.append(info);runs[key]={(r['internal_id'],r['raw_start']):int(label) for r,label in zip(rows,labels)}
                target=config['main_parameters']
                if (width,dim,min_size)==(target['window_length'],target['pca_dimension'],target['min_cluster_size']):
                    main_labels=labels;main_rows=rows;main_transformed=x
                    joblib.dump({'pca':pca,'clusterer':model},args.output_dir/'frozen_transform_proposal.joblib')
                    explained=[{'component':i+1,'explained_variance_ratio':float(v)} for i,v in enumerate(pca.explained_variance_ratio_)]
                    baseline=cluster_fit(composition,min(20,composition.shape[1]),min_size,config)[3]
                    write_csv(args.output_dir/'simple_baseline.csv',[{'baseline':'amino_acid_composition','n_clusters':len(set(baseline)-{-1}),
                               'noise_fraction':float(np.mean(baseline==-1)),'ARI_vs_ESM':float(adjusted_rand_score(labels,baseline))}])
                    fields={'length_bin':[r['sequence_length']//20 for r in rows],
                            'position_bin':[int(r['normalized_start']*5) for r in rows],
                            'organism':[meta[r['internal_id']].get('organism','missing') for r in rows],
                            'C_terminal_extension':[int(not r['hmm_states'] and r['normalized_start']>.5) for r in rows]}
                    confounds=[{'factor':factor,'adjusted_mutual_information':float(adjusted_mutual_info_score(labels,values)),
                                'interpretation':'association diagnostic, not causal adjustment; primary contains no partial/type-conflict'} for factor,values in fields.items()]
                    write_csv(args.output_dir/'candidate_windows.csv',[dict(r,cluster=int(label),probability=float(model.probabilities_[i])) for i,(r,label) in enumerate(zip(rows,labels))])
                print('CLUSTER',info,flush=True)
    if main_labels is None:raise ValueError('Main configuration missing')
    main_key='w30-p20-m50';reference=runs[main_key]
    for info in sensitivity:
        current=runs[info['run']];common=sorted(reference.keys()&current.keys())
        info['common_window_starts']=len(common)
        info['ARI_common_starts']=float(adjusted_rand_score([reference[k] for k in common],[current[k] for k in common]))
        info['ARI_limitation']='different widths compare anchored starts, not identical fragments; noise included'
    candidates=candidate_table(main_rows,main_labels,config)
    for name,rows in [('parameter_sensitivity.csv',sensitivity),('pca_explained_variance.csv',explained),
                      ('confound_checks.csv',confounds),('candidate_regions.csv',candidates)]:write_csv(args.output_dir/name,rows)
    fig,ax=plt.subplots();ax.scatter(main_transformed[:,0],main_transformed[:,1],c=main_labels,s=2,cmap='tab20');ax.set(xlabel='PC1',ylabel='PC2',title='Discovery window PCA; exploratory clusters');fig.savefig(args.output_dir/'window_pca.png',dpi=160);plt.close(fig)
    # UMAP display is a deterministic subsample; never used in PCA/HDBSCAN fitting.
    import umap
    rng=np.random.default_rng(config['seed']);index=np.sort(rng.choice(len(main_rows),min(3000,len(main_rows)),replace=False))
    u=umap.UMAP(random_state=config['seed'],n_jobs=1).fit_transform(main_transformed[index])
    fig,ax=plt.subplots();ax.scatter(u[:,0],u[:,1],c=main_labels[index],s=3,cmap='tab20');ax.set(title='UMAP display only; fixed subsample');fig.savefig(args.output_dir/'window_umap.png',dpi=160);plt.close(fig)
    write_json(args.output_dir/'candidate_freeze_proposal.json',{'status':'proposed_not_A_approved','fit_split':'discovery','parameters':config['main_parameters'],
               'candidate_table_sha256':sha256_file(args.output_dir/'candidate_regions.csv'),
               'transform_sha256':sha256_file(args.output_dir/'frozen_transform_proposal.joblib'),
               'validation_enabled':False,'limits':['Candidate clusters are not functional proof','No validation inference before A approves this concrete freeze proposal','Random diagnostics are exploratory with overlapping windows']})
    return {'n_sequences':len(records),'n_windows':len(main_rows),'n_candidate_clusters':len(candidates),
            'n_parameter_runs':len(sensitivity),'validation_used':False,'primary_only':True,
            'sensitivity_cohorts':'partial/type-conflict not fitted; require separate follow-up cohort run',
            'software':{n:importlib.metadata.version(n) for n in ['numpy','scikit-learn','hdbscan','umap-learn','matplotlib']}}


def cohort(args,config,records):
    import numpy as np
    import joblib,hdbscan
    from sklearn.metrics import adjusted_rand_score
    if not args.include_sensitivity or args.discovery_dir is None:raise ValueError('Cohort transfer requires all discovery cohorts and a primary discovery result')
    frozen=json.loads((args.discovery_dir/'candidate_freeze_proposal.json').read_text())
    artifact=args.discovery_dir/'frozen_transform_proposal.joblib'
    if sha256_file(artifact)!=frozen['transform_sha256']:raise ValueError('Primary transform hash mismatch')
    manifests=read_csv(args.embeddings/'embedding_manifest.csv')
    if len(manifests)!=len(records) or {m['internal_id'] for m in manifests}!={r.internal_id for r in records}:raise ValueError('Sensitivity embedding identities mismatch')
    embeddings={}
    for m in manifests:
        if m['split']!='discovery' or sha256_file(Path(m['embedding_path']))!=m['embedding_sha256']:raise ValueError('Sensitivity hash/split mismatch')
        with np.load(m['embedding_path'],allow_pickle=False) as z:embeddings[m['internal_id']]=z['residue_embedding']
    mapping=checked_coordinates(args.coordinates,records)
    rows,vectors,_=window_rows(records,embeddings,mapping,30,5)
    model=joblib.load(artifact)
    labels,strength=hdbscan.approximate_predict(model['clusterer'],model['pca'].transform(vectors))
    result=[]
    for name in sorted({r.analysis_cohort for r in records}):
        index=[i for i,row in enumerate(rows) if row['analysis_cohort']==name]
        counts=Counter(int(labels[i]) for i in index)
        result.append({'analysis_cohort':name,'n_sequences':sum(r.analysis_cohort==name for r in records),
                       'n_windows':len(index),'noise_fraction':counts[-1]/len(index),
                       'mean_mapped_fraction':float(np.mean([rows[i]['mapped_fraction'] for i in index])),
                       'cluster_counts':json.dumps(dict(sorted(counts.items()))),
                       'interpretation':'fixed primary transform transfer on discovery only; no refitting or functional labels'})
    write_csv(args.output_dir/'cohort_sensitivity.csv',result)
    original={(r['internal_id'],int(r['raw_start'])):int(r['cluster']) for r in read_csv(args.discovery_dir/'candidate_windows.csv')}
    primary_index=[i for i,row in enumerate(rows) if row['analysis_cohort']=='primary']
    ari=float(adjusted_rand_score([original[(rows[i]['internal_id'],rows[i]['raw_start'])] for i in primary_index],[int(labels[i]) for i in primary_index]))
    return {'n_sequences':len(records),'n_windows':len(rows),'cohorts':result,'validation_used':False,'refit':False,
            'primary_fit_vs_approximate_assignment_ARI':ari,
            'primary_transform_sha256':sha256_file(artifact),
            'software':{n:importlib.metadata.version(n) for n in ['numpy','scikit-learn','hdbscan']}}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['embed','analyze','cohort']);p.add_argument('--config',type=Path,default=ROOT/'configs/m3_discovery.json')
    for name in ['metadata','split-manifest','coordinates','output-dir']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--model-dir',type=Path);p.add_argument('--embeddings',type=Path);p.add_argument('--limit',type=int)
    p.add_argument('--batch-size',type=int,default=8);p.add_argument('--threads',type=int,default=4)
    p.add_argument('--include-sensitivity',action='store_true')
    p.add_argument('--discovery-dir',type=Path)
    args=p.parse_args();config=json.loads(args.config.read_text());args.output_dir.mkdir(parents=True,exist_ok=True)
    if (args.output_dir/'run_receipt.json').exists():raise ValueError('Use a new experiment output directory')
    provenance=clean_provenance(args.config);started=time.monotonic()
    try:
        if config['fit_split']!='discovery' or not config['primary_only']:raise ValueError('Discovery primary only')
        if args.batch_size<1 or args.threads<1 or (args.limit is not None and args.limit<1):raise ValueError('Positive counts required')
        records=checked_inputs(args,config)
        if args.stage=='analyze' and args.include_sensitivity:raise ValueError('Never refit primary discovery with sensitivity cohorts')
        result=embed(args,config,records) if args.stage=='embed' else cohort(args,config,records) if args.stage=='cohort' else analyze(args,config,records)
        files={str(path.relative_to(args.output_dir)):sha256_file(path) for path in args.output_dir.rglob('*') if path.is_file()}
        write_json(args.output_dir/'run_receipt.json',{'status':'completed','stage':args.stage,'provenance':provenance,
                  'config':config,'result':result,'seconds':time.monotonic()-started,'outputs_sha256':files})
    except Exception as exc:
        write_json(args.output_dir/'failed_run.json',{'status':'failed','stage':args.stage,'provenance':provenance,
                                                    'error':repr(exc),'seconds':time.monotonic()-started})
        raise

if __name__=='__main__':main()
