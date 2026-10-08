#!/usr/bin/env python3
"""Validate a specifically A-approved proposal without refitting parameters."""
import argparse
import hashlib
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from features.esm2_embed import read_embedding_records,iter_windows,sha256_file
from scripts.m3_discovery import clean_provenance,read_csv,write_csv,write_json


def approved_artifacts(discovery,approval_path):
    proposal=discovery/'candidate_freeze_proposal.json'
    receipt=discovery/'run_receipt.json'
    approval=json.loads(approval_path.read_text(encoding='utf-8'))
    if approval.get('approved') is not True or not approval.get('approved_by') or not approval.get('approval_evidence'):
        raise ValueError('A approval with genuine identity/evidence is required')
    if approval.get('proposal_sha256')!=sha256_file(proposal) or approval.get('discovery_receipt_sha256')!=sha256_file(receipt):
        raise ValueError('Approval must bind this exact proposal and discovery receipt')
    frozen=json.loads(proposal.read_text());run=json.loads(receipt.read_text())
    if run['status']!='completed' or run['result']['validation_used'] or frozen['fit_split']!='discovery':
        raise ValueError('Invalid discovery-only freeze source')
    for name,key in [('candidate_regions.csv','candidate_table_sha256'),('frozen_transform_proposal.joblib','transform_sha256')]:
        if sha256_file(discovery/name)!=frozen[key]:raise ValueError('Frozen artifacts changed after approval')
    selected=approval.get('selected_clusters')
    available={int(row['cluster']) for row in read_csv(discovery/'candidate_regions.csv')}
    if not isinstance(selected,list) or not selected or not set(selected)<=available:
        raise ValueError('A must explicitly select existing clusters for frozen validation')
    return approval,frozen,run


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['discovery-dir','approval','metadata','split-manifest','coordinates','model-dir','output-dir']:
        p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--batch-size',type=int,default=8);args=p.parse_args()
    # Approval checked BEFORE reading validation sequences or loading a model.
    approval,frozen,discovery=approved_artifacts(args.discovery_dir,args.approval)
    cfg=discovery['config'];provenance=clean_provenance(ROOT/'configs/m3_discovery.json')
    if provenance['config_sha256']!=discovery['provenance']['config_sha256']:raise ValueError('Discovery config changed')
    for path,key in [(args.metadata,'metadata_sha256'),(args.split_manifest,'split_manifest_sha256'),(args.coordinates,'coordinates_sha256')]:
        if sha256_file(path)!=cfg[key]:raise ValueError('Frozen validation input mismatch')
    records=read_embedding_records(args.metadata,args.split_manifest,selected_split='validation')
    if len(records)!=519:raise ValueError('Expected 519 primary validation sequences')
    if args.output_dir.exists():raise ValueError('Choose a fresh validation output directory')
    args.output_dir.mkdir(parents=True)
    try:
        import numpy as np
        import torch,joblib,hdbscan
        from transformers import AutoTokenizer,AutoModel
        model_receipt=json.loads((args.model_dir/'download_receipt.json').read_text())
        if model_receipt['revision']!=cfg['model_revision'] or model_receipt['model']!=cfg['model_id']:raise ValueError('Model revision mismatch')
        for item in discovery['result'].get('model_files',[]):
            if sha256_file(args.model_dir/item['name'])!=item['sha256']:raise ValueError('Model hash mismatch')
        # Discovery analysis binds model files through its embedding input receipt;
        # additionally verify every locally downloaded weight/tokenizer file here.
        for item in model_receipt['files']:
            if sha256_file(args.model_dir/item['name'])!=item['sha256']:raise ValueError('Model file mismatch')
        torch.set_num_threads(4);torch.manual_seed(cfg['seed']);torch.use_deterministic_algorithms(True)
        tokenizer=AutoTokenizer.from_pretrained(str(args.model_dir),local_files_only=True)
        model=AutoModel.from_pretrained(str(args.model_dir),local_files_only=True,use_safetensors=True,add_pooling_layer=False).eval()
        by_id={r.internal_id:r for r in records};mapping={k:{} for k in by_id}
        for row in read_csv(args.coordinates):
            if row['internal_id'] not in by_id:continue
            r=by_id[row['internal_id']];pos=int(row['raw_position'])
            if row['split']!='validation' or row['sequence_sha256']!=r.sequence_sha256 or not 1<=pos<=len(r.sequence) or r.sequence[pos-1]!=row['residue']:
                raise ValueError('Validation coordinate mismatch')
            if pos in mapping[r.internal_id]:raise ValueError('Ambiguous coordinate')
            mapping[r.internal_id][pos]=int(row['hmm_match_state']) if row['hmm_match_state'] else None
        rows=[];vectors=[]
        for offset in range(0,len(records),args.batch_size):
            batch=records[offset:offset+args.batch_size]
            encoded=tokenizer([r.sequence for r in batch],return_tensors='pt',padding=True,truncation=False)
            for i,r in enumerate(batch):
                if ''.join(tokenizer.convert_ids_to_tokens(encoded['input_ids'][i].tolist())[1:len(r.sequence)+1])!=r.sequence:
                    raise ValueError('Token/residue mismatch')
            with torch.inference_mode():output=model(**encoded).last_hidden_state
            for i,r in enumerate(batch):
                matrix=output[i,1:len(r.sequence)+1].cpu().numpy().astype(np.float32)
                if matrix.shape!=(len(r.sequence),480) or not np.isfinite(matrix).all():raise ValueError('Invalid validation embedding')
                for start,end in iter_windows(len(r.sequence),frozen['parameters']['window_length'],cfg['window_step']):
                    states=sorted({mapping[r.internal_id].get(pos) for pos in range(start,end+1)}-{None})
                    rows.append({'internal_id':r.internal_id,'split':'validation','raw_start':start,'raw_end':end,'hmm_states':';'.join(map(str,states))})
                    vectors.append(matrix[start-1:end].mean(axis=0))
        fitted=joblib.load(args.discovery_dir/'frozen_transform_proposal.joblib')
        transformed=fitted['pca'].transform(np.stack(vectors))
        labels,strength=hdbscan.approximate_predict(fitted['clusterer'],transformed)
        write_csv(args.output_dir/'validation_windows.csv',[dict(row,cluster=int(label),strength=float(score)) for row,label,score in zip(rows,labels,strength)])
        comparisons=[]
        for candidate in read_csv(args.discovery_dir/'candidate_regions.csv'):
            if int(candidate['cluster']) not in approval['selected_clusters']:continue
            cluster=int(candidate['cluster']);indices=np.flatnonzero(labels==cluster)
            start=candidate['majority_hmm_start'];end=candidate['majority_hmm_end'];fraction=[]
            if start and end:
                target=set(range(int(start),int(end)+1))
                for i in indices:
                    states={int(s) for s in rows[i]['hmm_states'].split(';') if s}
                    fraction.append(len(states&target)/max(len(states),1))
            comparisons.append({'cluster':cluster,'assigned_windows':len(indices),
                'validation_sequences':len({rows[i]['internal_id'] for i in indices}),
                'frozen_hmm_start':start,'frozen_hmm_end':end,
                'mean_fraction_mapped_states_in_frozen_region':float(np.mean(fraction)) if fraction else '',
                'interpretation':'fixed transfer coverage, not functional accuracy'})
        write_csv(args.output_dir/'validation_region_summary.csv',comparisons)
        write_json(args.output_dir/'run_receipt.json',{'status':'completed','provenance':provenance,
                   'A_approval':approval,'n_sequences':len(records),'n_windows':len(rows),'fit_performed':False,
                   'noise_fraction':float(np.mean(labels==-1)),
                   'outputs_sha256':{p.name:sha256_file(p) for p in args.output_dir.iterdir() if p.is_file()}})
    except Exception as exc:
        write_json(args.output_dir/'failed_run.json',{'status':'failed','provenance':provenance,'error':repr(exc)})
        raise

if __name__=='__main__':main()
