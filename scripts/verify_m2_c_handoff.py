"""Verify the downloaded C handoff bytes and recompute its scientific counts."""
import csv, hashlib, io, json, subprocess, zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
if not __debug__: raise RuntimeError('Do not disable verification assertions with python -O')
archive = ROOT / 'data/processed/m2_handoff_received/MLP_C_M2_scan_002_handoff.zip'
manifest = json.loads((ROOT / 'reports/M2/c_scan_002_handoff.json').read_text(encoding='utf-8'))
sha = lambda b: hashlib.sha256(b).hexdigest()
def rows(b): return list(csv.DictReader(io.StringIO(b.decode('utf-8'))))
assert sha(archive.read_bytes()) == manifest['archive_sha256']
assert len(archive.read_bytes()) == manifest['archive_bytes']
checks=[]
with zipfile.ZipFile(archive) as z:
    assert len(z.namelist()) == len(set(z.namelist()))
    listed = {item['name'] for item in manifest['files']}
    assert set(z.namelist()) == listed | {'manifest.json'}
    for item in manifest['files']:
        b=z.read(item['name'])
        assert len(b)==item['bytes'] and sha(b)==item['sha256'], item['name']
        checks.append(dict(name=item['name'],bytes=len(b),actual_sha256=sha(b),expected_sha256=item['sha256'],status='match'))
    internal=json.loads(z.read('manifest.json'))
    assert internal['files']==manifest['files']
    scan=json.loads(z.read('C_scan_002/scan_summary.json'))
    coordinate=json.loads(z.read('C_scan_002/coordinate_summary.json'))
    assert scan == json.loads((ROOT/'results/bioinformatics/pf00741_scan_002_summary.json').read_text(encoding='utf-8'))
    assert coordinate == json.loads((ROOT/'results/bioinformatics/hmm_coordinate_002_summary.json').read_text(encoding='utf-8'))
    assert scan['git']['dirty'] is False and scan['git']['commit']==manifest['run_commit']
    for name,h in scan['raw_output_sha256'].items(): assert sha(z.read('C_scan_002/raw/'+name))==h
    assert sha(z.read('reference/PF00741.hmm'))==scan['hmm']['sha256']
    assert sha(z.read('reference/PF00741_download_receipt.json'))==scan['hmm_download_receipt_sha256']
    with zipfile.ZipFile(io.BytesIO(z.read('B_input/gvpa_v1_reproduction.zip'))) as bz:
        bfiles = {'metadata.csv':'metadata_sha256','sequences_for_clustering.fasta':'fasta_sha256',
                  'split/split_manifest.csv':'split_manifest_sha256'}
        for name,key in bfiles.items():
            b=bz.read('gvpa_v1_reproduction/'+name)
            assert sha(b)==scan['input'][key], name
            checks.append(dict(name='B_input/'+name,actual_sha256=sha(b),expected_sha256=scan['input'][key],status='match'))
        # Reproduction summaries have different source provenance; do not call them scan inputs.
        for name,key,path in [('audit_summary.json','audit_summary_sha256','results/data_audit/gvpa_v1_audit_summary.json'),('split/split_summary.json','split_summary_sha256','results/data_audit/gvpa_v1_split_summary.json')]:
            packaged=bz.read('gvpa_v1_reproduction/'+name)
            original_blob=subprocess.check_output(['git','show',manifest['run_commit']+':'+path])
            windows_checkout=original_blob.replace(b'\r\n',b'\n').replace(b'\n',b'\r\n')
            assert sha(windows_checkout)==scan['input'][key]
            a=json.loads(packaged); b=json.loads(original_blob)
            relevant=['candidate_records','sequence_qc_eligible','cohort_counts','dataset_version','output_sha256'] if key.startswith('audit') else ['dataset_version','n_clusters','split_counts']
            for field in relevant:
                if field in a or field in b: assert a.get(field)==b.get(field), field
            checks.append(dict(name='B_input/'+name,actual_sha256=sha(packaged),scan_recorded_sha256=scan['input'][key],reconstructed_windows_checkout_sha256=sha(windows_checkout),packaged_split_version=a.get('split_version'),original_split_version=b.get('split_version'),status='reproduction_summary_differs_from_scan_input; core_statistics_match; original_CRLF_checkout_reconstructed_not_obtained'))
        metadata=rows(bz.read('gvpa_v1_reproduction/metadata.csv'))
        split=rows(bz.read('gvpa_v1_reproduction/split/split_manifest.csv'))
        audit=json.loads(bz.read('gvpa_v1_reproduction/audit_summary.json'))
        split_summary=json.loads(bz.read('gvpa_v1_reproduction/split/split_summary.json'))
        fasta={}
        for line in bz.read('gvpa_v1_reproduction/sequences_for_clustering.fasta').decode().splitlines():
            if line.startswith('>'):
                identity=line[1:].split()[0]
                assert identity not in fasta
                fasta[identity]=''
            elif line.strip(): fasta[identity]+=line.strip()
        assert len(metadata)==2078 and len({r['internal_id'] for r in metadata})==2078
        meta={r['internal_id']:r for r in metadata}
        eligible={k for k,v in meta.items() if v['sequence_qc_eligible'].lower()=='true'}
        assert len(eligible)==2076 and set(fasta)==eligible
        for k in eligible:
            assert fasta[k]==meta[k]['sequence'] and sha(fasta[k].encode())==meta[k]['sequence_sha256']
        split_by_id={r['internal_id']:r for r in split}
        assert len(split)==len(split_by_id)==2076 and set(split_by_id)==eligible
        clusters=defaultdict(set)
        for r in split: clusters[r['homology_cluster']].add(r['split'])
        assert len(clusters)==478 and all(len(v)==1 for v in clusters.values())
        split_counts=Counter(r['split'] for r in split)
        primary=Counter(r['split'] for r in split if meta[r['internal_id']]['primary_analysis_eligible'].lower()=='true')
        cohorts=Counter(r['analysis_cohort'] for r in metadata)
        assert dict(cohorts)==audit['cohort_counts']
        assert split_counts=={'discovery':1453,'validation':623} and primary=={'discovery':1202,'validation':519}
        # Materialize only checked inputs in a new ignored directory, preserving bytes.
        for name in bz.namelist():
            if name.endswith('/'): continue
            safe=PurePosixPath(name)
            assert not safe.is_absolute() and '..' not in safe.parts and ':' not in name
            dest=ROOT/'data/processed/m2_handoff_received/unpacked/B_input'/Path(name)
            dest.parent.mkdir(parents=True,exist_ok=True)
            b=bz.read(name)
            if dest.exists(): assert dest.read_bytes()==b
            else: dest.write_bytes(b)
    statuses=rows(z.read('C_scan_002/raw/per_sequence_status.csv'))
    status_by_id={r['internal_id']:r for r in statuses}
    assert len(statuses)==len(status_by_id)==2076 and set(status_by_id)==eligible
    assert Counter(r['status'] for r in statuses)==scan['counts']=={'accepted':2076}
    for r in statuses:
        k=r['internal_id']
        assert r['sequence_sha256']==meta[k]['sequence_sha256']
        assert r['split']==split_by_id[k]['split'] and r['analysis_cohort']==meta[k]['analysis_cohort']
    coords=rows(z.read('C_scan_002/raw/hmm_coordinate_map.csv'))
    assert len(coords)==coordinate['n_coordinate_rows']==80219
    coord_by_id=defaultdict(list)
    seen=set()
    for r in coords:
        k=r['internal_id']; pos=int(r['raw_position']); state=r['hmm_match_state']
        assert k in eligible and meta[k]['sequence'][pos-1]==r['residue']
        assert r['split']==split_by_id[k]['split'] and r['sequence_sha256']==meta[k]['sequence_sha256']
        if state: assert 1<=int(state)<=39
        t=(k,r['domain_number'],pos); assert t not in seen; seen.add(t)
        coord_by_id[k].append(r)
    assert set(coord_by_id)==eligible
    for k,rs in coord_by_id.items():
        assert len(rs)==int(status_by_id[k]['n_mapped_residues'])
        assert len({r['hmm_match_state'] for r in rs if r['hmm_match_state']})==int(status_by_id[k]['n_match_states'])
    raw=[json.loads(line) for line in z.read('C_scan_002/raw/raw_hits.jsonl').splitlines()]
    assert all(r['internal_id'] in eligible for r in raw)
    gathering=[r for r in raw if r['search']=='gathering']
    assert len(gathering)==2076 and {r['internal_id'] for r in gathering}==eligible
    assert all(r['hit'] and r['hit']['included'] for r in gathering)
    assert json.loads(z.read('C_scan_002/raw/failures.json'))['records']==[]
    structure=rows(z.read('reference/7r1c_residue_hmm_map.csv'))
    assert len(structure)==88 and sum(r['modeled']=='True' for r in structure)==65
    assert len({r['hmm_match_state'] for r in structure if r['hmm_match_state']})==39
    sr=json.loads((ROOT/'results/bioinformatics/structure_7r1c_summary.json').read_text(encoding='utf-8'))
    assert sha(z.read('reference/7R1C.pdb'))==sr['pdb_sha256']
    assert sha(z.read('reference/7r1c_residue_hmm_map.csv'))==sr['coordinate_map_sha256']
    for name in listed-{'B_input/gvpa_v1_reproduction.zip'}:
        dest=ROOT/'data/processed/m2_handoff_received/unpacked'/Path(name)
        dest.parent.mkdir(parents=True,exist_ok=True)
        b=z.read(name)
        if dest.exists(): assert dest.read_bytes()==b
        else: dest.write_bytes(b)
implementation=[]
for key,path in [('workflow','bioinformatics/m2_pf00741.py'),('coordinates','bioinformatics/m2_coordinates.py')]:
    b=subprocess.check_output(['git','show',manifest['run_commit']+':'+path])
    h=sha(b); assert h==scan['implementation_sha256'][key]
    assert b==subprocess.check_output(['git','show','HEAD:'+path])
    implementation.append(dict(path=path,recorded_commit=manifest['run_commit'],git_blob_sha256=h,status='match'))
out={'kind':'assistant_technical_verification_not_A_acceptance_or_member_signature',
     'checked_utc':datetime.now(timezone.utc).isoformat(),'checked_sha':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
     'shared_uri':manifest['shared_uri'],'access':'actually downloaded from public Release URL by assistant',
     'archive_sha256_actual':sha(archive.read_bytes()),'file_checks':checks,'implementation_checks':implementation,
     'recomputed':{'metadata':len(metadata),'qc_eligible':len(eligible),'cohorts':dict(cohorts),'split':dict(split_counts),
                   'primary_split':dict(primary),'clusters':len(clusters),'cross_split_clusters':0,'accepted':2076,
                   'coordinate_rows':len(coords),'coordinate_residue_identity_checked':'all rows',
                   'structure_deposited':88,'structure_modeled':65,'structure_hmm_states':39},
     'old_dirty_run':'historical receipt preserved; clean scan 002 provides replacement evidence, not reconstruction of old dirty environment',
     'pending':['A personal acceptance','team confirmation of public redistribution conditions','legacy E001 provenance remains unknown']}
dest=ROOT/'reports/M2/c_handoff_technical_verification.json'
dest.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('PASS: archive + 12 files + 3 exact B inputs; 2 reproduction summaries core statistics checked; all coordinate residues and split counts; clean-code provenance')
