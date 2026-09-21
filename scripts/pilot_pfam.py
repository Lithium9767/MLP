"""Local-only Pfam pilot. No network calls; no training; stdlib Python >=3.10."""
import argparse
import collections
import csv
import datetime as dt
import difflib
import hashlib
import json
import math
import random
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

GROUPS = ('clear_gvpa', 'conflicting_gvp', 'other_gvp')
AA = set('ACDEFGHIKLMNPQRSTVWY')
MANIFEST = 'internal_id original_accession source_file source_record_index organism sequence_length sample_group original_annotation partial_status sequence_sha256 selection_reason qc_status input_json source_types'.split()
QC = 'original_accession input_json source_record_index reason evidence related_accession'.split()
STATUS = 'internal_id submitted scan_status result_present parse_status failure_reason tool_name tool_version database_name database_version run_date'.split()
HITS = 'internal_id signature_database signature_accession signature_name interpro_accession start end score evalue domain_evalue sequence_length hit_length sequence_coverage hmm_coverage hit_status review_reason raw_result_reference'.split()
LABELS = 'internal_id sample_group scan_status pf00741_detected label label_reason needs_review sequence_sha256'.split()

def stamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def sha(data):
    return hashlib.sha256(data).hexdigest()

def write_json(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

def csv_write(path, fields, rows):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

def csv_read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def load_rules(path):
    rules=json.loads(Path(path).read_text(encoding='utf-8'))
    if rules.get('signature_database')!='Pfam' or rules.get('signature_accession')!='PF00741':
        raise ValueError('This pilot only accepts Pfam PF00741 rules')
    for key in ('max_sequence_evalue','max_domain_evalue','min_hit_length','min_hmm_coverage'):
        value=rules.get(key)
        if not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
            raise ValueError('Invalid frozen rule: '+key)
    if rules['min_hmm_coverage']>1:
        raise ValueError('HMM coverage must not exceed 1')
    return rules

def fresh(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise ValueError(f'Refusing to overwrite non-empty output: {path}; use a new run directory')
    path.mkdir(parents=True, exist_ok=True)
    return path

def input_bytes(spec):
    if '::' in spec:
        archive, member = spec.split('::', 1)
        with zipfile.ZipFile(archive) as z:
            return z.read(member)
    return Path(spec).read_bytes()

def ann_types(a):
    values = a.get('gvp_types', [])
    if isinstance(values, str):
        values = [values]
    text = ' '.join(str(v) for v in values) + ' ' + ' '.join(str(a.get(k, '')) for k in ('gene', 'description'))
    return set('Gvp' + m.upper() for m in re.findall(r'gvp([a-z])(?![a-z])', text, re.I))

def parse_input(spec, background=False):
    data = input_bytes(spec)
    raw = json.loads(data.decode('utf-8-sig'))
    if not isinstance(raw, list) or not raw:
        raise ValueError(f'{spec}: expected non-empty JSON list')
    rows, exclusions = [], []
    for i, x in enumerate(raw, 1):
        if not isinstance(x, dict) or 'sequence' not in x or 'representative_annotation' not in x:
            raise ValueError(f'{spec} record {i}: unsupported schema; sequence/representative_annotation required')
        a = x['representative_annotation']
        members = x.get('redundant_members_details', [])
        if not isinstance(a, dict) or not isinstance(members, list) or not all(isinstance(m, dict) for m in members):
            raise ValueError(f'{spec} record {i}: invalid annotation object/list')
        annotations = [a] + members
        seq = x['sequence']
        if not isinstance(seq, str):
            raise ValueError(f'{spec} record {i}: sequence must be string')
        # Preserve source exactly: inputs requiring normalization are quarantined.
        acc = x.get('unique_sequence_id') or a.get('id', '')
        all_types = set().union(*(ann_types(m) for m in annotations))
        partial = any(re.search(r'\bpartial\b', str(m.get('description', '')) + ' ' + str(m.get('full_header', '')), re.I) for m in annotations)
        reasons = []
        if not seq:
            reasons.append('empty_sequence')
        if set(seq) - AA:
            reasons.append('invalid_residues:' + repr(''.join(sorted(set(seq)-AA))))
        if len(seq) < 40 or len(seq) > 800:
            reasons.append('pilot_length_outside_40_800')
        if partial:
            reasons.append('partial_in_representative_or_member')
        if not all((acc, a.get('organism'), a.get('source_file'), all_types)):
            reasons.append('insufficient_annotation')
        if seq and max(collections.Counter(seq).values()) / len(seq) > 0.6:
            reasons.append('low_complexity_single_residue_over_60_percent')
        group = ('clear_gvpa' if all_types == {'GvpA'} else 'conflicting_gvp')
        if background:
            group = 'other_gvp'
            if 'GvpA' in all_types or len(all_types) != 1:
                reasons.append('background_not_unambiguous_nonA')
        row = dict(original_accession=acc, input_json=spec, source_record_index=i, sequence=seq,
                   source_file=a.get('source_file', ''), organism=a.get('organism', ''),
                   sequence_length=len(seq), sample_group=group, original_annotation=a.get('description', ''),
                   partial_status='detected' if partial else 'not_detected_not_proof_of_completeness',
                   sequence_sha256=sha(seq.encode()), source_types='|'.join(sorted(all_types)), qc_status='pass_pilot_screen')
        if reasons:
            exclusions.append({**row, 'reason': '|'.join(reasons), 'evidence': row['original_annotation']})
        else:
            rows.append(row)
    return rows, exclusions, dict(path=spec, sha256=sha(data), records=len(raw), top_level_fields=sorted(set().union(*(x.keys() for x in raw))), annotation_fields=sorted(set().union(*(x['representative_annotation'].keys() for x in raw))))

def near(a, b):
    """Conservative approximate global matching-block screen, not biological alignment."""
    if min(len(a), len(b)) / max(len(a), len(b)) < 0.8:
        return False
    if sum((collections.Counter(a) & collections.Counter(b)).values()) / min(len(a), len(b)) < 0.9:
        return False
    matched = sum(m.size for m in difflib.SequenceMatcher(None, a, b, autojunk=False).get_matching_blocks())
    return matched / min(len(a), len(b)) >= 0.9

def select_rows(rows, targets, seed):
    rng = random.Random(seed)
    rows = sorted(rows, key=lambda r: (r['input_json'], r['source_record_index']))
    rng.shuffle(rows)
    rank = {id(r): i for i,r in enumerate(rows)}
    seen, unique, excluded = {}, [], []
    # Prefer explicit primary source over background when exact sequences repeat.
    for r in sorted(rows, key=lambda r: (r['sample_group']=='other_gvp', rank[id(r)])):
        digest = r['sequence_sha256']
        if digest in seen:
            excluded.append({**r, 'reason':'exact_duplicate', 'related_accession':seen[digest]['original_accession'], 'evidence':digest})
        else:
            seen[digest] = r
            unique.append(r)
    selected = []
    for group, target in zip(GROUPS, targets):
        pool = [r for r in unique if r['sample_group']==group]
        org, source, bins, types = collections.Counter(), collections.Counter(), collections.Counter(), collections.Counter()
        while pool and sum(r['sample_group']==group for r in selected) < target:
            def priority(r):
                return (org[r['organism']], types[r['source_types']], bins[r['sequence_length']//25], source[r['source_file']], rank[id(r)])
            r = min(pool, key=priority)
            pool.remove(r)
            close = next((x for x in selected if near(r['sequence'],x['sequence'])), None)
            if close:
                excluded.append({**r, 'reason':'near_duplicate_of_selected_approximate', 'related_accession':close['original_accession'], 'evidence':'matching_blocks/min_length>=0.9;min_length/max_length>=0.8'})
                continue
            selected.append(r)
            org[r['organism']] += 1
            source[r['source_file']] += 1
            bins[r['sequence_length']//25] += 1
            types[r['source_types']] += 1
        excluded.extend({**r,'reason':'not_selected_quota_diversity','evidence':'near-duplicate status not exhaustively assessed'} for r in pool)
    for i,r in enumerate(selected,1):
        r['internal_id'] = f'PILOT_{i:04d}'
        r['selection_reason'] = 'seeded diversity across organism/type/25aa bins/source; passed approximate near-duplicate screen'
    return selected, excluded

def select(args):
    inputs = [args.input_json] + args.background_json
    if len(inputs) != len(set(inputs)):
        raise ValueError('Duplicate input specifications')
    all_rows, qc, inspected = [], [], []
    for i,spec in enumerate(inputs):
        rows, errors, info = parse_input(spec, i>0)
        all_rows.extend(rows); qc.extend(errors); inspected.append(info)
    chosen, rejected = select_rows(all_rows, args.targets, args.seed)
    if not chosen:
        raise ValueError('No usable sequences selected')
    qc.extend(rejected)
    out = fresh(args.output_dir)
    csv_write(out/'pilot_manifest.csv', MANIFEST, chosen)
    csv_write(out/'excluded_or_qc_sequences.csv', QC, qc)
    (out/'pilot_sequences.fasta').write_text(''.join('>'+r['internal_id']+'\n'+r['sequence']+'\n' for r in chosen),encoding='utf-8')
    try:
        commit = subprocess.check_output(['git','rev-parse','HEAD'], cwd=Path(__file__).resolve().parents[1], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = ''
    summary = dict(inputs=inspected, run_date=stamp(), seed=args.seed, target_counts=dict(zip(GROUPS,args.targets)),
                   actual_counts=dict(collections.Counter(r['sample_group'] for r in chosen)),
                   passed_basic_qc=len(all_rows), qc_or_unselected_count=len(qc),
                   exclusion_reasons=dict(collections.Counter(r['reason'] for r in qc)),
                   script_sha256=sha(Path(__file__).read_bytes()), git_commit=commit,
                   rules='40..800 aa pilot screen; canonical uppercase residues unchanged; any partial excluded; exact dedup; seeded organism/type/length-bin/source diversity; approximate matching-block near screen',
                   limitations='Not a random prevalence sample. Existing CD-HIT thresholds lack provenance; no verified similarity clusters reused. Near screen is approximate, not a full biological clustering. No partial text does not prove completeness.',
                   selected_details={g:dict(organisms=len({r['organism'] for r in chosen if r['sample_group']==g}), lengths=sorted(r['sequence_length'] for r in chosen if r['sample_group']==g),types=dict(collections.Counter(r['source_types'] for r in chosen if r['sample_group']==g))) for g in GROUPS})
    write_json(out/'selection_summary.json',summary)
    (out/'selection.log').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    validate_pilot(out)
    print(json.dumps(summary['actual_counts']))

def validate_pilot(path):
    path = Path(path)
    rows = csv_read(path/'pilot_manifest.csv')
    seqs = {}
    current = None
    for line in (path/'pilot_sequences.fasta').read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            current = line[1:]
            if current in seqs or not re.fullmatch(r'PILOT_\d{4,}',current):
                raise ValueError('Duplicate or invalid FASTA ID')
            seqs[current] = ''
        elif line:
            if current is None:
                raise ValueError('FASTA sequence before header')
            seqs[current] += line
    if len(rows)!=len(seqs) or len({r['internal_id'] for r in rows})!=len(rows) or set(seqs)!={r['internal_id'] for r in rows}:
        raise ValueError('Manifest/FASTA ID mismatch')
    hashes = set()
    for r in rows:
        s=seqs[r['internal_id']]
        if not s or set(s)-AA or len(s)!=int(r['sequence_length']) or sha(s.encode())!=r['sequence_sha256']:
            raise ValueError('FASTA hash, length or alphabet mismatch: '+r['internal_id'])
        if r['sequence_sha256'] in hashes:
            raise ValueError('Duplicate selected sequence')
        hashes.add(r['sequence_sha256'])
    return rows,seqs

def run(args):
    pilot=Path(args.pilot_dir).resolve()
    rows,seqs=validate_pilot(pilot)
    rule_path=Path(args.rules).resolve()
    rules=load_rules(rule_path)
    out=fresh(args.output_dir).resolve()
    for name in ('raw','logs','parsed','reports'):
        (out/name).mkdir()
    shutil.copyfile(rule_path,out/'rules.json')
    receipt=dict(run_date=stamp(),scan_status='scan_not_run',submitted=False,exit_code=None,
                 tool_name='InterProScan',tool_version='',database_name='Pfam',database_version='',
                 fasta_sha256=sha((pilot/'pilot_sequences.fasta').read_bytes()),
                 manifest_sha256=sha((pilot/'pilot_manifest.csv').read_bytes()),
                 rules_sha256=sha(rule_path.read_bytes()), ids=[r['internal_id'] for r in rows],command=[],
                 script_sha256=sha(Path(__file__).read_bytes()),
                 reason='Local scanner unavailable or prepare-only; remote submission not authorized')
    exe=shutil.which(args.interproscan) if args.interproscan else None
    if args.interproscan and Path(args.interproscan).is_file():
        exe=str(Path(args.interproscan).resolve())
    if args.execute and exe:
        exe=str(Path(exe).resolve())
        try:
            help_run=subprocess.run([exe,'--help'],capture_output=True,text=True,timeout=120,cwd=str(Path(exe).parent))
            help_text=help_run.stdout+'\n'+help_run.stderr
            (out/'logs/help.log').write_text(help_text,encoding='utf-8')
            version_run=subprocess.run([exe,'--version'],capture_output=True,text=True,timeout=120,cwd=str(Path(exe).parent))
            version_text=version_run.stdout+'\n'+version_run.stderr
            (out/'logs/version.log').write_text(version_text,encoding='utf-8')
            tool=re.search(r'(5\.\d+-\d+\.\d+)',version_text+'\n'+help_text)
            pfam=re.search(r'Pfam[-\s]+(\d+\.\d+)',help_text,re.I)
            if not tool or not pfam or not all(s in help_text for s in ('--disable-precalc','--applications','--formats','--output-file-base')):
                raise ValueError('Cannot verify InterProScan5 version/Pfam version/required flags from installed help; no scan submitted')
            receipt.update(tool_version=tool.group(1),database_version=pfam.group(1))
            command=[exe,'-i',str(pilot/'pilot_sequences.fasta'),'-appl','Pfam','-f','TSV,XML','-b',str(out/'raw/pilot'),'-dp']
            receipt.update(submitted=True,command=command,scan_status='running',reason='')
            write_json(out/'run_receipt.json',receipt)
            with (out/'logs/stdout.log').open('w',encoding='utf-8') as stdout,(out/'logs/stderr.log').open('w',encoding='utf-8') as stderr:
                completed=subprocess.run(command,cwd=str(Path(exe).parent),stdout=stdout,stderr=stderr)
            receipt.update(exit_code=completed.returncode,scan_status='completed' if completed.returncode==0 else 'failed')
            xml=out/'raw/pilot.xml'
            if xml.exists():
                receipt['xml_sha256']=sha(xml.read_bytes())
        except (OSError,ValueError,subprocess.TimeoutExpired) as e:
            receipt['reason']=str(e)
            receipt['scan_status']='failed' if receipt['submitted'] else 'scan_not_run'
    elif args.execute:
        receipt['reason']='Requested local executable not found; no sequence submitted'
    write_json(out/'run_receipt.json',receipt)
    (out/'logs/preflight.log').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    # Empty header-only hits, blank pending labels: never fabricated scanner results.
    initial=[]
    for r in rows:
        initial.append(dict(internal_id=r['internal_id'],submitted=receipt['submitted'],scan_status='not_submitted' if not receipt['submitted'] else 'result_missing',result_present=False,parse_status='not_run',failure_reason=receipt['reason'],tool_name='InterProScan',tool_version=receipt['tool_version'],database_name='Pfam',database_version=receipt['database_version'],run_date=receipt['run_date']))
    csv_write(out/'parsed/pilot_scan_status.csv',STATUS,initial)
    csv_write(out/'parsed/pilot_pfam_hits.csv',HITS,[])
    csv_write(out/'parsed/pilot_labels.csv',LABELS,[dict(internal_id=r['internal_id'],sample_group=r['sample_group'],scan_status=s['scan_status'],pf00741_detected='',label='',label_reason='scan_not_run' if not receipt['submitted'] else 'awaiting_verified_parse',needs_review=True,sequence_sha256=r['sequence_sha256']) for r,s in zip(rows,initial)])
    print(receipt['scan_status'])

def local_name(tag):
    return tag.rsplit('}',1)[-1]

def child(el,name):
    return next((x for x in el if local_name(x.tag)==name),None)

def numeric(value):
    if value is None or value=='':
        return None
    number=float(value)
    if not math.isfinite(number) or number<0:
        raise ValueError('Invalid nonnegative numeric field')
    return number

def parse_protein(protein,row,expected_seq,rules,db_version,reference):
    seqnode=child(protein,'sequence')
    seq=''.join(seqnode.itertext()).strip() if seqnode is not None else ''
    if seq!=expected_seq or (seqnode.get('md5') and seqnode.get('md5').lower()!=hashlib.md5(seq.encode()).hexdigest()):
        raise ValueError('XML sequence or MD5 differs from submitted input')
    matches=child(protein,'matches')
    if matches is None:
        raise ValueError('Missing matches container; completeness unknown')
    hits=[]
    for match in matches:
        signature=child(match,'signature')
        if signature is None:
            raise ValueError('Missing signature')
        release=child(signature,'signature-library-release')
        if release is None or release.get('library','').upper()!='PFAM' or release.get('version')!=db_version:
            raise ValueError('Unexpected or mixed signature database/version')
        if signature.get('ac')!=rules['signature_accession']:
            continue
        locations=child(match,'locations')
        if locations is None or len(locations)==0:
            raise ValueError('Target signature has no locations')
        entry=child(signature,'entry')
        for loc in locations:
            start,end=int(loc.get('start')),int(loc.get('end'))
            if not 1<=start<=end<=len(seq):
                raise ValueError('Invalid residue coordinates')
            ev=numeric(match.get('evalue'))
            dev=numeric(loc.get('evalue'))
            reasons=[]
            if ev is None or dev is None:
                reasons.append('significance_missing')
            elif ev>=rules['max_sequence_evalue'] or dev>=rules['max_domain_evalue']:
                reasons.append('weak_or_boundary_evalue')
            hmm_cov=None
            if all(loc.get(k) is not None for k in ('hmm-start','hmm-end','hmm-length')):
                hs,he,hl=(int(loc.get(k)) for k in ('hmm-start','hmm-end','hmm-length'))
                if not 1<=hs<=he<=hl:
                    raise ValueError('Invalid HMM coordinates')
                hmm_cov=(he-hs+1)/hl
            if hmm_cov is None:
                reasons.append('hmm_coverage_missing')
            elif hmm_cov<rules['min_hmm_coverage']:
                reasons.append('low_hmm_coverage')
            if end-start+1<rules['min_hit_length']:
                reasons.append('short_hit')
            hits.append(dict(internal_id=row['internal_id'],signature_database='Pfam',signature_accession=signature.get('ac'),signature_name=signature.get('name',''),interpro_accession=entry.get('ac','') if entry is not None else '',start=start,end=end,score=match.get('score',''),evalue=match.get('evalue',''),domain_evalue=loc.get('evalue',''),sequence_length=len(seq),hit_length=end-start+1,sequence_coverage=(end-start+1)/len(seq),hmm_coverage='' if hmm_cov is None else hmm_cov,hit_status='needs_review' if reasons else 'accepted',review_reason='|'.join(reasons),raw_result_reference=reference))
    return hits

def parse(args):
    pilot=Path(args.pilot_dir)
    rows,seqs=validate_pilot(pilot)
    run_dir=Path(args.run_dir)
    receipt=json.loads((run_dir/'run_receipt.json').read_text(encoding='utf-8'))
    rules=load_rules(run_dir/'rules.json')
    for key,path in [('fasta_sha256',pilot/'pilot_sequences.fasta'),('manifest_sha256',pilot/'pilot_manifest.csv'),('rules_sha256',run_dir/'rules.json')]:
        if receipt.get(key)!=sha(path.read_bytes()):
            raise ValueError('Receipt input/config hash mismatch: '+key)
    out=fresh(args.output_dir)
    xml=run_dir/'raw/pilot.xml'
    elements={}
    global_error=''
    if receipt.get('submitted') is not True:
        global_error='not_submitted'
    elif receipt.get('exit_code')!=0 or receipt.get('scan_status')!='completed':
        global_error='failed'
    elif not xml.exists():
        global_error='result_missing'
    else:
        try:
            if receipt.get('xml_sha256')!=sha(xml.read_bytes()):
                raise ValueError('Raw XML hash mismatch')
            root=ET.parse(xml).getroot()
            if root.get('interproscan-version')!=receipt.get('tool_version') or not receipt.get('database_version'):
                raise ValueError('Missing or mismatched tool/database version')
            if set(receipt.get('ids',[]))!=set(seqs):
                raise ValueError('Receipt ID list mismatch')
            for p in root:
                if local_name(p.tag)!='protein':
                    raise ValueError('Unexpected XML record type')
                xrefs=[x.get('id') for x in p if local_name(x.tag)=='xref']
                if not xrefs:
                    raise ValueError('Protein missing ID')
                for xid in xrefs:
                    if xid not in seqs or xid in elements:
                        raise ValueError('Unknown/duplicate XML ID')
                    elements[xid]=p
        except (ValueError,ET.ParseError,OSError) as e:
            global_error='parse_failed'
            (out/'parse_error.log').write_text(str(e),encoding='utf-8')
    statuses,labels,hits=[],[],[]
    for row in rows:
        iid=row['internal_id'];state=global_error;reason=global_error;these=[]
        if not state:
            if iid not in elements:
                state=reason='result_missing'
            else:
                try:
                    these=parse_protein(elements[iid],row,seqs[iid],rules,receipt['database_version'],str(xml))
                    state='needs_review' if any(h['hit_status']=='needs_review' for h in these) else 'success'
                    reason='uncertain_target_match' if state=='needs_review' else ''
                    if not these and (row['qc_status']!='pass_pilot_screen' or row['partial_status']!='not_detected_not_proof_of_completeness'):
                        state='needs_review';reason='input_not_eligible_for_negative'
                except (ValueError,TypeError,KeyError) as e:
                    state='parse_failed';reason=str(e)
        hits.extend(these)
        label=''
        if state=='success':
            label='1' if these else '0'
            reason='PF00741 accepted under frozen computational rules' if these else '在当前工具、数据库版本和固定计算规则下未检出合格 PF00741 命中'
        statuses.append(dict(internal_id=iid,submitted=receipt.get('submitted',False),scan_status=state,result_present=iid in elements,parse_status='success' if state in ('success','needs_review') else state,failure_reason='' if state=='success' else reason,tool_name='InterProScan',tool_version=receipt.get('tool_version',''),database_name='Pfam',database_version=receipt.get('database_version',''),run_date=receipt.get('run_date','')))
        labels.append(dict(internal_id=iid,sample_group=row['sample_group'],scan_status=state,pf00741_detected=label,label=label,label_reason=reason,needs_review=state!='success',sequence_sha256=row['sequence_sha256']))
    csv_write(out/'pilot_scan_status.csv',STATUS,statuses)
    csv_write(out/'pilot_pfam_hits.csv',HITS,hits)
    csv_write(out/'pilot_labels.csv',LABELS,labels)
    write_json(out/'parse_provenance.json',dict(run_receipt_sha256=sha((run_dir/'run_receipt.json').read_bytes()),rules_sha256=receipt['rules_sha256'],script_sha256=sha(Path(__file__).read_bytes()),run_date=stamp()))
    print(dict(collections.Counter(r['scan_status'] for r in statuses)))

def summarize(args):
    parsed=Path(args.parsed_dir)
    labels=csv_read(parsed/'pilot_labels.csv'); statuses=csv_read(parsed/'pilot_scan_status.csv')
    sm={r['internal_id']:r for r in statuses}
    if len(sm)!=len(statuses) or len({r['internal_id'] for r in labels})!=len(labels) or set(sm)!={r['internal_id'] for r in labels}:
        raise ValueError('Status/label ID mismatch')
    for r in labels:
        s=sm[r['internal_id']]
        if r['scan_status']!=s['scan_status'] or r['label']!=r['pf00741_detected']:
            raise ValueError('Status/label inconsistency')
        if r['label'] not in ('','0','1') or (r['label']!='' and (s['scan_status']!='success' or s['submitted']!='True' or s['result_present']!='True' or s['parse_status']!='success')):
            raise ValueError('Label without verified successful result')
        if s['scan_status']=='success' and r['label']=='':
            raise ValueError('Successful classification cannot have blank label')
    table=[]
    for group in GROUPS:
        selected=[r for r in labels if r['sample_group']==group]
        submitted=[r for r in selected if sm[r['internal_id']]['submitted']=='True']
        pos=sum(r['label']=='1' for r in submitted);neg=sum(r['label']=='0' for r in submitted)
        pending=sum(r['label']=='' for r in submitted)
        row=dict(sample_group=group,selected=len(selected),submitted=len(submitted),scan_success=pos+neg,pf00741_positive=pos,not_detected=neg,pending_or_failed=pending,not_submitted=len(selected)-len(submitted),needs_review=sum(r['scan_status']=='needs_review' for r in submitted))
        assert row['submitted']==row['scan_success']+pending
        assert row['scan_success']==pos+neg
        table.append(row)
    out=fresh(args.output_dir)
    csv_write(out/'pilot_group_summary.csv',list(table[0]),table)
    total=sum(r['submitted'] for r in table)
    decision='scan_not_run：等待真实扫描结果；不能判断正负类别可行性。' if not total else '依据成功、失败和待复核数量进一步评估；正负标签出现不等于二分类任务成立。'
    lines=['# Pfam pilot report','',decision,'','| group | selected | submitted | success | positive | not detected | pending/failed | not submitted | review |','| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    lines.extend('| '+' | '.join(str(v) for v in r.values())+' |' for r in table)
    lines.extend(['','统计口径：success 为已成功扫描且可完成标签判定；完成计算但边界命中待复核的序列计入 pending/failed 和 review，不计入 success。not_submitted 单列，不放入已提交分母。','', '本次试扫描采用刻意分层选样，不代表原始数据分布，因此不能用于估计全量数据的 PF00741 阳性比例。','', '标签为空不等于阴性。原始命中不存在时不得报告“全部未检出”。工具和数据库版本见 run_receipt.json，固定阈值见 rules.json。','', '当前不得以本批次支持全量训练或宣称模型可行性。没有训练模型或上传序列。'])
    (out/'pilot_scan_report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(decision)

def cli():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='cmd',required=True)
    s=sub.add_parser('select');s.add_argument('--input-json',required=True);s.add_argument('--background-json',action='append',default=[]);s.add_argument('--output-dir',required=True);s.add_argument('--seed',type=int,default=20260908);s.add_argument('--targets',type=int,nargs=3,default=[30,15,15]);s.set_defaults(func=select)
    r=sub.add_parser('run');r.add_argument('--pilot-dir',required=True);r.add_argument('--output-dir',required=True);r.add_argument('--rules',required=True);r.add_argument('--interproscan');r.add_argument('--execute',action='store_true');r.set_defaults(func=run)
    a=sub.add_parser('parse');a.add_argument('--pilot-dir',required=True);a.add_argument('--run-dir',required=True);a.add_argument('--output-dir',required=True);a.set_defaults(func=parse)
    s=sub.add_parser('summarize');s.add_argument('--parsed-dir',required=True);s.add_argument('--output-dir',required=True);s.set_defaults(func=summarize)
    args=p.parse_args()
    try:
        if hasattr(args,'targets') and (any(n<0 for n in args.targets) or not sum(args.targets)):
            raise ValueError('Targets must be nonnegative with nonzero total')
        args.func(args)
    except (ValueError,OSError,KeyError,json.JSONDecodeError,zipfile.BadZipFile) as e:
        p.exit(2,f'ERROR: {e}\n')

if __name__=='__main__':
    cli()
