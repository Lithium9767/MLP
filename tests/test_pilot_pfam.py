"""Synthetic parser fixtures only: never written to real pilot raw results."""
import argparse
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
import xml.etree.ElementTree as ET

SPEC=importlib.util.spec_from_file_location('pilot',Path(__file__).resolve().parents[1]/'scripts/pilot_pfam.py')
pilot=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pilot)
RULES=json.loads((Path(__file__).resolve().parents[1]/'configs/pilot_pfam_rules.json').read_text())
SEQ='ACDEFGHIKLMNPQRSTVWY'*4

def source(seq=SEQ,partial=False):
    return {'unique_sequence_id':'TEST_ONLY_1','sequence':seq,'representative_annotation':{'id':'TEST_ONLY_1','organism':'synthetic test','source_file':'test_only.fasta','description':'GvpA partial' if partial else 'GvpA','gvp_types':['GvpA']},'redundant_members_details':[]}

def protein(signature='PF00741',ev='1e-20',end=60):
    return ET.fromstring(f'''<protein><sequence>{SEQ}</sequence><xref id="PILOT_0001"/><matches><hmmer3-match score="60" evalue="{ev}"><signature ac="{signature}" name="gas vesicle"><entry ac="IPR000638"/><signature-library-release library="PFAM" version="TEST"/></signature><locations><hmmer3-location start="1" end="{end}" evalue="1e-20" hmm-start="1" hmm-end="60" hmm-length="70"/></locations></hmmer3-match></matches></protein>''')

class Selection(unittest.TestCase):
    def test_schema_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'input.json';p.write_text('[{"wrong":1}]')
            with self.assertRaisesRegex(ValueError,'unsupported schema'):
                pilot.parse_input(str(p))

    def test_qc_empty_illegal_partial(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'input.json';p.write_text(json.dumps([source(),source(''),source('X'*80),source(partial=True)]))
            rows,qc,info=pilot.parse_input(str(p))
            self.assertEqual(len(rows),1);self.assertEqual(len(qc),3)
            self.assertEqual(rows[0]['sequence_sha256'],pilot.sha(SEQ.encode()))

    def test_seeded_selection_dedup(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'input.json';p.write_text(json.dumps([source(),source(),source(SEQ[:-1]+'A')]))
            rows,_,_=pilot.parse_input(str(p))
            a,qc=pilot.select_rows(rows,[3,0,0],42)
            b,_=pilot.select_rows(rows,[3,0,0],42)
            self.assertEqual([r['sequence_sha256'] for r in a],[r['sequence_sha256'] for r in b])
            self.assertEqual(len(a),1)
            self.assertTrue(any(r['reason']=='exact_duplicate' for r in qc))
            self.assertTrue(any(r['reason'].startswith('near_duplicate') for r in qc))

    def test_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d)/'existing').write_text('keep')
            with self.assertRaises(ValueError):pilot.fresh(d)
            self.assertEqual((Path(d)/'existing').read_text(),'keep')

    def test_reject_wrong_target(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'rules.json';rules=dict(RULES);rules['signature_accession']='PF01132'
            p.write_text(json.dumps(rules))
            with self.assertRaisesRegex(ValueError,'PF00741'):pilot.load_rules(p)

class Parser(unittest.TestCase):
    def hit(self,p):
        return pilot.parse_protein(p,{'internal_id':'PILOT_0001'},SEQ,RULES,'TEST','unit-test-only')

    def test_exact_signature_required(self):
        self.assertEqual(self.hit(protein('PF01132')),[])
        self.assertEqual(self.hit(protein('IPR000638')),[])
        self.assertEqual(self.hit(protein())[0]['hit_status'],'accepted')

    def test_boundary_and_missing_score_review(self):
        self.assertEqual(self.hit(protein(ev='1e-5'))[0]['hit_status'],'needs_review')
        p=protein();p.find('matches/hmmer3-match').attrib.pop('evalue')
        h=self.hit(p)[0];self.assertEqual(h['evalue'],'');self.assertEqual(h['hit_status'],'needs_review')

    def test_missing_hmm_length_review(self):
        p=protein();p.find('matches/hmmer3-match/locations/hmmer3-location').attrib.pop('hmm-length')
        self.assertEqual(self.hit(p)[0]['hit_status'],'needs_review')

    def test_mismatch_and_invalid_coordinates(self):
        with self.assertRaises(ValueError):self.hit(protein(end=999))
        p=protein();p.find('sequence').text='AAAA'
        with self.assertRaises(ValueError):self.hit(p)

    def test_missing_matches_not_negative(self):
        p=protein();p.remove(p.find('matches'))
        with self.assertRaisesRegex(ValueError,'completeness'):self.hit(p)

    def test_mixed_database_rejected(self):
        p=protein();p.find('matches/hmmer3-match/signature/signature-library-release').set('version','OTHER')
        with self.assertRaises(ValueError):self.hit(p)

    def test_preserve_all_locations(self):
        p=protein();loc=p.find('matches/hmmer3-match/locations')
        loc.append(ET.fromstring(ET.tostring(loc[0])))
        self.assertEqual(len(self.hit(p)),2)

    def setup_run(self,d):
        d=Path(d);p=d/'pilot';p.mkdir();r=d/'run';(r/'raw').mkdir(parents=True)
        row=dict(internal_id='PILOT_0001',sequence_length=len(SEQ),sample_group='clear_gvpa',sequence_sha256=pilot.sha(SEQ.encode()),qc_status='pass_pilot_screen',partial_status='not_detected_not_proof_of_completeness')
        pilot.csv_write(p/'pilot_manifest.csv',pilot.MANIFEST,[row]);(p/'pilot_sequences.fasta').write_text('>PILOT_0001\n'+SEQ+'\n')
        pilot.write_json(r/'rules.json',RULES)
        receipt=dict(submitted=True,exit_code=0,scan_status='completed',tool_version='TEST',database_version='TEST',ids=['PILOT_0001'],fasta_sha256=pilot.sha((p/'pilot_sequences.fasta').read_bytes()),manifest_sha256=pilot.sha((p/'pilot_manifest.csv').read_bytes()),rules_sha256=pilot.sha((r/'rules.json').read_bytes()))
        return p,r,receipt

    def parsed_label(self,p,r,d,receipt,xml=None):
        if xml is not None:
            (r/'raw/pilot.xml').write_text(xml)
            receipt['xml_sha256']=pilot.sha((r/'raw/pilot.xml').read_bytes())
        pilot.write_json(r/'run_receipt.json',receipt)
        out=Path(d)/'parsed'
        pilot.parse(argparse.Namespace(pilot_dir=p,run_dir=r,output_dir=out))
        return pilot.csv_read(out/'pilot_labels.csv')[0]

    def test_verified_empty_matches_negative(self):
        with tempfile.TemporaryDirectory() as d:
            p,r,receipt=self.setup_run(d)
            xml=f'<protein-matches interproscan-version="TEST"><protein><sequence>{SEQ}</sequence><xref id="PILOT_0001"/><matches/></protein></protein-matches>'
            self.assertEqual(self.parsed_label(p,r,d,receipt,xml)['label'],'0')

    def test_missing_protein_not_negative(self):
        with tempfile.TemporaryDirectory() as d:
            p,r,receipt=self.setup_run(d)
            row=self.parsed_label(p,r,d,receipt,'<protein-matches interproscan-version="TEST"/>')
            self.assertEqual(row['label'],'');self.assertEqual(row['scan_status'],'result_missing')

    def test_not_submitted_and_failed_blank(self):
        for submitted,code in [(False,None),(True,1)]:
            with self.subTest(submitted=submitted),tempfile.TemporaryDirectory() as d:
                p,r,receipt=self.setup_run(d);receipt.update(submitted=submitted,exit_code=code)
                self.assertEqual(self.parsed_label(p,r,d,receipt)['label'],'')

    def test_truncated_xml_blank(self):
        with tempfile.TemporaryDirectory() as d:
            p,r,receipt=self.setup_run(d)
            row=self.parsed_label(p,r,d,receipt,'<protein-matches')
            self.assertEqual(row['label'],'');self.assertEqual(row['scan_status'],'parse_failed')

    def test_summary_rejects_unsubmitted_negative(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            pilot.csv_write(p/'pilot_labels.csv',pilot.LABELS,[dict(internal_id='PILOT_0001',sample_group='clear_gvpa',scan_status='success',label='0',pf00741_detected='0')])
            pilot.csv_write(p/'pilot_scan_status.csv',pilot.STATUS,[dict(internal_id='PILOT_0001',scan_status='success',submitted=False,result_present=False,parse_status='success')])
            with self.assertRaisesRegex(ValueError,'without verified'):
                pilot.summarize(argparse.Namespace(parsed_dir=p,output_dir=p/'report'))

if __name__=='__main__':unittest.main()
