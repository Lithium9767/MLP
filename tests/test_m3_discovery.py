import csv
import hashlib
import tempfile
import unittest
from pathlib import Path
import numpy as np
from features.esm2_embed import SequenceRecord
from scripts.m3_discovery import checked_coordinates, window_rows, candidate_table


class DiscoveryCoordinateTests(unittest.TestCase):
    def test_unmapped_cluster_does_not_invent_hmm_state_one(self):
        rows=[{'internal_id':str(i),'hmm_states':'','sequence_length':40,'normalized_start':0.,'mapped_fraction':0.} for i in range(2)]
        result=candidate_table(rows,np.array([0,0]),{'seed':42,'random_baseline_repeats':2})
        self.assertEqual(result[0]['dominant_hmm_state'],'')
        self.assertEqual(result[0]['mean_mapped_fraction'],0)
    def record(self):
        seq='ACDEFGHIKLMNPQRSTVWY'*2
        return SequenceRecord('a','seq',seq,hashlib.sha256(seq.encode()).hexdigest(),'discovery','primary')

    def test_context_windows_pool_exact_residues_and_preserve_hmm_insertions(self):
        r=self.record();a=np.arange(40*480,dtype=np.float32).reshape(40,480)
        mapping={'a':{p:(p if p<30 else None) for p in range(1,41)}}
        rows,vectors,baseline=window_rows([r],{'a':a},mapping,30,5)
        self.assertEqual([(x['raw_start'],x['raw_end']) for x in rows],[(1,30),(6,35),(11,40)])
        np.testing.assert_array_equal(vectors[1],a[5:35].mean(axis=0))
        self.assertAlmostEqual(rows[0]['mapped_fraction'],29/30)
        self.assertEqual(rows[0]['hmm_states'].split(';')[-1],'29')
        np.testing.assert_allclose(baseline.sum(axis=1),1,atol=1e-6)

    def test_invalid_embedding_values_are_rejected(self):
        r=self.record();a=np.ones((40,480),dtype=np.float32);a[0,0]=np.nan
        with self.assertRaisesRegex(ValueError,'finite'):
            window_rows([r],{'a':a},{'a':{}},30,5)

    def test_coordinate_mismatch_and_validation_poison_are_rejected(self):
        r=self.record()
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'coordinates.csv'
            row={'internal_id':'a','sequence_sha256':r.sequence_sha256,'split':'validation',
                 'raw_position':'1','residue':'A','hmm_match_state':'1'}
            def save():
                with path.open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=list(row));w.writeheader();w.writerow(row)
            save()
            with self.assertRaisesRegex(ValueError,'split'):checked_coordinates(path,[r])
            row['split']='discovery';row['residue']='Y';save()
            with self.assertRaisesRegex(ValueError,'residue'):checked_coordinates(path,[r])

    def test_duplicate_coordinates_fail_instead_of_silently_shifting(self):
        r=self.record()
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'coordinates.csv';row={'internal_id':'a','sequence_sha256':r.sequence_sha256,
                'split':'discovery','raw_position':'1','residue':'A','hmm_match_state':'1'}
            with p.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(row));w.writeheader();w.writerows([row,row])
            with self.assertRaisesRegex(ValueError,'Ambiguous'):checked_coordinates(p,[r])
