import json
import tempfile
import unittest
from pathlib import Path
from scripts.m3_validate_frozen import approved_artifacts


class FreezeGateTests(unittest.TestCase):
    def test_unapproved_proposal_fails_before_reading_data_or_model(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'approval.json'
            p.write_text(json.dumps({'approved':False}))
            with self.assertRaisesRegex(ValueError,'A approval'):
                approved_artifacts(root,p)

    def test_approval_cannot_unlock_changed_proposal(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'approval.json'
            p.write_text(json.dumps({'approved':True,'approved_by':'test fixture, not real A',
                       'approval_evidence':'synthetic unit test','proposal_sha256':'0'*64,'discovery_receipt_sha256':'0'*64}))
            (root/'candidate_freeze_proposal.json').write_text('{}')
            (root/'run_receipt.json').write_text('{}')
            with self.assertRaisesRegex(ValueError,'exact proposal'):
                approved_artifacts(root,p)
