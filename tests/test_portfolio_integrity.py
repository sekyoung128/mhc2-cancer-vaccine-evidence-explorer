"""Independent checks of the one distributed matrix and historical manifest."""
import csv,hashlib,json,re,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
NAME='mhc2_tcell_evidence_matrix_post_validation_portfolio_v1.csv'
HASH='cfcc5915888465ff00795cf12dddfd330c8fdebfad30ff93c649ee862819f135'
class PortfolioIntegrityTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.m=json.loads((ROOT/'docs/portfolio_release_manifest.json').read_text())
  with (ROOT/'data/processed'/NAME).open(newline='') as f:
   r=csv.DictReader(f);cls.schema=r.fieldnames;cls.rows=list(r)
 def test_hash_schema_inventory(self):
  self.assertEqual(hashlib.sha256((ROOT/'data/processed'/NAME).read_bytes()).hexdigest(),HASH)
  self.assertEqual(self.m['distributed_dataset']['derivative_sha256'],HASH)
  self.assertEqual(self.schema,self.m['distributed_dataset']['column_order'])
  self.assertEqual(len(self.schema),74)
  self.assertEqual([p.name for p in (ROOT/'data/processed').iterdir()],[NAME])
  expected = {
   '.github/workflows/portal-tests.yml',
   '.gitignore',
   'ATTRIBUTION.md',
   'LICENSE',
   'LICENSES/CC-BY-4.0.txt',
   'LICENSES/MIT.txt',
   'README.md',
   'data/LICENSE.md',
   'data/processed/mhc2_tcell_evidence_matrix_post_validation_portfolio_v1.csv',
   'docs/case_study.md',
   'docs/methods.md',
   'docs/portfolio_release_manifest.json',
   'docs/validation.md',
   'portal/app.py',
   'portal/portal_utils.py',
   'portal/requirements.txt',
   'requirements.txt',
   'tests/test_mhc2_portal_utils.py',
   'tests/test_portfolio_integrity.py',
  }
  if (ROOT/'.git').exists():
   import subprocess
   actual=set(subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().rstrip('\0').split('\0'))
  else:
   actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
  self.assertEqual(actual,expected)

 def test_pairs_assays_and_attribution(self):
  self.assertEqual(len(self.rows),7141)
  self.assertEqual(len({(r['reference_pmid'],r['epitope_iri']) for r in self.rows}),7141)
  assays=[]
  for r in self.rows:
   a=json.loads(r['assay_evidence_records']);self.assertEqual(len(a),int(r['tcell_assay_row_count']))
   for x in a:
    self.assertEqual((x['reference_pmid'],x['epitope_iri']),(r['reference_pmid'],r['epitope_iri']))
   assays.extend(a)
  self.assertEqual(len(assays),9145);self.assertEqual(len({a['assay_id'] for a in assays}),9145)
 def test_historical_seven_decisions_and_pending(self):
  decisions=self.m['historical_lineage']['seven_decisions'];self.assertEqual(len(decisions),7)
  expected={(r['reference_pmid'],r['epitope_iri']) for r in decisions}
  actual={(r['reference_pmid'],r['epitope_iri']) for r in self.rows if r['staged_resolution_state']=='resolved_by_human_post_validation_literature_review'}
  self.assertEqual(actual,expected)
  for r in self.rows:
   if (r['reference_pmid'],r['epitope_iri']) in expected:
    self.assertEqual((r['staged_epitope_status'],r['staged_epitope_role'],r['ground_truth_state']),('include','cancer_target','resolved_cancer_target'))
   if r['staged_resolution_state']=='unresolved_pending_review':self.assertEqual((r['staged_epitope_status'],r['staged_epitope_role']),('',''))
 def test_private_inputs_and_manifest(self):
  for name in ['data/raw','data/interim','data/manual','local_literature']:self.assertFalse((ROOT/name).exists())
  text=(ROOT/'docs/portfolio_release_manifest.json').read_text()
  for pattern in [r'/Users/',r'"reviewer"',r'"reviewed_at"',r'MDS-\d{3}',r'\bPr\d{2}\b']:self.assertIsNone(re.search(pattern,text))
  self.assertEqual(self.m['distributed_dataset']['modified_rows'],38)
 def test_approved_note_regressions(self):
  for rows in [self.rows]:
   hpv=[r for r in rows if r['reference_pmid']=='32716518'];self.assertEqual(len(hpv),28)
   for r in hpv:
    self.assertNotIn('Proposed draft only; researcher approval pending.',r['staged_review_note'])
    self.assertIn('machine-assisted review proposal',r['staged_review_note']);self.assertIn('researcher-approved in the VR03 freeze',r['staged_review_note']);self.assertIn('no additional individual sign-off is asserted',r['staged_review_note'])
   targets={'28408606':{'2248624','2248626','2248630'},'33262074':{'1539513','1539704','1539946','1542517'},'39774774':{'2257694'},'15837732':{'1872156','1948755'}}
   for r in rows:
    if r['epitope_iri'].split('/')[-1] in targets.get(r['reference_pmid'],set()):self.assertIsNone(re.search(r'patient \d{4}|MDS-\d{3}|\bPr\d{2}\b|result in P\d{2}',r['staged_review_note']))
  post=self.rows
  for pm,e,labels in [('16000955','1714010',['OT-3.10','OT-W']),('25449106','1319898',['P1.17'])]:
   row=next(r for r in post if r['reference_pmid']==pm and r['epitope_iri'].endswith('/'+e))
   for label in labels:self.assertIn(label,row['staged_review_note'])
