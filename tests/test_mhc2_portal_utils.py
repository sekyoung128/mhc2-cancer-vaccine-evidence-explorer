"""Portal utilities preserve the current cohort and independent evidence dimensions."""
import hashlib
import io
import sys
import unittest
import tempfile
from unittest.mock import patch
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'portal'))
import portal_utils as u

class PortalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.digest = u.file_hash(u.canonical_path())
        cls.frame = u.load_matrix()

    def test_current_invariants(self):
        self.assertEqual(self.frame.shape, (7141, 74))
        self.assertEqual(len(set(u.pair_keys(self.frame))), 7141)
        self.assertEqual(int(u.numeric_counts(self.frame).tcell_assay_row_count.sum()), 9145)
        self.assertEqual(int(self.frame.staged_resolution_state.eq('unresolved_pending_review').sum()), 422)

    def test_role_totals(self):
        self.assertEqual(self.frame.staged_epitope_role.value_counts().to_dict(),
            {'cancer_target': 4517, 'relevant_comparator': 2188, 'non_cancer_antigen': 12, 'control': 2, '': 422})
        self.assertEqual(len(u.filter_pairs(self.frame, {'staged_epitope_role': ['']})), 422)

    def test_response_totals(self):
        self.assertEqual(self.frame.cd4_exact_response_pattern.value_counts().to_dict(),
            {'negative_only': 5403, 'positive_only': 970, 'mixed_positive_negative': 46, 'no_exact_cd4_functional_evidence': 722})
        self.assertEqual(len(u.filter_pairs(self.frame, {'cd4_exact_response_pattern': ['negative_only']})), 5403)

    def test_substring_case_and_literal_matching(self):
        result = u.search(self.frame, ' iqqllmgtv ')
        self.assertGreater(len(result), 0)
        self.assertTrue(result.exact_sequence.str.contains('IQQLLMGTV', regex=False).all())
        self.assertEqual(len(u.search(self.frame, '.*')), 0)
        sample = self.frame.iloc[:2].copy()
        sample['reference_title'] = ['[literal marker]', 'ordinary title']
        self.assertEqual(len(u.search(sample, '[')), 1)

    def test_exact_sequence(self):
        result = u.search(self.frame, 'iqqllmgtvnivvpt', True)
        self.assertGreater(len(result), 0)
        self.assertTrue(result.exact_sequence.eq('IQQLLMGTVNIVVPT').all())
        self.assertEqual(len(u.search(self.frame, 'IQQLLMGTV', True)), 0)

    def test_numeric_exact_ids(self):
        result = u.search(self.frame, '10864677')
        self.assertTrue(result.reference_pmid.eq('10864677').all())
        self.assertEqual(len(u.search(self.frame, '108646')), 0)
        result = u.search(self.frame, '110906')
        self.assertTrue(result.epitope_iri.str.endswith('/110906').all())
        self.assertGreater(len(result), 0)

    def test_full_iri_exact_search(self):
        iri = 'http://www.iedb.org/epitope/110906'
        self.assertTrue(u.search(self.frame, iri).epitope_iri.eq(iri).all())
        self.assertGreater(len(u.search(self.frame, 'https://www.iedb.org/epitope/110906/')), 0)

    def test_blanks_and_unknown_display(self):
        self.assertEqual(u.display('antigen_label', ''), 'Source antigen not recorded')
        self.assertEqual(u.display('staged_epitope_role', ''), 'Unresolved role')
        self.assertEqual(u.display('staged_epitope_status', ''), 'Pending')
        self.assertIn('Not recoverable', u.display('blocking_supported_class2', 'unknown'))
        self.assertIn('No qualifying', u.display('cd4_exact_response_pattern', 'no_exact_cd4_functional_evidence'))

    def test_negative_target_preserved(self):
        row = u.choose_pair(self.frame, ('10864677', 'http://www.iedb.org/epitope/110906'))
        self.assertEqual(row.exact_sequence, 'IQQLLMGTVNIVVPT')
        self.assertEqual(row.staged_epitope_role, 'cancer_target')
        self.assertEqual(row.cd4_exact_response_pattern, 'negative_only')
        self.assertEqual(u.display('cd4_exact_response_pattern', row.cd4_exact_response_pattern), 'Negative-only exact-CD4 observations')

    def test_comparator_preserved(self):
        row = u.choose_pair(self.frame, ('29731959', 'http://www.iedb.org/epitope/1075974'))
        self.assertEqual(row.exact_sequence, 'AATSHPKPTTGHIKP')
        self.assertEqual(row.staged_epitope_role, 'relevant_comparator')

    def test_boundary_wording(self):
        self.assertEqual(u.display('tumour_cell_assay_object_present', 'true'), 'Observed')
        self.assertIn('not recoverable', u.display('native_tumour_recognition_structured_support', 'unavailable'))
        self.assertIn('does not establish biological absence', u.display('natural_presentation_structured_support', 'not_observed_in_structured_data'))
        self.assertIn('not a comprehensive', u.display('staged_individually_human_reviewed', 'no'))

    def test_same_epitope_multiple_publications(self):
        result = u.search(self.frame, 'http://www.iedb.org/epitope/72906')
        self.assertGreater(result.reference_pmid.nunique(), 1)
        self.assertEqual(len(result), len(set(u.pair_keys(result))))

    def test_nonmutation(self):
        before = self.frame.copy(deep=True)
        u.search(self.frame, 'WT1'); u.filter_pairs(self.frame, {'staged_epitope_role': ['cancer_target']})
        u.ordered(self.frame); u.table_view(self.frame.iloc[:10]); u.numeric_counts(self.frame)
        pd.testing.assert_frame_equal(before, self.frame)
        self.assertEqual(self.digest, u.file_hash(u.canonical_path()))

    def test_filters_and_empty_selection(self):
        result = u.filter_pairs(self.frame, {'staged_epitope_role': ['cancer_target', 'relevant_comparator'],
            'cd4_exact_response_pattern': ['negative_only']})
        self.assertEqual(len(result), 5234)  # Current matrix includes two later negative-only target decisions.
        self.assertEqual(len(u.filter_pairs(self.frame, {'staged_epitope_role': []})), 7141)
        with self.assertRaises(ValueError): u.filter_pairs(self.frame, {'unknown': []})

    def test_pagination_and_ordering(self):
        sorted_frame = u.ordered(self.frame)
        page, number, total = u.paginate(sorted_frame, 1)
        self.assertEqual((len(page), number, total), (50, 1, 143))
        last, number, total = u.paginate(sorted_frame, 999)
        self.assertEqual((len(last), number, total), (41, 143, 143))
        with self.assertRaises(ValueError): u.paginate(self.frame, 0)

    def test_json_and_errors(self):
        row = self.frame.iloc[0]
        records = u.parse_assays(row.assay_evidence_records)
        self.assertEqual(len(records), int(row.tcell_assay_row_count))
        for value in ['bad', '{}', '[]', '[{}]']:
            with self.assertRaises(ValueError): u.parse_assays(value)

    def test_validation_rejects_corruption(self):
        with self.assertRaises(ValueError): u.validate_matrix(self.frame.iloc[:-1])
        bad = self.frame.copy(); bad.loc[0, 'staged_epitope_role'] = 'invented'
        with self.assertRaises(ValueError): u.validate_matrix(bad)
        bad = self.frame.copy(); bad.loc[0, 'tcell_assay_row_count'] = ''
        with self.assertRaises(ValueError): u.validate_matrix(bad)
        with self.assertRaises(ValueError): u.parse_boolean('unknown')
        self.assertFalse(u.parse_boolean('false'))

    def test_safe_urls_and_plain_text(self):
        self.assertIsNotNone(u.safe_external_url('https://www.iedb.org/epitope/110906', 'epitope'))
        for url in ['javascript:alert(1)', 'https://www.iedb.org.evil/epitope/110906',
                'https://user@www.iedb.org/epitope/110906', 'https://www.iedb.org:443/epitope/110906']:
            self.assertIsNone(u.safe_external_url(url, 'epitope'))
        self.assertEqual(u.display('reference_title', '&lt;sup&gt;CD4&lt;/sup&gt;'), '<sup>CD4</sup>')

    def test_download_preserves_columns_and_strings(self):
        filtered = u.search(self.frame, '10864677')
        downloaded = pd.read_csv(io.BytesIO(u.download_csv(filtered)), dtype=str, keep_default_na=False)
        pd.testing.assert_frame_equal(downloaded, filtered.reset_index(drop=True))
        self.assertEqual(list(downloaded), list(self.frame))

    def test_empty_defaults_and_resolution_scope(self):
        result = u.filter_pairs(self.frame, {field: [] for field in u.FILTERS})
        self.assertEqual(len(result), 7141)
        self.assertEqual(len(u.resolution_scope(result, 'Resolved')), 6719)
        self.assertEqual(len(u.resolution_scope(result, 'Pending')), 422)
        self.assertEqual(len(u.filter_pairs(u.resolution_scope(result, 'Pending'),
            {'staged_epitope_role': ['cancer_target']})), 0)

    def test_selection_is_explicit_and_retained_by_key(self):
        first, _, _ = u.paginate(self.frame)
        options, selected = u.selection_options(first, self.frame)
        self.assertIsNone(selected)
        key = options[0]
        second, _, _ = u.paginate(self.frame, 2)
        options, selected = u.selection_options(second, self.frame, key)
        self.assertEqual(selected, key)
        self.assertIn(key, options)
        options, selected = u.selection_options(second, second, key)
        self.assertIsNone(selected)

    def test_root_and_canonical_only(self):
        self.assertEqual(u.repository_root(), ROOT)
        with self.assertRaises(ValueError): u.load_matrix(ROOT / 'data/interim/other.csv')

class CurrentDatasetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frame=u.load_matrix()

    def test_one_dataset_and_no_alternative(self):
        self.assertEqual(list(u.DATASETS),['post_validation_v1'])
        self.assertEqual(u.DEFAULT_VERSION,'post_validation_v1')
        with self.assertRaises(ValueError):u.load_matrix(version='frozen_level2')
        self.assertEqual(u.dataset_counts(self.frame),dict(pairs=7141,assays=9145,resolved=6719,pending=422,targets=4517))

    def test_post_decisions_keep_negative_assays(self):
        for pm,ep in [('15837732','1948755'),('22529020','1350619')]:
            row=u.choose_pair(self.frame,(pm,'http://www.iedb.org/epitope/'+ep))
            self.assertEqual(row.staged_epitope_role,'cancer_target')
            self.assertTrue(all(r['qualitative_result']=='Negative' for r in u.parse_assays(row.assay_evidence_records)))
        self.assertEqual(len(u.filter_pairs(self.frame,{'staged_resolution_state':[u.POST_STATE]})),7)

    def test_missing_corrupt_no_fallback(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(u,'repository_root',return_value=Path(tmp)):
            path=u.canonical_path();path.parent.mkdir(parents=True)
            (path.parent/'alternative.csv').write_text('never a fallback')
            with self.assertRaises(FileNotFoundError):u.load_matrix()
            path.write_text('corrupt')
            with self.assertRaisesRegex(ValueError,'SHA-256'):u.load_matrix()

    def test_key_response_and_provenance_corruption_rejected(self):
        for case in ['duplicate','response','provenance','hla']:
            bad=self.frame.copy()
            if case=='duplicate':bad.loc[0,u.KEYS]=bad.loc[1,u.KEYS].to_numpy()
            elif case=='response':bad.at[0,'cd4_exact_response_pattern']='invented'
            elif case=='hla':bad.at[0,'restriction_support_summary']='invented'
            else:
                i=bad.index[bad.staged_resolution_state.eq(u.POST_STATE)][0]
                bad.at[i,'staged_decision_source']='manual_machine_human_validation'
            with self.assertRaises(ValueError):u.validate_matrix(bad)

    def test_app_search_details_export_and_no_selector(self):
        import streamlit as st
        from streamlit.testing.v1 import AppTest
        with patch.object(st,'download_button',wraps=st.download_button) as download:
            at=AppTest.from_file(str(ROOT/'portal/app.py'),default_timeout=30).run()
            self.assertFalse(at.exception)
            self.assertEqual([x.value for x in at.metric],['7,141','9,145','6,719','422','4,517'])
            self.assertFalse(any(x.label=='Dataset version' for x in at.selectbox))
            at.text_input[0].set_value('1948755').run()
            at.selectbox(key='selected_pair').set_value(('15837732','http://www.iedb.org/epitope/1948755')).run()
            self.assertFalse(at.exception)
            text=' '.join(x.value for x in at.text)
            self.assertIn('Scientific review note:',text);self.assertIn('Negative',text)
            self.assertNotIn('Researcher sign-off',text)
            saved=pd.read_csv(io.BytesIO(download.call_args.args[1]),dtype=str,keep_default_na=False)
            pd.testing.assert_frame_equal(saved,u.ordered(u.search(self.frame,'1948755')).reset_index(drop=True))
            self.assertEqual(download.call_args.kwargs['file_name'],u.download_filename())
            captions=' '.join(x.value for x in at.caption)
            self.assertIn('not',captions);self.assertIn('54 confirmed',captions)

    def test_app_pending_and_role_filters(self):
        from streamlit.testing.v1 import AppTest
        at=AppTest.from_file(str(ROOT/'portal/app.py'),default_timeout=30).run()
        at.radio(key='resolution_scope').set_value('Pending').run()
        self.assertFalse(at.exception);self.assertTrue(any('422 matching pairs' in s.value for s in at.subheader))
        at.radio(key='resolution_scope').set_value('All').run()
        at.multiselect(key='filter_staged_epitope_role').set_value(['cancer_target']).run()
        self.assertFalse(at.exception);self.assertTrue(any('4,517 matching pairs' in s.value for s in at.subheader))

    def test_app_missing_asset_fails_closed(self):
        from streamlit.testing.v1 import AppTest
        with tempfile.TemporaryDirectory() as tmp,patch.object(u,'canonical_path',return_value=Path(tmp)/'missing.csv'):
            at=AppTest.from_file(str(ROOT/'portal/app.py'),default_timeout=30).run()
            self.assertFalse(at.exception);self.assertFalse(at.metric)
            self.assertIn('No substitute dataset',at.error[0].value)

    def test_isolated_four_asset_runtime(self):
        import shutil,subprocess
        assets=['portal/app.py','portal/portal_utils.py','portal/requirements.txt',u.SOURCE]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for asset in assets:
                target=root/asset;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/asset,target)
            runner="from streamlit.testing.v1 import AppTest; a=AppTest.from_file('portal/app.py',default_timeout=30).run(); assert not a.exception; assert [x.value for x in a.metric]==['7,141','9,145','6,719','422','4,517']"
            result=subprocess.run([sys.executable,'-B','-c',runner],cwd=root,capture_output=True,text=True,timeout=90)
            self.assertEqual(result.returncode,0,result.stderr)

if __name__=='__main__':unittest.main()
