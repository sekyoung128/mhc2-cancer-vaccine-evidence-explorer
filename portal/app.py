"""Read-only Streamlit view over the current evidence matrix."""
import sys
from pathlib import Path
import pandas as pd
import streamlit as st
sys.path.insert(0, str(Path(__file__).resolve().parent))
import portal_utils as u

st.set_page_config(page_title='MHC-II Cancer Evidence Explorer', layout='wide')

@st.cache_data(show_spinner='Validating the evidence matrix…')
def cached_matrix(path, digest):
    frame = u.load_matrix(path)
    u.require(u.file_hash(path) == digest, 'Source changed while loading')
    return frame


def fields(row, names):
    for name in names:
        st.text(name.replace('_', ' ').capitalize() + ': ' + u.display(name, row[name]))


def detail(row):
    st.divider()
    st.subheader('Selected publication–epitope pair')
    for column, label, field in zip(st.columns(6),
            ['Relevance', 'Role', 'CD4 response', 'Restriction', 'Attribution', 'Provenance'],
            ['staged_epitope_status', 'staged_epitope_role', 'cd4_exact_response_pattern',
             'restriction_support_summary', 'evidence_attribution_granularity', 'staged_resolution_state']):
        with column:
            with st.container(border=True):
                st.caption(label)
                st.text(u.display(field, row[field]))
    st.subheader('Identity')
    fields(row, ['exact_sequence', 'antigen_label', 'reference_title', 'reference_pmid', 'epitope_iri', 'reference_iedb_iri'])
    for column, (label, value, kind) in zip(st.columns(3), [
            ('PubMed', 'https://pubmed.ncbi.nlm.nih.gov/' + row.reference_pmid + '/', 'pubmed'),
            ('IEDB epitope', row.epitope_iri, 'epitope'), ('IEDB reference', row.reference_iedb_iri, 'reference')]):
        url = u.safe_external_url(value, kind)
        if url:
            column.link_button(label, url)
    st.caption('Role describes the publication-specific intended biological or experimental purpose. It does not indicate assay positivity or vaccine efficacy. Negative response does not assign comparator or control biology.')
    st.text(u.display('cd4_exact_response_pattern', row.cd4_exact_response_pattern))
    st.caption('Exact-CD4 observations require reviewed functional modality, explicitly named CD4 cells, and exact assay-antigen identity. No qualifying evidence is not a negative response.')
    st.caption(u.RESTRICTION_WARNING)
    with st.expander('Assay counts / broader functional observations'):
        fields(row, ['tcell_assay_row_count', 'cd4_functional_exact_row_count', 'cd4_exact_positive_row_count',
            'cd4_exact_negative_row_count', 'cd4_exact_other_or_unknown_row_count',
            'cd4_functional_nonexact_present', 'cd4_functional_nonexact_row_count',
            'any_functional_assay_present', 'any_functional_row_count', 'any_positive_functional_present',
            'any_negative_functional_present', 'any_mixed_functional_present'])
    with st.expander('Restriction code details'):
        fields(row, ['restriction_support_summary', 'functional_specific_hla_code_support_present',
            'mhc_subset_code_present', 'mhc_binding_assay_code_present', 'generic_class2',
            'cd4_functional_without_exact_restriction', 'prediction_restriction_present',
            'prediction_only_restriction', 'linked_or_reference_class2', 'insufficient_restriction',
            'mhc_restriction_names', 'mhc_restriction_codes'] + u.UNKNOWN_FIELDS[:3])
    with st.expander('Attribution counts and scope'):
        fields(row, ['evidence_attribution_granularity', 'has_exact_pair_evidence', 'has_subcluster_linked_evidence',
            'has_reference_linked_evidence', 'has_unresolved_attribution', 'assay_antigen_relation',
            'exact_assay_antigen_row_count', 'nonexact_assay_antigen_row_count', 'unknown_assay_antigen_relation_row_count'])
        st.caption('Linked observations may be valid but are not automatically exact-pair evidence. Assay attribution and pending ground truth are distinct.')
    with st.expander('Review provenance details'):
        fields(row, ['staged_resolution_state', 'staged_decision_source', 'staged_relationship_rule',
            'staged_human_resolution', 'staged_review_note', 'staged_individually_human_reviewed', 'staged_rule_human_validated'])
        st.caption('This portfolio dataset contains note-only privacy and historical-provenance edits. These edits are not new biological validation; original approval records remain private.')
        st.caption('Historical spotcheck flags do not enumerate all human review. A “no” flag must not be interpreted as never human reviewed.')
    with st.expander('Downstream evidence boundaries'):
        fields(row, ['immunogen_relation', 'epitope_immunogen_mismatch_present', 'protein_or_construct_assay_present',
            'cell_based_assay_present'])
        st.text('Protein/cell processing-candidate context: ' + u.display('linked_antigen_processing_evidence_present', row.linked_antigen_processing_evidence_present))
        st.text('Tumour-cell assay object: ' + u.display('tumour_cell_assay_object_present', row.tumour_cell_assay_object_present))
        fields(row, ['native_tumour_recognition_structured_support', 'natural_presentation_structured_support'])
        st.caption('Immunogen mismatch is not automatically an error. Object context alone does not confirm endogenous processing or native tumour recognition.')
    with st.expander('Show assay-level evidence'):
        try:
            records = u.parse_assays(row.assay_evidence_records)
            st.dataframe(pd.DataFrame(records)[u.ASSAY_FIELDS], hide_index=True)
        except ValueError as exc:
            st.error('Data error: ' + str(exc))
        st.caption('Each row preserves same-assay associations. Aggregated names/codes must not be cross-combined into inferred experiments.')


def main():
    st.title('MHC-II Cancer Evidence Explorer')
    st.caption('A view of Class II T-cell evidence for cancer-vaccine research')
    st.write('Search publication-specific MHC-II T-cell evidence while keeping cohort relevance, biological role, functional response, restriction annotation, evidence attribution and provenance separate.')
    with st.expander('What this portal shows — and its interpretation boundaries', expanded=True):
        st.write('What this portal shows')
        st.write('Publication-specific MHC-II T-cell evidence is displayed across separate dimensions: cohort relevance, biological role, exact-antigen CD4 response, restriction annotation, evidence attribution and provenance. Linked evidence retains its scope; unresolved decisions remain visible.')
        st.write('What it does not show')
        st.write('Target role does not imply positive CD4 immunogenicity, experimental allele-specific restriction, natural presentation, native tumour recognition or vaccine efficacy. Negative response does not assign comparator or control biology. Missing evidence is not negative evidence. Restriction categories are structured proxies; native tumour recognition and natural presentation are not established atlas-wide.')
    metadata = u.dataset()
    st.subheader(metadata['label'])
    st.caption(metadata['description'])
    st.caption(u.VERSION_BOUNDARY)
    st.caption(u.PRIMARY_HISTORY)
    st.caption(u.PRIMARY_FOLLOWUP)
    try:
        path = u.canonical_path()
        digest = u.file_hash(path)
        frame = cached_matrix(str(path), digest)
    except (OSError, ValueError, pd.errors.ParserError):
        st.error('Post-validation v1 is unavailable or failed validation. No substitute dataset is shown.')
        st.caption('Required dataset: ' + metadata['path'])
        st.stop()
    counts = u.dataset_counts(frame)
    for col, label, value in zip(st.columns(5), ['Evidence pairs', 'Assay rows', 'Resolved pairs', 'Pending pairs', 'Cancer-target pairs'],
            [counts[k] for k in ['pairs', 'assays', 'resolved', 'pending', 'targets']]):
        col.metric(label, f'{value:,}')
    selections = {}
    with st.sidebar:
        st.header('Filter evidence pairs')
        scope = st.radio('Resolution scope', ['All', 'Resolved', 'Pending'], horizontal=True, key='resolution_scope')
        if st.button('Reset filters'):
            for field in u.FILTERS:
                st.session_state['filter_' + field] = []
        placeholders = ['All resolution pathways', 'All roles', 'All response patterns',
                        'All restriction categories', 'All attribution categories']
        for (field, label), placeholder in zip(u.FILTERS.items(), placeholders):
            options = sorted(frame[field].unique())
            selections[field] = st.multiselect(label, options, key='filter_' + field,
                placeholder=placeholder, help='Leave empty to include all values.',
                format_func=lambda value, field=field: u.display(field, value))
        st.caption('Empty means all values. OR within each filter; AND across filters and resolution scope.')
    query = st.text_input('Search peptide sequence, source antigen, PMID, IEDB epitope ID or publication title')
    st.caption('Try a peptide sequence, PMID, antigen name or publication title.')
    exact = st.checkbox('Exact peptide sequence')
    matching = u.ordered(u.filter_pairs(u.search(u.resolution_scope(frame, scope), query, exact), selections))
    st.subheader(f'{len(matching):,} matching pairs')
    st.caption('Pairs are publication PMID × IEDB epitope. The same peptide in different publications remains separate.')
    st.download_button('Download filtered pairs', u.download_csv(matching),
        file_name=u.download_filename(), mime='text/csv')
    st.caption('Filtered export of ' + metadata['label'] + '; includes all matching pairs and original columns/values, not only the visible page. Review limitations and raw assay evidence are retained.')
    if matching.empty:
        st.info('No pairs match the current search and filters.')
        st.session_state.pop('selected_pair', None)
    else:
        pages = max(1, (len(matching) + 49) // 50)
        st.session_state['page'] = min(int(st.session_state.get('page', 1)), pages)
        page_number = st.number_input('Page', min_value=1, max_value=pages, step=1, key='page')
        page, current, total = u.paginate(matching, int(page_number))
        st.caption(f'Page {current} of {total} · up to 50 pairs per page')
        st.dataframe(u.table_view(page), hide_index=True, height=400,
            column_config={
                'Restriction': st.column_config.TextColumn('Restriction', width=380),
                'PMID': st.column_config.TextColumn('PMID', width=90),
                'Inclusion': st.column_config.TextColumn('Inclusion', width=130),
            })
        options, prior = u.selection_options(page, matching, st.session_state.get('selected_pair'))
        st.session_state['selected_pair'] = prior
        def selection_label(key):
            row = u.choose_pair(matching, key)
            return f'PMID {key[0]} · IEDB {key[1].rsplit("/", 1)[-1]} · {row.exact_sequence}'
        selected = st.selectbox('Select an evidence pair', options, index=None,
            key='selected_pair', format_func=selection_label, placeholder='Choose a pair to inspect')
        if selected is not None:
            detail(u.choose_pair(matching, selected))
        else:
            st.info('Search or select an evidence pair to inspect its evidence dimensions.')
    st.divider()
    st.caption('Read-only view. Binding evidence is not integrated. No vaccine design, clinical decision or predictive recommendations.')
    st.caption('Active matrix SHA-256: ' + digest)


if __name__ == '__main__':
    main()
