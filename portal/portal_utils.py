"""Read-only utilities for the current MHC-II matrix; no Streamlit dependency."""
from pathlib import Path
import hashlib
import html
import json
import re
from urllib.parse import urlsplit
import pandas as pd

KEYS = ['reference_pmid', 'epitope_iri']
SOURCE = 'data/processed/mhc2_tcell_evidence_matrix_post_validation_portfolio_v1.csv'
DEFAULT_VERSION = 'post_validation_v1'
DATASETS = {'post_validation_v1': {
    'label': 'Post-validation v1 — Portfolio v1', 'path': SOURCE,
    'sha256': 'cfcc5915888465ff00795cf12dddfd330c8fdebfad30ff93c649ee862819f135',
    'resolved': 6719, 'pending': 422, 'targets': 4517,
    'description': 'Current evidence matrix incorporating seven documented researcher-approved pair-level decisions. Note-only portfolio edits preserve observations and scientific limitations.',
}}
POST_STATE = 'resolved_by_human_post_validation_literature_review'
POST_SOURCE = 'manual_post_validation_literature_review'
POST_KEYS = {(pm, 'http://www.iedb.org/epitope/' + ep) for pm, ep in [
    ('15837732', '1872156'), ('15837732', '1948755'), ('16000955', '1714010'),
    ('22529020', '1350619'), ('25449106', '1319898'), ('29405561', '2257779'), ('31216557', '2188420')]}
POST_COLUMNS = ['staged_epitope_status', 'staged_epitope_role', 'staged_human_resolution',
    'staged_review_note', 'staged_decision_source', 'staged_resolution_state', 'ground_truth_state']
VERSION_BOUNDARY = ('The total evidence universe, raw assay values, and CD4 and HLA evidence classifications '
    'are unchanged. The newer view updates seven pair-level biological-role decisions; remaining unresolved cases remain visible.')
PRIMARY_HISTORY = 'Frozen PRIMARY validation: 54 confirmed / 6 corrected / 15 unresolved out of 75.'
PRIMARY_FOLLOWUP = ('The separate post-validation re-audit resolved seven of the original fifteen unresolved pairs, '
    'leaving eight in the newer view. This follow-up layer does not alter the original validation experiment.')
SEARCH_FIELDS = ['exact_sequence', 'antigen_label', 'reference_pmid', 'epitope_iri', 'reference_title']
FILTERS = {
    'staged_resolution_state': 'Resolution pathway',
    'staged_epitope_role': 'Biological role',
    'cd4_exact_response_pattern': 'Exact-CD4 response',
    'restriction_support_summary': 'Restriction annotation',
    'evidence_attribution_granularity': 'Evidence attribution',
}
CATEGORIES = {
    'staged_resolution_state': {
        'resolved_by_human_validated_rule', 'resolved_by_human_literature_review',
        'resolved_by_human_cluster_rule_review', 'resolved_by_human_machine_validation',
        'resolved_by_individual_rule_spotcheck', 'resolved_by_human_validated_machine_subgroup',
        'resolved_by_human_individual_review', 'unresolved_pending_review'},
    'staged_epitope_role': {'cancer_target', 'relevant_comparator', 'non_cancer_antigen', 'control', ''},
    'staged_epitope_status': {'include', 'exclude', ''},
    'cd4_exact_response_pattern': {'positive_only', 'negative_only', 'mixed_positive_negative',
        'other_or_unknown_only', 'mixed_with_unknown', 'no_exact_cd4_functional_evidence'},
    'restriction_support_summary': {'functional_specific_hla_code_support', 'class2_subset_code_support',
        'binding_assay_code_support', 'prediction_only', 'prediction_with_other_support',
        'cd4_function_without_specific_hla', 'generic_class2_context',
        'cited_reference_or_linked_context', 'insufficient_or_unclassified'},
    'evidence_attribution_granularity': {'exact_pair', 'subcluster_linked', 'reference_linked', 'unresolved'},
    'assay_antigen_relation': {'exact_only', 'exact_and_nonexact', 'nonexact_only', 'unknown', 'no_assay_rows'},
    'immunogen_relation': {'exact_epitope', 'overlapping_or_related_peptide', 'peptide_pool',
        'protein_or_construct', 'tumour_or_cell_material', 'other_antigen', 'mixed', 'unknown', 'not_available'},
    'staged_decision_source': {'human_validated_relationship_rule', 'manual_literature_review',
        'manual_cluster_rule_review', 'manual_machine_human_validation',
        'manual_relationship_rule_spotcheck', 'human_validated_machine_subgroup_rule', 'manual_individual_review', ''},
    'staged_relationship_rule': {'RULE_MUTANT_TARGET', 'RULE_WT_COMPARATOR', ''},
    'staged_individually_human_reviewed': {'yes', 'no'},
    'staged_rule_human_validated': {'yes', 'no'},
    'ground_truth_state': {'resolved_cancer_target', 'resolved_relevant_comparator',
        'resolved_non_cancer_antigen', 'resolved_control', 'resolved_other', 'unresolved_pending_review'},
    'native_tumour_recognition_structured_support': {'unavailable'},
    'natural_presentation_structured_support': {'not_observed_in_structured_data'},
}
COUNT_FIELDS = ['tcell_assay_row_count', 'exact_assay_antigen_row_count', 'nonexact_assay_antigen_row_count',
    'unknown_assay_antigen_relation_row_count', 'cd4_functional_exact_row_count',
    'cd4_functional_nonexact_row_count', 'cd4_exact_positive_row_count', 'cd4_exact_negative_row_count',
    'cd4_exact_other_or_unknown_row_count', 'any_functional_row_count']
BOOL_FIELDS = ['cd4_functional_exact_present', 'cd4_functional_nonexact_present',
    'any_functional_assay_present', 'any_positive_functional_present', 'any_negative_functional_present',
    'any_mixed_functional_present', 'functional_specific_hla_code_support_present', 'mhc_subset_code_present',
    'mhc_binding_assay_code_present', 'generic_class2', 'cd4_functional_without_exact_restriction',
    'prediction_restriction_present', 'prediction_supported_restriction', 'prediction_only_restriction',
    'linked_or_reference_class2', 'insufficient_restriction', 'epitope_immunogen_mismatch_present',
    'protein_or_construct_assay_present', 'cell_based_assay_present',
    'linked_antigen_processing_evidence_present', 'linked_processing_candidate', 'tumour_cell_assay_object_present',
    'has_exact_pair_evidence', 'has_subcluster_linked_evidence', 'has_reference_linked_evidence', 'has_unresolved_attribution']
UNKNOWN_FIELDS = ['experimental_allele_specific_functional', 'blocking_supported_class2',
    'allele_specific_nonfunctional', 'tumour_lysate_assay_present', 'engineered_cell_assay_object_present',
    'peptide_pulsed_target_assay_present', 'native_tumour_recognition_structurally_confirmed']
TEXT_FIELDS = ['reference_pmid', 'epitope_iri', 'epitope_name', 'reference_title', 'reference_iedb_iri',
    'source_molecule', 'exact_sequence', 'antigen_label', 'staged_human_resolution', 'staged_review_note',
    'mhc_restriction_names', 'mhc_restriction_codes', 'assay_ids', 'assay_methods',
    'assay_response_types', 'assay_evidence_records']
REQUIRED = set(CATEGORIES) | set(COUNT_FIELDS) | set(BOOL_FIELDS) | set(UNKNOWN_FIELDS) | set(TEXT_FIELDS)
RESPONSE_LABELS = {
    'positive_only': 'Positive-only exact-CD4 observations',
    'negative_only': 'Negative-only exact-CD4 observations',
    'mixed_positive_negative': 'Mixed positive and negative exact-CD4 observations',
    'other_or_unknown_only': 'Other or unknown exact-CD4 observations',
    'mixed_with_unknown': 'Mixed exact-CD4 observations including unknown results',
    'no_exact_cd4_functional_evidence': 'No qualifying exact-CD4 evidence recorded',
}
RESTRICTION_LABELS = dict(zip([
    'functional_specific_hla_code_support', 'class2_subset_code_support', 'binding_assay_code_support',
    'prediction_only', 'prediction_with_other_support', 'cd4_function_without_specific_hla',
    'generic_class2_context', 'cited_reference_or_linked_context', 'insufficient_or_unclassified'], [
    'Specific-HLA functional-code support', 'Class-II subset-code support', 'Binding-code support in T-cell metadata',
    'Prediction-only restriction annotation', 'Prediction with other structured support',
    'CD4 function without specific-HLA structured support', 'Generic class-II context',
    'Cited/reference-linked class-II context', 'Insufficient or unclassified structured restriction support']))
ATTRIBUTION_LABELS = {'exact_pair': 'Exact-pair observation', 'subcluster_linked': 'Subcluster-linked observation',
    'reference_linked': 'Reference-linked observation', 'unresolved': 'Unresolved assay attribution'}
BOUNDARY_NATIVE = 'Native tumour recognition is not recoverable from the current structured metadata.'
BOUNDARY_NATURAL = 'Natural presentation was not observed in the structured data; this does not establish biological absence.'
RESTRICTION_WARNING = 'Restriction categories are structured annotation proxies and do not by themselves establish experimentally confirmed allele restriction.'
TABLE_COLUMNS = dict(zip(['exact_sequence', 'antigen_label', 'reference_pmid', 'staged_epitope_status',
    'staged_epitope_role', 'cd4_exact_response_pattern', 'restriction_support_summary'],
    ['Peptide', 'Antigen', 'PMID', 'Inclusion', 'Role', 'CD4 response', 'Restriction']))
ASSAY_FIELDS = ['assay_id', 'functional_group', 'effector_cell', 'antigen_name', 'antigen_relation',
    'antigen_object_type', 'qualitative_result', 'mhc_name', 'mhc_evidence_code', 'reference_location']


def repository_root():
    return Path(__file__).resolve().parents[1]


def dataset(version=DEFAULT_VERSION):
    require(version in DATASETS, 'Unknown dataset version')
    return DATASETS[version]


def canonical_path(version=DEFAULT_VERSION):
    return repository_root() / dataset(version)['path']


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def parse_boolean(value):
    require(value in ('true', 'false'), 'Invalid literal observation flag: ' + str(value))
    return value == 'true'


def parse_assays(value):
    try:
        records = json.loads(value)
    except (json.JSONDecodeError, TypeError) as exc:
        raise ValueError('Malformed assay evidence JSON') from exc
    require(isinstance(records, list) and bool(records), 'Assay evidence must be a nonempty JSON array')
    for record in records:
        require(isinstance(record, dict) and set(ASSAY_FIELDS + KEYS + ['exact_antigen_state']) <= set(record),
                'Assay record schema is incomplete')
        require(all(isinstance(v, str) for v in record.values()), 'Assay values must remain source strings')
    return records


def validate_matrix(frame, version=DEFAULT_VERSION):
    metadata = dataset(version)
    require(frame.shape == (7141, 74) and set(frame) == REQUIRED, 'Canonical matrix schema/shape drift')
    require(not frame.isna().any().any(), 'Missing values must be preserved as source strings')
    require(frame[KEYS].ne('').all().all() and not frame.duplicated(KEYS).any(), 'Blank or duplicate pair key')
    for col, allowed in CATEGORIES.items():
        if version == 'post_validation_v1':
            if col == 'staged_resolution_state': allowed = allowed | {POST_STATE}
            if col == 'staged_decision_source': allowed = allowed | {POST_SOURCE}
        require(set(frame[col]) <= allowed, 'Unknown category in ' + col)
    for col in COUNT_FIELDS:
        require(frame[col].str.fullmatch(r'\d+').all(), 'Invalid nonnegative count in ' + col)
    for col in BOOL_FIELDS:
        require(frame[col].isin(['true', 'false']).all(), 'Invalid observation flag in ' + col)
    for col in UNKNOWN_FIELDS:
        require(frame[col].eq('unknown').all(), 'Unexpected unsupported endpoint in ' + col)
    n = frame[COUNT_FIELDS].astype('int64')
    require(int(n.tcell_assay_row_count.sum()) == 9145, 'Assay total drift')
    require((n.exact_assay_antigen_row_count + n.nonexact_assay_antigen_row_count +
             n.unknown_assay_antigen_relation_row_count).eq(n.tcell_assay_row_count).all(), 'Antigen partition drift')
    require((n.cd4_exact_positive_row_count + n.cd4_exact_negative_row_count +
             n.cd4_exact_other_or_unknown_row_count).eq(n.cd4_functional_exact_row_count).all(), 'CD4 partition drift')
    pending = frame.staged_resolution_state.eq('unresolved_pending_review')
    require(int(pending.sum()) == metadata['pending'] and int((~pending).sum()) == metadata['resolved'], 'Resolution totals drift')
    require(frame.loc[pending, ['staged_epitope_status', 'staged_epitope_role']].eq('').all().all(), 'Pending decisions changed')
    require(frame.loc[~pending, 'staged_epitope_status'].isin(['include', 'exclude']).all(), 'Resolved inclusion missing')
    require(frame.staged_epitope_role.value_counts().to_dict() == {
        'cancer_target': metadata['targets'], 'relevant_comparator': 2188, 'non_cancer_antigen': 12,
        'control': 2, '': metadata['pending']}, 'Role totals drift')
    if version == 'post_validation_v1':
        resolved_new = frame.staged_resolution_state.eq(POST_STATE)
        require(set(pair_keys(frame.loc[resolved_new])) == POST_KEYS and
                frame.staged_decision_source.eq(POST_SOURCE).equals(resolved_new), 'Post-validation provenance/key drift')
        require(frame.loc[resolved_new, 'staged_epitope_status'].eq('include').all() and
                frame.loc[resolved_new, 'staged_epitope_role'].eq('cancer_target').all() and
                frame.loc[resolved_new, 'ground_truth_state'].eq('resolved_cancer_target').all(), 'Post-validation labels drift')
    require(frame.cd4_exact_response_pattern.value_counts().to_dict() == {
        'negative_only': 5403, 'positive_only': 970, 'mixed_positive_negative': 46,
        'no_exact_cd4_functional_evidence': 722}, 'Response totals drift')
    require(frame.exact_sequence.equals(frame.epitope_name) and frame.antigen_label.equals(frame.source_molecule), 'Identity alias drift')
    assay_ids = set()
    for row in frame.to_dict('records'):
        records = parse_assays(row['assay_evidence_records'])
        require(len(records) == int(row['tcell_assay_row_count']), 'Assay JSON count differs')
        for record in records:
            require(all(record[k] == row[k] for k in KEYS), 'Assay assigned to wrong pair')
            require(record['assay_id'] not in assay_ids, 'Repeated assay identifier')
            assay_ids.add(record['assay_id'])
    return frame


def load_matrix(path=None, version=DEFAULT_VERSION):
    metadata = dataset(version)
    path = canonical_path(version) if path is None else Path(path)
    require(path.resolve() == canonical_path(version).resolve(), 'Only the selected canonical matrix may be loaded')
    require(file_hash(path) == metadata['sha256'], 'Protected dataset SHA-256 differs')
    frame = validate_matrix(pd.read_csv(path, dtype=str, keep_default_na=False), version)
    return frame


def dataset_counts(frame):
    pending = frame.staged_resolution_state.eq('unresolved_pending_review')
    return {'pairs': len(frame), 'assays': int(numeric_counts(frame).tcell_assay_row_count.sum()),
            'resolved': int((~pending).sum()), 'pending': int(pending.sum()),
            'targets': int(frame.staged_epitope_role.eq('cancer_target').sum())}


def download_filename(version=DEFAULT_VERSION):
    dataset(version)
    return 'mhc2_' + version + '_portfolio_v1_filtered_evidence_pairs.csv'


def numeric_counts(frame):
    """Separate numeric view; retain original canonical strings for downloads."""
    return frame[COUNT_FIELDS].astype('int64').copy()


def search(frame, query='', exact_sequence=False):
    query = query.strip()
    if not query:
        return frame.copy()
    if exact_sequence:
        mask = frame.exact_sequence.str.casefold().eq(query.casefold())
    elif query.isascii() and query.isdigit():
        mask = frame.reference_pmid.eq(query) | frame.epitope_iri.str.rsplit('/', n=1).str[-1].eq(query)
    elif safe_external_url(query, 'epitope'):
        mask = frame.epitope_iri.str.rsplit('/', n=1).str[-1].eq(query.rstrip('/').rsplit('/', 1)[-1])
    else:
        mask = pd.Series(False, index=frame.index)
        for col in SEARCH_FIELDS:
            mask |= frame[col].str.casefold().str.contains(query.casefold(), regex=False)
    return frame.loc[mask].copy()


def filter_pairs(frame, selections):
    result = frame
    for col, values in selections.items():
        require(col in FILTERS, 'Unsupported filter: ' + col)
        allowed = CATEGORIES[col] | ({POST_STATE} if col == 'staged_resolution_state' else set())
        require(set(values) <= allowed, 'Invalid filter values: ' + col)
        if values:
            result = result.loc[result[col].isin(values)]
    return result.copy()


def resolution_scope(frame, scope='All'):
    require(scope in ['All', 'Resolved', 'Pending'], 'Invalid resolution scope')
    pending = frame.staged_resolution_state.eq('unresolved_pending_review')
    if scope == 'All':
        return frame.copy()
    return frame.loc[pending if scope == 'Pending' else ~pending].copy()


def selection_options(page, matching, prior=None):
    """No implicit selection; preserve a full key only while still matched."""
    options = pair_keys(page)
    selected = prior if prior in set(pair_keys(matching)) else None
    if selected is not None and selected not in options:
        options = [selected] + options
    return options, selected


def ordered(frame):
    result = frame.assign(_pmid=pd.to_numeric(frame.reference_pmid),
                          _epitope=pd.to_numeric(frame.epitope_iri.str.rsplit('/', n=1).str[-1]))
    return result.sort_values(['_pmid', '_epitope'], kind='stable').drop(columns=['_pmid', '_epitope']).copy()


def paginate(frame, page=1, page_size=50):
    require(page_size > 0 and page >= 1, 'Invalid pagination')
    total = max(1, (len(frame) + page_size - 1) // page_size)
    page = min(page, total)
    return frame.iloc[(page - 1) * page_size:page * page_size].copy(), page, total


def pair_keys(frame):
    return list(frame[KEYS].itertuples(index=False, name=None))


def choose_pair(frame, key):
    match = frame.loc[frame.reference_pmid.eq(key[0]) & frame.epitope_iri.eq(key[1])]
    require(len(match) == 1, 'Selected pair is not uniquely represented')
    return match.iloc[0]


def display(col, value):
    if col in ['antigen_label', 'source_molecule'] and value == '':
        return 'Source antigen not recorded'
    if col == 'staged_epitope_status':
        return {'': 'Pending', 'include': 'Included in evidence cohort', 'exclude': 'Excluded from evidence cohort'}[value]
    if col == 'staged_epitope_role' and value == '':
        return 'Unresolved role'
    if col == 'cd4_exact_response_pattern':
        return RESPONSE_LABELS[value]
    if col == 'restriction_support_summary':
        return RESTRICTION_LABELS[value]
    if col == 'evidence_attribution_granularity':
        return ATTRIBUTION_LABELS[value]
    if col == 'native_tumour_recognition_structured_support':
        return BOUNDARY_NATIVE
    if col == 'natural_presentation_structured_support':
        return BOUNDARY_NATURAL
    if col in UNKNOWN_FIELDS:
        return 'Not recoverable from the current structured metadata'
    if col in BOOL_FIELDS:
        return 'Observed' if parse_boolean(value) else 'Not observed in the structured data'
    if col == 'staged_individually_human_reviewed':
        return 'Historical spotcheck flag: ' + value + ' (not a comprehensive human-review indicator)'
    if col == 'staged_rule_human_validated':
        return 'Historical rule-validation flag: ' + value
    if not value:
        return 'Not recorded'
    return value.replace('_', ' ').capitalize() if col in CATEGORIES else html.unescape(value)


def table_view(frame):
    return pd.DataFrame({label: frame[col].map(lambda value: display(col, value))
                         for col, label in TABLE_COLUMNS.items()}, index=frame.index)


def safe_external_url(value, kind):
    try:
        u = urlsplit(value)
        if u.scheme not in ['http', 'https'] or u.username or u.password or u.port:
            return None
        if kind == 'pubmed' and u.hostname == 'pubmed.ncbi.nlm.nih.gov' and re.fullmatch(r'/\d+/?', u.path):
            return value if not u.query and not u.fragment else None
        if kind in ['epitope', 'reference', 'assay'] and u.hostname in ['www.iedb.org', 'iedb.org']:
            return value if re.fullmatch('/' + kind + r'/\d+/?', u.path) and not u.query and not u.fragment else None
    except ValueError:
        pass
    return None


def download_csv(frame):
    return frame.to_csv(index=False, lineterminator='\n').encode('utf-8')
