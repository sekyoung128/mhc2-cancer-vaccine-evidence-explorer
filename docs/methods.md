# Methods and evidence boundaries

The evidence unit is publication PMID × IEDB epitope IRI. Identical peptides in different publications remain distinct. Each pair aggregates linked assay records while retaining same-assay associations in JSON; names, HLA codes and measurements must not be cross-combined into invented experiments.

## Research workflow

1. **IEDB exports:** acquire experimental records and retain publication, epitope and assay identifiers, source metadata and original measurements.
2. **Reference and epitope curation:** assess cancer context, then distinguish target, relevant comparator and control roles within each publication. Machine-assisted screening and literature assessment produce proposals, not approval authority.
3. **Documented decisions:** apply recorded researcher-approved decisions and human-validated deterministic rules through exact pair keys. Rule-based resolution is distinguished from individual review; unresolved cases retain explicit pending state.
4. **Canonical evidence matrix:** aggregate linked assays by publication–epitope pair, preserving same-assay provenance and computing counts, response patterns, restriction-support categories and attribution flags.
5. **Streamlit exploration:** load the current approved matrix with hash and schema checks, then expose search, filters, detail and exports without changing scientific values.

The researcher directed the curation purpose and authorised documented decision freezes; coding and literature assessment were AI-assisted. This record does not establish independently completed human verification of every assay or proposal.

Source-derived identifiers, measurements and metadata remain distinct from machine-assisted proposals, researcher-approved annotations and computed evidence fields. Computation organises observations; it does not add experimental validation. This checkout distributes the current matrix and portal utilities, **not the full matrix-building scripts or original exports**.

## Independent dimensions

- **Inclusion and role:** cohort relevance and intended publication-specific purpose (target, relevant comparator, control or other). Neither follows automatically from assay positivity. Pending roles remain blank.
- **CD4 response:** exact-CD4 functional observations require qualifying functional modality, explicitly CD4 effectors and exact assay-antigen attribution. Negative-only and mixed results remain distinct. No qualifying observation is not biological absence. PBMC or DTH findings must not become purified-CD4 proof.
- **Restriction:** specific-HLA functional codes, subset codes, binding codes, prediction and linked context are reported separately. These are metadata proxies, not uniform experimental allele validation.
- **Attribution:** exact pair, subcluster, reference-linked or unresolved evidence. Protein/construct or tumour-cell objects do not by themselves prove native tumour recognition or naturally presented peptide length.
- **Provenance:** resolution pathway, manual source and preserved scientific review notes remain visible. Historical individual-review flags have limited scope.

## Why binding remains a separate arm

The original research inventory found limited provenance-supported linkage between binding and T-cell evidence: 17 T-cell pairs shared both PMID and epitope IRI with binding records, 51 shared an epitope across publications, and eight shared publication, specific HLA and exact-antigen CD4 functional evidence. These are distinct overlap definitions, not sequential stages of a demonstrated biological chain.

Binding was therefore retained separately rather than assuming a unified binding → CD4 function → tumour recognition pathway. Cross-publication molecular overlap is not the same experiment; unavailable binding evidence is not a negative result. These inventory findings are documented in the original research evidence-hierarchy audit, not reconstructed by this public checkout. The portal does not integrate the binding arm.

## Dataset lineage

The distributed dataset is Post-validation v1 only. Its historical lineage includes seven exact pair-level role decisions and their provenance/derived ground-truth fields, with assay observations, CD4 and restriction evidence preserved. Frozen Level-2 remains a historical validation reference and is not distributed.

The portfolio derivative modifies only `staged_review_note` in 38 rows to clarify historical proposal wording and minimise unnecessary subject linkage, while retaining experimental clone identifiers and scientific limitations. The release manifest records the current source and derivative hashes, note-edit keys and historical seven-decision lineage. Original canonical approval authority and original datasets remain privately retained.

The loader fails closed on a missing or corrupted current asset, validates its hash, schema, counts and seven documented decision keys, and never loads another matrix. Historical comparisons and full source reconstruction require retained research inputs outside this checkout.
