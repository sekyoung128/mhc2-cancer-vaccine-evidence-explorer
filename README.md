# MHC-II Cancer Vaccine Evidence Atlas

A research portfolio for exploring publication-specific MHC-II T-cell evidence in cancer-vaccine research, with assay observations, biological roles and interpretation limits kept separate.

## Overview

A peptide described in a cancer study is not necessarily an intended tumour target. It may be a native comparison reagent, a control, or a related peptide used in a different experiment. Binding annotations and CD4 responses also answer different questions. Combining them without tracing the experimental context can overstate the evidence for a vaccine candidate.

This project developed from an IEDB evidence inventory into a curated publication–epitope atlas. Machine-assisted screening and literature assessment supported the work; documented researcher decisions and validated rules supplied decision authority. The resulting explorer makes the reasoning and its limitations available alongside the original assay observations.

## What the atlas contains

The atlas contains **7,141 publication–epitope pairs**, linked to **9,145 T-cell assay records** from the Immune Epitope Database (IEDB). There are 6,719 resolved pairs, 422 pending pairs and 4,517 cancer-target pairs.

Each pair retains peptide and antigen identity, publication and IEDB links, biological role, exact-antigen CD4 response, restriction annotations, assay attribution and review notes. Target role does not imply a positive response or vaccine efficacy. Pending classifications remain visible.

## Validation and key findings

The original frozen PRIMARY comparison assessed 75 deterministically selected pairs: **54 confirmed, six corrected and 15 unresolved**. All six corrections changed target role to relevant comparator while preserving inclusion. Agreement among the 60 adjudicable pairs was 54/60 (90%); this conditional result is not an atlas-wide accuracy estimate.

A separate primary-source follow-up reviewed 12 of the 15 unresolved pairs. Seven later include/cancer-target decisions received documented researcher approval. Five reviewed and three unreviewed original cases remained unresolved. These later decisions do not replace the original validation results or add experimental evidence.

Native vimentin illustrates a target-to-comparator correction based on experimental purpose. The shorter PSA determinant illustrates a later target-role decision that preserves Negative IEDB observations and the paper/export discrepancy. See the [case studies](docs/case_study.md) and [validation account](docs/validation.md).

## Streamlit Evidence Explorer

The single-dataset application supports peptide, antigen, PMID, epitope and title searches; role, CD4, restriction, attribution and resolution filters; assay-level details; source links; and filtered CSV downloads. It verifies the current matrix checksum before loading and shows no substitute if validation fails.

This project is prepared for local use. No public hosted deployment is claimed.

## Local setup and tests

Use Python 3.11 or later in a virtual environment:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
streamlit run portal/app.py
```

GitHub Actions is configured for Python 3.11 and 3.13. The public tests use the distributed matrix and do not require private exports, PDFs or approval records.

## Methods and limitations

The [methods](docs/methods.md) distinguish cohort relevance, biological role, CD4 response, HLA annotation and assay attribution. Restriction categories are structured metadata proxies. PBMC or DTH responses are not purified-CD4 proof; linked observations must not be reassigned to exact peptides. Natural presentation, native tumour recognition and vaccine efficacy are not established atlas-wide. Binding remains a separate research arm.

The matrix is a note-sanitised portfolio derivative: scientific classifications and assay observations are unchanged. Original approval authority remains private. Historical comparisons and full source reconstruction require research inputs outside this distribution. The [release manifest](docs/portfolio_release_manifest.json) records the current checksum and historical seven-decision lineage.

## Licensing and attribution

Original software is covered by [MIT](LICENSES/MIT.txt). Original annotations and documentation are covered by [CC BY 4.0](LICENSES/CC-BY-4.0.txt) to the extent the owner holds the rights. IEDB-derived records and other third-party material retain their applicable rights; these terms do not claim ownership of entire mixed datasets or article content.

See [component scopes](LICENSE), [data licensing](data/LICENSE.md) and [source attribution](ATTRIBUTION.md). Automated checks and note review do not establish complete legal or privacy clearance.
