# Validation

## Frozen PRIMARY comparison

The full pre-specified sample comprised **75 deterministically selected publication–epitope pairs**. Of these, **60 were adjudicable** and 15 remained unresolved.

| Frozen PRIMARY outcome | Pairs | Fraction of all 75 |
| --- | ---: | ---: |
| Machine status/role confirmed | 54 | 72% |
| Biological role corrected | 6 | 8% |
| Unresolved pair-level ground truth | 15 | 20% |
| Total | 75 | 100% |

Among the 60 adjudicable pairs, conditional role agreement was **54/60 (90%)**. All six corrections changed `cancer_target` to `relevant_comparator`; inclusion remained unchanged. No confirmed inclusion-status discordance was observed among those 60 pairs. Unresolved cases are not counted as agreement or excluded silently from the full-sample result.

The 90% conditional figure is **not an atlas-wide model accuracy estimate**. Selection was deterministic, non-random and risk-aware; peptides were clustered within references, so observations are not independent. Fifteen pairs lacked defensible final ground truth, and review was not demonstrated to be blinded. These are literature-based curation comparisons, not new laboratory verification.

Experimental purpose determines target/comparator role; source naming, title overlap and response polarity are insufficient. The [native vimentin case](case_study.md) illustrates one original PRIMARY correction. The PSA case illustrates a later post-validation resolution, not a second example of those six corrections.

Machine-assisted primary-literature assessment produced evidence summaries and proposals. Canonical researcher approval and rule-validation records are separate authority. Historical individual-review flags are not comprehensive sign-off lists. Supplementary/enrichment observations were separate from the frozen PRIMARY denominator.

## Post-validation follow-up

A separate audit revisited 12 of the original 15 unresolved pairs using primary-source evidence. Seven proposed include/cancer_target decisions received documented researcher approval and were propagated by exact key; five reviewed and three unreviewed cases remain unresolved. This follow-up does not replace the frozen 54/6/15 result, change its denominator, or demonstrate improved experimental evidence. Original IEDB outcomes, including Negative values discrepant with a paper, were retained.

## Public reproduction scope

The public suite checks the current matrix hash, schema, unique keys, counts, seven documented decisions, functional/restriction evidence, note regressions and portal interactions. Only Post-validation v1 is distributed. Historical export-to-matrix reconstruction and frozen/current cell comparisons require inputs outside this checkout; they are not reproduced by these tests.

The earlier two-dataset portfolio passed 40 tests locally under Python 3.13 and in GitHub Actions under Python 3.11 and 3.13. Those results do not validate this adapted single-dataset project. Fresh local results for this project are recorded below. Its GitHub Actions workflow is configured for Python 3.11 and 3.13 but has not run remotely. No hosted deployment is claimed.

## Final-only local technical validation

A newly created isolated Python 3.13.14 environment installed the public requirements successfully; dependency checking reported no broken requirements. The adapted suite passed **35 tests, zero failures, zero errors and zero skips**. Tests cover current-data integrity and portal behaviour, including an isolated four-asset runtime. The checksum matched the approved Post-validation portfolio v1 value.

The headless Streamlit server returned HTTP 200 from its health endpoint and was terminated after checking. This verifies server startup; AppTest separately verifies dataset loading and interactions. It is not a full browser usability test. Python 3.11 and 3.13 are configured in CI, but this new project's workflow has not run on GitHub. No deployment is claimed.
