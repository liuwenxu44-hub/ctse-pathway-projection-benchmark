# Executed validation evidence — v0.3.0

This is an engineering reproduction record, not evidence of a new model fit or
journal readiness. Date: 2026-09-10. Manuscripts and figures were not changed.

## Real frozen-data rebuilds

- Two separate, network-disabled, non-root container runs completed.
- 129 scientific tables per run (simulation: 25; reference: 55; diagnostics: 29; cellbench: 20).
- All 129 A/B tables are byte-identical.
- Frozen-table reconciliation passed for both builds: counts, keys, status,
  eligibility and missingness are exact; only predeclared continuous fields may
  use the unchanged 1e-12 absolute plus 1e-12 relative roundoff tolerance.
- Largest absolute difference from the historical frozen tables:
  1.0436096431476471e-14. Original tables were not edited.
- 2526 payload files and the separate input manifest
  were verified before/after reconstruction and after archive extraction.
- All 2321 public RDS files passed recursive metadata/content
  privacy inspection. The text/source scans also passed.

The runner can copy explicitly identified inherited historical tables. This
does not claim that every table was freshly computed from a model. The data
dictionary and reproduction guide identify the inherited products.

## Artificial tests and environment

Eight actual test entrypoints passed, including 30
semantic checks, 16 comparison/archive/privacy
tool checks, 6 atomic-publication checks, original
core tests, reference-summary fixtures and source-only model-interface parsing.
Intentional unit mismatches, gene-key changes, silent missing-value filtering,
oracle/proxy access in practical fitting, post-subset CPM denominators, duplicate
publication, and falsely passing constant/shifted chains are negative cases.

R 4.5.3, 20 exact R sources and six exact Python wheels were installed in the
isolated image. A separate invocation of the public preparation helper downloaded
and verified all 27 pinned R/package/wheel build inputs. No host R libraries,
home directory or method packages were mounted into the execution environment.

## Explicit limits

New model fits = 0; NNLS fits = 0. The tested public route is archived scientific
output to endpoint and summary. Source interfaces, prepared fitting inputs,
configuration/seed identities and corrected scale rules are supplied, but a new
end-to-end four-method refit or cross-machine stochastic identity is not certified.

GitHub Actions is not activated because the publishing credential lacks workflow
permission. The optional repro/ci-template.yml is not an executed CI workflow.
The attached evidence describes the real independent containers, not a hosted
green badge. This limitation is not hidden by a generic success claim.

See RELEASE_VERIFICATION.json for the exact image ID and executed scientific
source hashes; the TSV audits expose every table and comparison status.
