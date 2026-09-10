# Reproduction coverage and limits

The four components are simulation, reference, cellbench and diagnostics.
The closed comparison inventory contains 129 frozen scientific tables: long-form
records, eligibility/coverage, paired summaries, baseline controls and diagnostics.
Historical BayesPrism direction tables and the fraction-library table are inherited
verbatim and verified, not falsely described as new model calculations.

Exports retain numeric arrays, IDs, statuses and supports. Objects were read back
and compared with R identical() on their scientific payload. Private provenance
fields are omitted; changed serialization is not byte identity of the original
entire object. The public manifest identifies the exported generation.

The container installs R 4.5.3 and pinned dependencies from verified upstream
sources. It mounts no host R library or home, and has no network during execution.
Scientific method packages are absent. The fixed comparison policy distinguishes
byte identity, exact values with different order/serialization, and predeclared
roundoff. Counts, statuses, keys, eligibility and missingness must match exactly.

Automated GitHub Actions is not activated: the publishing credential does not
have workflow permission. repro/ci-template.yml is an optional template, not
evidence of a hosted CI run. The release evidence is from the actually executed
independent containers, not a claimed green GitHub badge.

## Fitting source is not a freshly verified fit

study/model_source contains corrected EPIC calls, exact canonical helpers,
recorded method interfaces and package/commit identities. fit-inputs provides
frozen responses, fractions, references and posteriors separately from truth.
No model/NNLS was rerun for this public version. Cross-machine stochastic fit
identity and a complete new download-and-refit workflow are not certified.

The legacy R/method_wrappers.R is retained for compatibility. Its parameters
must not be substituted for later simulation calls. Current simulation and
CellBench configurations are explicitly separated in frozen_calls.R.

## Corrections and failures

The four erroneous historical EPIC CellBench stage2 outputs are excluded from
the valid input collection: they mixed a linear response with a log-scale prior.
They are not evidence of algorithm failure. The prior scientific correction
performed four stage2 corrective fits and reused all four stage1 posteriors.
This engineering release performs zero fits.

Original TCA null r002/r004 failures, state0-reference complete_cancellation r010
failure, returned-gene exclusions and CellBench availability limits are retained.
No unavailable entry is zero-filled. Complete paired subsets do not erase a
nonfailed method's remaining available units.

Manuscript/figure changes, submission readiness, additional models/data,
new endpoint definitions and acceptance predictions are outside this release.
