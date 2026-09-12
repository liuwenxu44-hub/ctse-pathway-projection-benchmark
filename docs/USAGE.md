# Reusable R API

The study reproduction uses study/ frozen kernels. R/ contains reusable helpers
and legacy compatibility interfaces. Load R/ctse_benchmark.R before other helpers.

## Explicit units and targets

Amplitude comparison now requires a declared matching unit and mathematical target:

```r
source("R/ctse_benchmark.R")
contract <- ctse_comparison_contract(
  estimate_unit="CPM", reference_unit="CPM",
  estimate_target="conditional_expression", reference_target="conditional_expression")
ctse_vector_metrics(estimate_vector, reference_vector,
                    comparison_contract=contract)
```

This declaration is not evidence by itself: callers must establish its truth.
Wrong units/targets, nonfinite values and mismatched named IDs fail closed.
SNE is RMSE / max(SD(reference), 1e-8), only for a compatible comparison.

run_ctse_benchmark(..., truth=truth_tensor, comparison_contract=contract) also
requires this contract. When truth is not supplied it calculates the declared
projection without a truth-performance comparison.

A group contrast is formed first on the appropriate expression scale; the signed
projection divides its pathway dot product by the weight norm and full supported
contrast norm. It is not a mean of sample-normalized projections. Direction
threshold0.05 is not a significance level.

## Tensor and scale guards

R/adapters.R validates gene × cell_type × sample arrays with unique ordered IDs.
No gene is filled or imputed. R/epic_scale_contract.R offers typed linear input,
one full-support CPM/log transform, and subset operations retaining the original
denominator. The fit-access guard rejects evaluation fields and oracle fields
in practical inputs.

## Legacy code

R/simulation.R exposes the historical generator; it is not run by repro/run.py.
R/method_wrappers.R remains for legacy use and is not the later four-method
simulation contract. See study/model_source/frozen_calls.R and MODEL_EXECUTION.md.
R/study_endpoints.R has reusable legacy aggregation helpers; do not substitute
its three-level convenience route for the later exact CellBench mass-stratum
contract. study/rebuild_cellbench.R supplies that frozen implementation.

## Tests

```bash
Rscript --vanilla tests/run_tests.R
Rscript --vanilla study/tests/test_semantics.R
python3 repro/privacy_guard.py .
```

See README for the independent-container test and real frozen-data reconstruction.
