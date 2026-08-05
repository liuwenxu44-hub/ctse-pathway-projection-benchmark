# CTSE pathway-projection benchmark

This repository provides the executable, failure-aware code companion for a controlled cell-type-specific expression (CTSE) benchmark and its downstream signed pathway-projection evaluation.

It includes the prespecified simulation generator, method wrappers, canonical adapters and endpoint calculations. It intentionally contains no manuscript, study data, figures, result tables, checkpoints, analysis-environment records, machine paths, author metadata, funding information or contribution statements.

## What the benchmark does

The benchmark evaluates whether a CTSE method preserves a prespecified signed pathway contrast after the method output has been converted to a common tensor representation.

The workflow has four stages:

1. Prepare CTSE estimates in a common gene-by-cell-type-by-sample layout.
2. Apply the same predefined pathway gene set and sample-group contrast to every method.
3. Calculate the signed pathway projection with a method-independent adapter.
4. Record endpoint values and structured failures without replacing, imputing or silently dropping failed tasks.

The study companion additionally provides:

- the six frozen controlled-simulation scenarios;
- deterministic reference and bulk generation;
- the shared NNLS fraction input;
- parameterized BayesPrism, Unico and TCA calls;
- gene, fraction, SNE, external-proxy and paired-bootstrap endpoints.

The projection for pathway weights `w` and a cell-type-specific group contrast `delta` is:

```text
sum(w * delta[pathway genes]) / (sqrt(sum(w^2)) * ||delta||2)
```

The denominator uses the full-gene contrast norm. When all pathway weights are `1`, this reduces to the size-normalized formula used by the accompanying study.

## How to use it

The implementation uses base R only.

```r
source("R/ctse_benchmark.R")
source("R/adapters.R")
source("R/simulation.R")
source("R/method_wrappers.R")
source("R/study_endpoints.R")

result <- run_ctse_benchmark(
  estimate = ctse_tensor,
  sample_metadata = sample_metadata,
  pathway = pathway_definition,
  group0 = "control",
  group1 = "case"
)
```

Run the included checks and synthetic example from the repository root:

```bash
Rscript tests/run_tests.R
Rscript examples/minimal_example.R
Rscript examples/simulation_case.R
```

See [docs/USAGE.md](docs/USAGE.md) for the input contract, validation rules, adapters, output schema and failure policy, and [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for the exact public-code boundary.

## Public data

No public data are redistributed. The experimental-mixture source datasets remain available from NCBI GEO under [GSE220605](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220605) and [GSE220606](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220606), subject to the source repository and original-study terms.

## Licence

The code in this repository is released under the [MIT License](LICENSE). That software licence does not relicense third-party GEO data.

## Repository boundary

This repository contains the study's reusable simulation, method-call and endpoint implementation. Reproducing data-dependent analyses requires user-supplied GEO downloads or canonical CTSE outputs. Data, scientific results, operational infrastructure, manuscript materials and submission assets are not distributed here.
