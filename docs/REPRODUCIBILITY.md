# Reproducibility boundary

This repository is the public code companion for the controlled CTSE benchmark and its downstream pathway-projection evaluation.

## Included

- The six-scenario, 2,000-gene, four-cell-type controlled simulation design.
- Deterministic reference and bulk generators with the prespecified RNG namespace.
- The shared NNLS fraction estimator used as input by Unico and TCA.
- Parameterized wrappers for BayesPrism 2.2.3, Unico 0.1.0 and TCA 1.2.1.
- Canonical tensor adapters and method-internal returned-gene handling.
- Gene-level, fraction, SNE, signed-projection, paired-bootstrap and external-proxy endpoint functions.
- Synthetic examples and executable checks.

## Deliberately excluded

- GEO source data and derived matrices.
- Method outputs, scientific results, checkpoints and source-data tables.
- Manuscript and submission files, figures and author declarations.
- Server orchestration, installation recovery, process monitoring, logs and environment-specific library paths.
- Usernames, machine paths, network addresses, credentials and internal asset registries.

The exclusions above are not required to understand the algorithms. They prevent redistribution of third-party data, scientific results and private operational metadata. Users supply their own public-data downloads or canonical method outputs under the schemas documented in [USAGE.md](USAGE.md).

## Software versions used in the study

- R 4.5 series
- BayesPrism 2.2.3
- Unico 0.1.0
- TCA 1.2.1
- `nnls` for the shared non-negative least-squares fraction input

The core endpoint and simulation functions use base R. Method wrappers require the corresponding packages and preserve thrown numerical or termination errors as structured failures; they do not tune or retry failed fits.

## Public data

The external experimental-mixture data are publicly available from NCBI GEO as GSE220605 and GSE220606. This repository does not copy or relicense those records. Download and use remain subject to GEO and the original study's terms.
