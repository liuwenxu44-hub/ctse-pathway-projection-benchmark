# Source-specific terms and attribution

Verified 2026-09-10. The MIT license of this repository covers its original code;
it is not a replacement license for third-party datasets, methods or dependencies.
The release does not claim that a GEO accession automatically makes software MIT.

## CellBench

Tian et al., A benchmark for single-cell transcriptome profiling with controlled
experimental mixtures, Nature Methods (2019),
[doi:10.1038/s41592-019-0425-8](https://www.nature.com/articles/s41592-019-0425-8).
Original experimental resources are under GEO GSE118767, including GSE117617
and GSE117618. Source processed resources/code:
[LuyiTian/sc_mixology](https://github.com/LuyiTian/sc_mixology).
We distribute selected numerical derivatives for reproduction, not new human
subject data, and preserve accession/study attribution. The source repository's
[MIT license](https://github.com/LuyiTian/sc_mixology/blob/master/LICENSE)
identifies Copyright (c) 2018 Luyi Tian. Its software license is not asserted as
a universal license for all GEO records. Source terms remain applicable.

## Reactome

The frozen V97 gene-membership data retain their version and identifiers. The
official [Reactome license](https://reactome.org/license) states that database
data and derived data are CC0. The separate illustration/icon license is not
used here; no Reactome figure assets are distributed in this release.

## Scientific methods

Pinned identities are in study/model_source/versions.tsv. Method implementations
are obtained from their official repositories, not relicensed as this software:

- [BayesPrism](https://github.com/Danko-Lab/BayesPrism), version 2.2.3;
  [original article](https://doi.org/10.1038/s43018-022-00356-3).
- [Unico](https://github.com/cozygene/Unico), version 0.1.0;
  [original article](https://doi.org/10.1186/s13059-025-03776-3).
- [TCA](https://github.com/cozygene/TCA), version 1.2.1;
  [original article](https://doi.org/10.1038/s41467-019-11052-9).
- [EPIC-unmix](https://github.com/quansun98/EPICunmix/tree/b512206e3e0bcb9f0720423293094935861fdcde),
  version 0.0.1; [original article](https://doi.org/10.1186/s13059-025-03847-5).
  Its pinned DESCRIPTION declares GPL without a version suffix; the historical
  environment manifest's more specific label is not an independent license grant.

The public fitting files are study-owned adapters/instrumentation around these
APIs. No full third-party method source or installed package directory is bundled.
Follow upstream terms when installing/redistributing method implementations.

## Runtime dependencies

R is obtained from official CRAN. Package source URLs, versions, SHA-256 and
declared licenses are in environment/r-source-lock.json. Python wheel versions
and hashes are in environment/python-wheels.tsv, downloaded from official PyPI.
Their copyright/license metadata remain within their source/wheel distributions.
The Docker base is official Ubuntu 24.04 pinned by digest. Build dependencies
come from its official package repositories; this does not impose MIT on them.

The earlier GSE220605/GSE220606 external case is cited for historical context;
this release does not redistribute its entire raw data archive. No nonexistent
repository address, release DOI or software-paper DOI is supplied.
