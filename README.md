# CTSE benchmark: cell-type gene expression and pathway changes

Cell-type-specific expression (CTSE) methods estimate gene expression for individual cell types within a mixed sample. This benchmark asks four separate questions:

1. Does a method recover the overall gene-expression profile?
2. Does it recover changes between samples?
3. Does it recover the direction of pathway changes?
4. How large are those changes?

Agreement on one question does not establish agreement on the others. Results measured on different scales are kept separate rather than combined into one method ranking.

## What is included in v0.3.0?

- Six simulation scenarios, each with 20 units, 2,000 genes, four cell types and 40 samples per unit.
- Saved results for BayesPrism, Unico, TCA and EPIC-unmix under two reference settings, pooled-reference and state0-reference. Reference-only, Bulk-copy and Always-zero are included as controls.
- Corrected CellBench CEL-seq2 and SORT-seq results, with the returned genes and Reactome V97 pathway membership recorded. The two proxy references, C1/C2, and the practical/oracle analysis settings remain separate.
- Saved EPIC diagnostic runs, calculation code, exact dependency versions and tests.

Large numerical inputs and expected result tables are provided as release downloads, not stored directly in Git. The repository does not include manuscripts, figures, execution logs or private machine settings.

For v0.3.0, two separate container runs produced the same 129 scientific tables, byte for byte. The [validation report](verification/VALIDATION_SUMMARY.md) lists what was checked and the limits of those checks. Some identified historical tables are copied unchanged, not newly calculated from model fits.

## Recalculate tables from the saved model outputs

The tested workflow starts with the supplied model outputs and recalculates evaluation measures and summaries. **It does not train the methods or fit non-negative least-squares (NNLS) models.** Training code and prepared inputs are included, but a complete new four-method training run was not tested for this release.

Linux x86-64, Docker and host Python 3.10+ are required. Allow approximately 60 GB
disk space for downloads, inputs and two rebuilds. Run from the repository root:

```bash
python3 repro/download_data.py --destination data-v0.3.0
python3 repro/prepare_runtime.py environment
docker build -t ctse-derived:0.3.0 environment
mkdir -p local-rebuilds
docker run --rm --network none --read-only --cap-drop=ALL \
  --security-opt=no-new-privileges --user "$(id -u):$(id -g)" \
  --tmpfs /tmp:rw,nosuid,size=8g \
  -v "$PWD:/work:ro" -v "$PWD/data-v0.3.0:/inputs:ro" \
  -v "$PWD/local-rebuilds:/outputs" -w /work ctse-derived:0.3.0 \
  python3 repro/run.py --data /inputs --output /outputs/A --jobs 4
```

Repeat with `/outputs/B` to make a second rebuild. The runner checks input-file hashes, reconstructs the tables, compares the results with the expected files, and checks that the inputs have not changed. Numerical tolerances are specified in the [comparison policy](repro/COMPARISON_POLICY.md).

## How to interpret the comparisons

- Gene-profile correlations compare the same available genes within each comparison. Results based on different gene sets are not counted as extra replicates.
- Unico/TCA expression on the linear counts-per-million (CPM) scale is not merged with EPIC's native-log secondary results. Rescaling the numbers does not make these targets equivalent.
- The historical BayesPrism direction results describe their original contribution/profile target. They are not reinterpreted as the size of changes in conditional cell-type expression.
- C1 and C2 are two proxy references, not absolute biological truth. CellBench summaries use mixture composition as the unit after aggregating the predefined RNA-mass groups; technical libraries are not extra independent samples.
- A failed method, insufficient evidence, an undefined correlation and an incompatible comparison are different outcomes. None is replaced with zero.
- The direction cutoff of **0.05 is not a p-value**. A change can have a nonzero direction score while being very small in magnitude.

## Documentation

- [Reproduction coverage and limits](docs/REPRODUCIBILITY.md)
- [Data dictionary](docs/DATA_DICTIONARY.md)
- [Model code, expression scales and comparison targets](docs/MODEL_EXECUTION.md)
- [Reusable R API](docs/USAGE.md)
- [Third-party terms](THIRD_PARTY_NOTICES.md)
- [Changelog](CHANGELOG.md)
- [Recorded tests and file-integrity checks](verification/VALIDATION_SUMMARY.md)

Experimental sources include [CellBench GSE118767](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118767). The earlier GSE220605/GSE220606 case is retained only as a limited historical assessment; it was not rerun as new validation for this release.

Original repository code uses the [MIT License](LICENSE). Third-party methods and datasets keep their own terms. This release has no assigned DOI.
