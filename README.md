# CTSE expression and signed-projection reproducibility companion

Source code and versioned scientific data for a cell-type-specific expression
(CTSE) benchmark. Static gene-profile agreement, sample-specific change recovery,
normalized direction and change magnitude are evaluated separately. This project
does not construct one cross-scale ranking of all methods.

## v0.3.0

- Six frozen scenarios × 20 units; 2,000 genes, four cell types and 40 samples/unit.
- Original pooled-reference and state0-reference results for BayesPrism, Unico,
  TCA and EPIC-unmix; Reference-only, Bulk-copy and Always-zero controls.
- Corrected CellBench CEL-seq2/SORT-seq outputs, separate C1/C2 eligibility,
  practical/oracle arms, returned genes and Reactome V97 membership.
- Saved chains from the bounded EPIC diagnostic design.
- Exact frozen endpoint kernels, scientific manifests, locked dependencies,
  negative tests and clean-container reproduction.

Large numeric inputs and expected tables are release assets, not Git objects.
Manuscripts, figures, logs, machine paths and private configuration are excluded.

Two independent isolated derived rebuilds passed: all 129 scientific tables are
byte-identical between builds. [Executed validation evidence](verification/VALIDATION_SUMMARY.md)
includes table-level hashes, frozen-result reconciliation and explicit test limits.

## Rebuild frozen evidence

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

Repeat with /outputs/B for an independent derived rebuild. The runner verifies
input hashes, reconstructs tables, checks the closed expected inventory, and
rechecks inputs. See the fixed [comparison policy](repro/COMPARISON_POLICY.md).

**No model or NNLS is fitted by this command.** The tested route is frozen
output → derived endpoint/summary. Fitting source and prepared inputs are also
supplied, but new four-method fits were not executed in this engineering release.

## Scientific boundaries

- Four-method gene-vector Spearman uses identical actual support within each
  comparison; different support contexts are never pooled into extra replicates.
- Unico/TCA linear CPM and EPIC native-log secondary endpoints remain separate.
  No clipping, forced exponentiation or arbitrary scaling creates comparability.
- Historical BayesPrism direction results retain their original contribution/
  profile target, not a common conditional-expression amplitude interpretation.
- C1/C2 are separate proxies, not absolute truth. Composition is the CellBench
  descriptive unit after frozen mass-stratum aggregation.
- Failure, abstention, undefined correlation and semantic non-comparability are
  distinct; none is imputed as zero.
- The 0.05 direction threshold is not a p-value. A nonzero normalized direction
  can coexist with very small change magnitude.

## Documentation

- [Reproduction coverage and limits](docs/REPRODUCIBILITY.md)
- [Data dictionary](docs/DATA_DICTIONARY.md)
- [Model source and scale contracts](docs/MODEL_EXECUTION.md)
- [Reusable R API](docs/USAGE.md)
- [Third-party terms](THIRD_PARTY_NOTICES.md)
- [Changelog](CHANGELOG.md)
- [Executed tests and data/source integrity](verification/VALIDATION_SUMMARY.md)

Experimental sources include [CellBench GSE118767](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118767).
The earlier GSE220605/GSE220606 case remains a restricted historical assessment,
not a newly rerun validation layer. Original repository code uses the
[MIT License](LICENSE), which does not relicense third-party methods or datasets.
No DOI is assigned or implied for this release.
