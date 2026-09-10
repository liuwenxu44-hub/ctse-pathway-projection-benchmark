# Scientific payload dictionary

Paths are relative to the extracted data root. input-manifest.json is the closed
file/SHA inventory. expected/expected-results.json defines table schemas, row
counts, hashes and continuous columns eligible for the fixed roundoff check.

| Directory | Contents | Boundary |
| --- | --- | --- |
| simulation | 120 unit bundles, seeds, original/state0 method tensors, truth, baselines, inherited directions | Evaluation only; not fitting input |
| cellbench | 17 valid canonical objects, proxies, eligibility, contexts, Reactome, inherited fractions | Platform, arm and C1/C2 separate |
| diagnostics | 24 jobs and 1,920 saved chains with schedules, seeds, bounds and warnings | Bounded design only |
| fit-inputs | 240 prepared simulation inputs, 12 reference objects, four corrected CellBench input/posterior sets | Separate from target CTSE truth |
| expected | 129 frozen scientific tables | Verification targets, not optimization targets |

Simulation scenarios: composition_only, intrinsic_only, concordant_mixed,
discordant_mixed, complete_cancellation and null, each with 20 independent units,
40 samples, four cell types and 2,000 genes. Tensor order is gene × type × sample.
Missing methods remain explicit failures. Missing returned genes remain absent.

Pooled references have 4,800 simulated cells including both states; state0
references have 2,400 cells. The comparison changes both state information and
cell count, not a single causal prior-strength factor. Donors are synthetic
design labels. Baselines use reference/bulk, not target cell-type truth.

CellBench umbrella GSE118767 includes RNA-mixture subseries GSE117617 (CEL-seq2)
and GSE117618 (SORT-seq). Input universe is 13,868 Ensembl genes. Corrected EPIC
returns 13,847/13,868 genes for CEL practical/oracle and 13,845/13,868 for SORT.
Unico's frozen execution exclusion ENSG00000105519 remains recorded.
Other returned universes are read from manifests, never guessed or filled.

C1 has four actually usable mixed compositions; C2 has five compositions with
the corresponding type present. Seven planned compositions remain in status
denominators. Mass strata 3.75/7.5/15/30 pg are aggregated before composition
summaries. Technical libraries, genes and pathways are not independent n.

Projection measures direction alignment, not absolute change or significance.
A nonzero label at 0.05 is not automatically a biological false positive.
Read contrast magnitude and baselines alongside direction. EPIC native-log truth
is the log of a frozen population profile; the original linear mixture generator
is not thereby declared log-additive. Linear/native-log rankings stay separate.
