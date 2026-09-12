# Model source, scale and execution boundary

The endpoint runner never sources model entrypoints. New public wrappers default
to allow_new_fit=FALSE. Source/commit maps are in study/model_source and
study/SOURCE_LINEAGE.json. Prepared task/seed registries are in release data.

Simulation Unico uses its recorded package defaults, C1/C2=NULL, parallel=FALSE.
Simulation TCA uses max_iters=10, refit_W=FALSE and tensor(scale=FALSE).
Both retain signed CPM including negative values. BayesPrism uses the count/
reference API and its explicitly declared simulation unit-sum type profile.
CellBench has its own separately recorded parameters, not these simulation defaults.

## EPIC two-stage scale

Simulation computes L=log2(CPM_full_2000+1) once before prior/returned-gene filtering.
First stage receives L with log2transf=FALSE; second receives the same response
on stage1-returned IDs. No denominator is recomputed after removal. Official
get_prior uses the frozen reference and original synthetic donor labels.

CellBench's preserved first stage internally normalized the full 13,868-gene
CPM response and applied log2(x+1). Corrected stage2 receives that identical stored
response and reuses the posterior. The interface expects the verified persisted
log response; it never guesses scale from magnitude or transforms it again.

Frozen EPIC: stage1 nu=50, nitt=1300, burnin=300, thin=1; stage2 nstop=1,
delta=0.1, nu0=50, nu1=50, one worker. The installed y=NULL branch does not
forward the outer seed; effective per-gene stage2 seed123 is preserved/disclosed.
The separate bounded diagnostic instrumentation records its actual seeds.

Output A is sample-specific conditional log-expression. Canonicalization only
performs approved ID reorder/renaming, not exp, clipping, normalization or filling.
Spearman eligibility does not imply amplitude or projection eligibility.

## What has and has not been tested

The public derived runtime contains no method packages, hence performs no fit.
Source-level interfaces, pinned identities and prepared inputs are supplied for
inspection and intentional future use. They do not certify a complete new
cross-machine stochastic refit. Third-party fitting dependencies require their
own installation and licenses. The old legacy wrappers are not the current
simulation configuration.

Prepared state0 inputs contain the frozen W estimated for that regime; pooled
inputs contain original W. This reproduction never reruns NNLS. Practical and
oracle data remain distinct; never give evaluation bundles to a fitting function.

repro/archive.py provides close/fsync → hash verification → atomic no-replace
rename → read-only publication. Any later deliberately requested fit must be a
new generation, not a replacement or performance-selected retry of this release.
