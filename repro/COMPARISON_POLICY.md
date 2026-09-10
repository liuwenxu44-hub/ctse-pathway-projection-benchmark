# Rebuild verification policy (fixed before the clean-environment run)

All registered units, methods, failures, undefined values, gene-support counts,
proxy/arm labels and endpoint definitions must remain present. A new successful
fit is not a substitute for a missing archived result. The rebuild has no fit
entrypoint and does not load method packages.

1. Every published input and expected-result file is SHA-256 verified before use.
2. Byte-identical rebuilt tables are the strongest match and are reported as such.
3. A different row order is allowed only after exact, unique scientific-key
   reconciliation. No row may be added, dropped or duplicated. All categorical
   fields, integer counts and missingness masks must match exactly.
4. Only predeclared continuous numeric fields may use cross-environment tolerance:
   `abs(new - frozen) <= 1e-12 + 1e-12 * abs(frozen)`. Maximum observed differences
   and affected cell counts are reported. This is a numerical transport/rebuild
   check, not an endpoint threshold or a change to published values.
5. Input bundle metadata hashes and relative source paths differ by design after
   privacy-safe transport. These are checked against the new public manifest;
   they are not described as byte-identical to private operational manifests.
6. Original failures and semantic non-comparability remain failures/NA. Any
   change in a label, eligibility flag or finite/undefined mask blocks release.
7. Two new isolated rebuilds must agree under this policy. A match is not evidence
   of universal convergence, methodological superiority or journal readiness.
