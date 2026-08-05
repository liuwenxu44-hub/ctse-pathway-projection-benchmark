# Usage guide

The repository exposes two base-R source files:

- `R/ctse_benchmark.R`: validation, projection, vector metrics, structured task execution and atomic output;
- `R/adapters.R`: deterministic conversion of common array/list layouts to `gene × cell_type × sample`.

Source both files from the repository root before use.

## 1. Prepare the inputs

Provide the following objects from your own analysis environment. No input data are bundled with this repository.

### CTSE estimate

A numeric tensor with dimensions:

```text
gene x cell_type x sample
```

Required identifiers:

- unique gene identifiers;
- unique cell-type labels;
- unique sample identifiers.

### Sample metadata

A table with one row per sample and at least:

```text
sample_id, group
```

The `sample_id` values must match the CTSE tensor exactly. The `group` column defines the prespecified comparison.

### Pathway definition

A table with one row per gene and at least:

```text
gene_id, weight
```

`weight` must be finite and numeric, with at least one non-zero value. Positive and negative weights encode the prespecified direction of the pathway contrast. The pathway definition must be fixed before evaluating method outputs.

### Optional simulation truth

For simulations, the truth tensor must use the same dimensions and identifiers as the estimate. It is used only for prespecified benchmark endpoints.

## 2. Convert method output

The benchmark does not call or tune a deconvolution method. Supply its completed output and convert it deterministically:

```r
# Already canonical: gene × cell type × sample
x <- canonicalize_ctse_array(raw_array, genes, cell_types, samples)

# Source × gene × sample, as returned by some tensor methods
x <- canonicalize_source_gene_sample(raw_array, genes, cell_types, samples)

# Named list of gene × sample matrices, one per cell type
x <- canonicalize_celltype_matrices(raw_list, genes, cell_types, samples)
```

Adapters validate dimensions, identifier order, duplicates and finite values. They never fill, impute, filter or rescale genes.

## 3. Validate before scoring

Reject a task with a structured reason if any of the following checks fails:

- duplicated or missing identifiers;
- non-finite numeric values;
- missing samples or cell types;
- incompatible tensor dimensions;
- no overlap with the fixed pathway gene set;
- missing comparison groups;
- an adapter output that cannot be mapped to the common tensor contract.

`run_ctse_benchmark()` converts validation or numerical errors into a structured terminal record. Do not impute failed outputs and do not remove failures from the task denominator.

## 4. Compute the signed projection

For each sample and cell type:

1. intersect the CTSE gene identifiers with the fixed pathway definition;
2. preserve the fixed pathway weights;
3. take the weighted pathway dot product and divide it by the pathway-weight norm and the full-gene contrast L2 norm;
4. keep sample and cell-type identifiers attached to every value.

Apply the prespecified group contrast to these projection values. Use the same transformation, contrast definition, and missing-value policy for every method.

Executable interface:

```r
result <- run_ctse_benchmark(
  estimate = x,
  sample_metadata = sample_metadata,
  pathway = pathway_definition,
  group0 = "control",
  group1 = "case",
  truth = optional_truth_tensor,
  direction_tolerance = 0.05,
  l2_tolerance = 1e-15
)

if (identical(result$status, "SUCCESS")) {
  print(result$projection)
  print(result$comparison) # present when truth is supplied
} else {
  print(result$failure)
}
```

For an estimated-versus-reference gene vector, `ctse_vector_metrics()` returns finite-pair count, Spearman correlation, Pearson correlation, RMSE, reference standard deviation and scale-normalized error:

```text
SNE = RMSE(estimate, reference) / max(SD(reference), 1e-8)
```

## 5. Record the result

The minimum task record should contain:

```text
task_id
method
scenario
replicate_id
status
started_at
finished_at
elapsed_seconds
failure_class
failure_message
output_checksum
```

Allowed terminal states should distinguish at least:

- success;
- structured method or input failure;
- infrastructure failure;
- timeout or interruption.

Use `write_ctse_benchmark_result(result, output_dir, task_id)` to write one RDS result, one SHA-256 checksum record and exactly one terminal marker. A success marker is created only after validation succeeds. Structured failures receive a separate failure marker and are never converted to success.

## 6. Compare methods

Compare methods only on endpoints fixed before the run. Report the complete task denominator, all failure classes, and the same aggregation rule for every method. Do not tune pathway definitions, contrasts, or exclusion rules after inspecting method performance.

## 7. What is not supplied here

This guide does not provide or expose:

- manuscripts or submission files;
- original, processed, or source data;
- benchmark result values or rankings;
- figures or tables;
- frozen project checkpoints;
- server paths, package-session dumps or machine information;
- author, affiliation, funding, or contribution metadata.

The MIT licence applies to this repository's code. It does not alter the terms of third-party datasets or third-party method packages.
