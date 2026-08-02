# Usage guide

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

`weight` must be numeric. Positive and negative weights encode the prespecified direction of the pathway contrast. The pathway definition must be fixed before evaluating method outputs.

### Optional simulation truth

For simulations, the truth tensor must use the same dimensions and identifiers as the estimate. It is used only for prespecified benchmark endpoints.

## 2. Validate before scoring

Reject a task with a structured reason if any of the following checks fails:

- duplicated or missing identifiers;
- non-finite numeric values;
- missing samples or cell types;
- incompatible tensor dimensions;
- insufficient overlap with the fixed pathway gene set;
- missing comparison groups;
- an adapter output that cannot be mapped to the common tensor contract.

Do not impute failed method outputs and do not remove failures from the task denominator.

## 3. Compute the signed projection

For each sample and cell type:

1. intersect the CTSE gene identifiers with the fixed pathway definition;
2. preserve the fixed pathway weights;
3. combine the CTSE values with those weights to obtain one signed projection value;
4. keep sample and cell-type identifiers attached to every value.

Apply the prespecified group contrast to these projection values. Use the same transformation, contrast definition, and missing-value policy for every method.

Conceptual pseudocode:

```text
for each method output:
    estimate = adapt_to_gene_celltype_sample_tensor(method_output)
    validate(estimate, sample_metadata, pathway_definition)

    for each cell_type:
        for each sample:
            projection = signed_weighted_projection(
                estimate[pathway_genes, cell_type, sample],
                fixed_pathway_weights
            )

        endpoint = prespecified_group_contrast(projection, sample_metadata)
        save(endpoint)

    save_success_marker_and_checksum()
```

## 4. Record the result

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

Write outputs atomically, then create the success marker and checksum only after validation succeeds.

## 5. Compare methods

Compare methods only on endpoints fixed before the run. Report the complete task denominator, all failure classes, and the same aggregation rule for every method. Do not tune pathway definitions, contrasts, or exclusion rules after inspecting method performance.

## 6. What is not supplied here

This guide does not provide or expose:

- manuscripts or submission files;
- original, processed, or source data;
- benchmark result values or rankings;
- figures or tables;
- frozen project checkpoints;
- server paths, package-session dumps, or machine information;
- author, affiliation, funding, or contribution metadata.
