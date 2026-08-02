# CTSE pathway-projection benchmark

This repository is a documentation-only guide for applying a downstream pathway-projection benchmark to cell-type-specific expression (CTSE) estimates.

It intentionally contains no manuscript, study data, source data, figures, result tables, checkpoints, analysis environment records, machine paths, author metadata, or contribution statements.

## What the benchmark does

The benchmark evaluates whether a CTSE method preserves a predefined signed pathway contrast after the method output has been converted to a common tensor representation.

The workflow has four stages:

1. Prepare CTSE estimates in a common gene-by-cell-type-by-sample layout.
2. Apply the same predefined pathway gene set and sample-group contrast to every method.
3. Calculate the signed pathway projection with a method-independent adapter.
4. Record endpoint values and structured failures without replacing or silently dropping failed tasks.

## How to use it

See [docs/USAGE.md](docs/USAGE.md) for the input contract, validation rules, benchmark steps, and the minimum output record.

## Repository boundary

This repository documents the benchmark interface only. Data, analyses, manuscript materials, scientific results, and submission assets are maintained outside this repository and are not distributed here.
