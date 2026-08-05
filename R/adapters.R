adapter_assert <- function(ok, code) {
  if (!isTRUE(ok)) stop(code, call. = FALSE)
  invisible(TRUE)
}

validate_adapter_keys <- function(values, label) {
  adapter_assert(is.character(values) && length(values) > 0L, paste0(label, "_KEY_TYPE"))
  adapter_assert(!anyNA(values) && all(nzchar(values)), paste0(label, "_KEY_MISSING"))
  adapter_assert(!anyDuplicated(values), paste0(label, "_KEY_DUPLICATED"))
}

canonicalize_ctse_array <- function(raw, genes, cell_types, samples) {
  validate_adapter_keys(genes, "GENE")
  validate_adapter_keys(cell_types, "CELL_TYPE")
  validate_adapter_keys(samples, "SAMPLE")
  adapter_assert(is.array(raw) && length(dim(raw)) == 3L, "RAW_ARRAY_DIMENSION")
  adapter_assert(identical(dim(raw), c(length(genes), length(cell_types), length(samples))),
                 "RAW_ARRAY_SHAPE")
  adapter_assert(!is.null(dimnames(raw)) && identical(dimnames(raw), list(genes, cell_types, samples)),
                 "RAW_ARRAY_KEYS")
  adapter_assert(is.numeric(raw) && all(is.finite(raw)), "RAW_ARRAY_VALUES")
  raw
}

canonicalize_source_gene_sample <- function(raw, genes, cell_types, samples) {
  validate_adapter_keys(genes, "GENE")
  validate_adapter_keys(cell_types, "CELL_TYPE")
  validate_adapter_keys(samples, "SAMPLE")
  adapter_assert(is.array(raw) && length(dim(raw)) == 3L, "SOURCE_GENE_SAMPLE_DIMENSION")
  adapter_assert(identical(dim(raw), c(length(cell_types), length(genes), length(samples))),
                 "SOURCE_GENE_SAMPLE_SHAPE")
  adapter_assert(!is.null(dimnames(raw)) && identical(dimnames(raw), list(cell_types, genes, samples)),
                 "SOURCE_GENE_SAMPLE_KEYS")
  adapter_assert(is.numeric(raw) && all(is.finite(raw)), "SOURCE_GENE_SAMPLE_VALUES")
  aperm(raw, c(2L, 1L, 3L))
}

canonicalize_celltype_matrices <- function(raw, genes, cell_types, samples) {
  validate_adapter_keys(genes, "GENE")
  validate_adapter_keys(cell_types, "CELL_TYPE")
  validate_adapter_keys(samples, "SAMPLE")
  adapter_assert(is.list(raw) && length(raw) == length(cell_types), "CELLTYPE_LIST_LENGTH")
  adapter_assert(!is.null(names(raw)) && identical(names(raw), cell_types), "CELLTYPE_LIST_KEYS")
  output <- array(
    NA_real_,
    dim = c(length(genes), length(cell_types), length(samples)),
    dimnames = list(genes, cell_types, samples)
  )
  for (index in seq_along(cell_types)) {
    matrix_value <- raw[[index]]
    adapter_assert(is.matrix(matrix_value) && is.numeric(matrix_value), "CELLTYPE_MATRIX_TYPE")
    adapter_assert(identical(dim(matrix_value), c(length(genes), length(samples))),
                   "CELLTYPE_MATRIX_SHAPE")
    adapter_assert(identical(rownames(matrix_value), genes) && identical(colnames(matrix_value), samples),
                   "CELLTYPE_MATRIX_KEYS")
    adapter_assert(all(is.finite(matrix_value)), "CELLTYPE_MATRIX_VALUES")
    output[, index, ] <- matrix_value
  }
  output
}
