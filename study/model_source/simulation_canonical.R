# Exact frozen helpers. Source R/ctse_benchmark.R first.
q2_assert <- bp_assert <- ctse_assert
bp_stop <- function(code,detail="") stop(paste(code,detail),call.=FALSE)
q2_canonicalize_unico <- function(raw, genes, celltypes, samples) {
  q2_assert(is.array(raw) && identical(dim(raw), c(4L, 2000L, 40L)), "FORMAL_UNICO_RAW_DIMENSION")
  q2_assert(identical(dimnames(raw)[[1L]], celltypes) && identical(dimnames(raw)[[2L]], genes) && identical(dimnames(raw)[[3L]], samples), "FORMAL_UNICO_RAW_KEYS")
  aperm(raw, c(2L, 1L, 3L))
}
q2_canonicalize_tca <- function(raw, genes, celltypes, samples) {
  q2_assert(is.list(raw) && length(raw) == 4L, "FORMAL_TCA_RAW_LIST")
  out <- array(NA_real_, dim = c(2000L, 4L, 40L), dimnames = list(genes, celltypes, samples))
  for (k in seq_len(4L)) {
    m <- raw[[k]]
    q2_assert(is.matrix(m) && identical(dim(m), c(2000L, 40L)) && identical(rownames(m), genes) && identical(colnames(m), samples), "FORMAL_TCA_RAW_KEYS")
    out[, k, ] <- m
  }
  out
}
q2_validate_tensor <- function(x, genes, celltypes, samples) {
  q2_assert(is.array(x) && identical(dim(x), c(2000L, 4L, 40L)), "FORMAL_TENSOR_DIMENSION")
  q2_assert(identical(dimnames(x)[[1L]], genes) && identical(dimnames(x)[[2L]], celltypes) && identical(dimnames(x)[[3L]], samples), "FORMAL_TENSOR_KEYS")
  q2_assert(all(is.finite(x)), "FORMAL_TENSOR_NONFINITE")
}
bp_coerce_public_matrix <- function(x, sample_ids, feature_ids, code) {
  m <- as.matrix(x)
  bp_assert(is.numeric(m) && !is.null(rownames(m)) && !is.null(colnames(m)), code)
  if (identical(rownames(m), sample_ids) && identical(colnames(m), feature_ids)) return(m)
  if (identical(rownames(m), feature_ids) && identical(colnames(m), sample_ids)) return(t(m))
  bp_stop(code)
}

bp_coerce_public_expression_matrix <- function(x, sample_ids, genes) {
  m <- as.matrix(x)
  if (!is.numeric(m) || is.null(rownames(m)) || is.null(colnames(m))) bp_stop("BAYESPRISM_TYPE_EXPRESSION_UNAVAILABLE")
  if (setequal(rownames(m), sample_ids) && !setequal(colnames(m), genes)) bp_stop("BAYESPRISM_GENE_UNIVERSE_CHANGED")
  if (setequal(colnames(m), sample_ids) && !setequal(rownames(m), genes)) bp_stop("BAYESPRISM_GENE_UNIVERSE_CHANGED")
  if (identical(rownames(m), sample_ids) && identical(colnames(m), genes)) return(m)
  if (identical(rownames(m), genes) && identical(colnames(m), sample_ids)) return(t(m))
  bp_stop("BAYESPRISM_TYPE_EXPRESSION_UNAVAILABLE")
}

bp_normalize_type_expression <- function(expression, sample_ids, celltypes, genes) {
  bp_assert(identical(dim(expression), c(length(genes), length(celltypes), length(sample_ids))), "BAYESPRISM_TYPE_EXPRESSION_DIMENSION_MISMATCH")
  bp_assert(identical(dimnames(expression)[[1L]], genes) && identical(dimnames(expression)[[2L]], celltypes) && identical(dimnames(expression)[[3L]], sample_ids), "BAYESPRISM_TYPE_EXPRESSION_IDENTIFIER_MISMATCH")
  bp_assert(all(is.finite(expression)) && all(expression >= 0), "BAYESPRISM_TYPE_EXPRESSION_INVALID")
  raw_sums <- apply(expression, c(2L, 3L), sum)
  bp_assert(all(is.finite(raw_sums)) && all(raw_sums > 0), "BAYESPRISM_TYPE_EXPRESSION_ZERO_PROFILE")
  normalized <- expression
  for (sample_index in seq_along(sample_ids)) {
    for (type_index in seq_along(celltypes)) {
      normalized[, type_index, sample_index] <- normalized[, type_index, sample_index] / raw_sums[type_index, sample_index]
    }
  }
  bp_assert(max(abs(apply(normalized, c(2L, 3L), sum) - 1)) < 1e-14, "BAYESPRISM_TYPE_EXPRESSION_NORMALIZATION_FAILURE")
  list(normalized = normalized, raw_sums = raw_sums)
}

