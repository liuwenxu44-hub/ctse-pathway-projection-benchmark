ctse_abort <- function(code, detail = "") {
  message <- if (nzchar(detail)) paste(code, detail, sep = ": ") else code
  condition <- structure(
    list(message = message, call = NULL, code = code, detail = detail),
    class = c("ctse_benchmark_error", "error", "condition")
  )
  stop(condition)
}

ctse_assert <- function(ok, code, detail = "") {
  if (!isTRUE(ok)) ctse_abort(code, detail)
  invisible(TRUE)
}

ctse_set_rng <- function(seed) {
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(as.integer(seed))
  invisible(seed)
}

ctse_validate_identifiers <- function(values, axis) {
  ctse_assert(is.character(values), "IDENTIFIER_TYPE", axis)
  ctse_assert(length(values) > 0L, "IDENTIFIER_EMPTY_AXIS", axis)
  ctse_assert(!anyNA(values) && all(nzchar(values)), "IDENTIFIER_MISSING", axis)
  ctse_assert(!anyDuplicated(values), "IDENTIFIER_DUPLICATED", axis)
  invisible(TRUE)
}

validate_ctse_inputs <- function(estimate, sample_metadata, pathway) {
  ctse_assert(is.array(estimate) && length(dim(estimate)) == 3L,
              "CTSE_TENSOR_DIMENSION", "expected gene x cell_type x sample")
  ctse_assert(is.numeric(estimate), "CTSE_TENSOR_TYPE", "estimate must be numeric")
  ctse_assert(all(is.finite(estimate)), "CTSE_TENSOR_NONFINITE")
  keys <- dimnames(estimate)
  ctse_assert(length(keys) == 3L && !any(vapply(keys, is.null, logical(1L))),
              "CTSE_TENSOR_DIMNAMES")
  ctse_validate_identifiers(keys[[1L]], "gene")
  ctse_validate_identifiers(keys[[2L]], "cell_type")
  ctse_validate_identifiers(keys[[3L]], "sample")

  ctse_assert(is.data.frame(sample_metadata), "SAMPLE_METADATA_TYPE")
  ctse_assert(all(c("sample_id", "group") %in% names(sample_metadata)),
              "SAMPLE_METADATA_COLUMNS")
  ctse_validate_identifiers(as.character(sample_metadata$sample_id), "sample_metadata.sample_id")
  ctse_assert(!anyNA(sample_metadata$group) && all(nzchar(as.character(sample_metadata$group))),
              "SAMPLE_GROUP_MISSING")
  ctse_assert(setequal(keys[[3L]], as.character(sample_metadata$sample_id)),
              "SAMPLE_IDENTIFIER_MISMATCH")

  ctse_assert(is.data.frame(pathway), "PATHWAY_TYPE")
  ctse_assert(all(c("gene_id", "weight") %in% names(pathway)), "PATHWAY_COLUMNS")
  ctse_validate_identifiers(as.character(pathway$gene_id), "pathway.gene_id")
  ctse_assert(is.numeric(pathway$weight) && all(is.finite(pathway$weight)),
              "PATHWAY_WEIGHT_INVALID")
  ctse_assert(any(pathway$weight != 0), "PATHWAY_WEIGHT_ALL_ZERO")
  ctse_assert(any(pathway$gene_id %in% keys[[1L]]), "PATHWAY_NO_GENE_OVERLAP")

  invisible(TRUE)
}

prepare_ctse_inputs <- function(estimate, sample_metadata, pathway) {
  validate_ctse_inputs(estimate, sample_metadata, pathway)
  samples <- dimnames(estimate)[[3L]]
  sample_metadata <- sample_metadata[match(samples, sample_metadata$sample_id), , drop = FALSE]
  pathway <- pathway[pathway$gene_id %in% dimnames(estimate)[[1L]], , drop = FALSE]
  pathway <- pathway[match(dimnames(estimate)[[1L]], pathway$gene_id, nomatch = 0L), , drop = FALSE]
  list(estimate = estimate, sample_metadata = sample_metadata, pathway = pathway)
}

signed_pathway_projection <- function(estimate, sample_metadata, pathway,
                                      group0, group1, l2_tolerance = 1e-15) {
  input <- prepare_ctse_inputs(estimate, sample_metadata, pathway)
  estimate <- input$estimate
  metadata <- input$sample_metadata
  pathway <- input$pathway

  ctse_assert(is.character(group0) && length(group0) == 1L && nzchar(group0),
              "GROUP0_INVALID")
  ctse_assert(is.character(group1) && length(group1) == 1L && nzchar(group1),
              "GROUP1_INVALID")
  ctse_assert(!identical(group0, group1), "GROUPS_IDENTICAL")
  ctse_assert(is.numeric(l2_tolerance) && length(l2_tolerance) == 1L &&
                is.finite(l2_tolerance) && l2_tolerance >= 0,
              "L2_TOLERANCE_INVALID")

  idx0 <- which(as.character(metadata$group) == group0)
  idx1 <- which(as.character(metadata$group) == group1)
  ctse_assert(length(idx0) > 0L && length(idx1) > 0L, "GROUP_EMPTY")

  delta <- apply(estimate[, , idx1, drop = FALSE], c(1L, 2L), mean) -
    apply(estimate[, , idx0, drop = FALSE], c(1L, 2L), mean)
  if (is.null(dim(delta))) {
    delta <- matrix(delta, ncol = 1L,
                    dimnames = list(dimnames(estimate)[[1L]], dimnames(estimate)[[2L]]))
  }

  pathway_index <- match(pathway$gene_id, rownames(delta))
  weights <- pathway$weight
  weight_norm <- sqrt(sum(weights^2))
  ctse_assert(is.finite(weight_norm) && weight_norm > 0, "PATHWAY_WEIGHT_NORM_INVALID")

  rows <- lapply(seq_len(ncol(delta)), function(cell_index) {
    contrast <- delta[, cell_index]
    full_l2_norm <- sqrt(sum(contrast^2))
    projection <- if (full_l2_norm <= l2_tolerance) {
      0
    } else {
      sum(weights * contrast[pathway_index]) / (weight_norm * full_l2_norm)
    }
    ctse_assert(is.finite(projection) && abs(projection) <= 1 + 1e-12,
                "PROJECTION_INVALID")
    data.frame(
      cell_type = colnames(delta)[[cell_index]],
      group0 = group0,
      group1 = group1,
      n_group0 = length(idx0),
      n_group1 = length(idx1),
      n_genes = nrow(delta),
      n_pathway_genes = length(pathway_index),
      full_l2_norm = full_l2_norm,
      signed_projection = projection,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

classify_projection_direction <- function(value, tolerance = 0.05) {
  ctse_assert(is.numeric(tolerance) && length(tolerance) == 1L &&
                is.finite(tolerance) && tolerance >= 0,
              "DIRECTION_TOLERANCE_INVALID")
  ifelse(is.na(value), NA_character_,
         ifelse(value < -tolerance, "negative",
                ifelse(value > tolerance, "positive", "zero")))
}

ctse_vector_metrics <- function(estimate, reference, sd_floor = 1e-8) {
  ctse_assert(is.numeric(estimate) && is.numeric(reference), "VECTOR_TYPE_INVALID")
  ctse_assert(length(estimate) == length(reference), "VECTOR_LENGTH_MISMATCH")
  ctse_assert(is.numeric(sd_floor) && length(sd_floor) == 1L &&
                is.finite(sd_floor) && sd_floor > 0, "SD_FLOOR_INVALID")
  keep <- is.finite(estimate) & is.finite(reference)
  n <- sum(keep)
  if (n < 2L) {
    return(data.frame(n = n, spearman = NA_real_, pearson = NA_real_,
                      rmse = NA_real_, reference_sd = NA_real_, sne = NA_real_))
  }
  x <- estimate[keep]
  y <- reference[keep]
  rmse <- sqrt(mean((x - y)^2))
  reference_sd <- stats::sd(y)
  safe_cor <- function(method) {
    if (stats::sd(x) == 0 || reference_sd == 0) NA_real_
    else suppressWarnings(stats::cor(x, y, method = method))
  }
  data.frame(
    n = n,
    spearman = safe_cor("spearman"),
    pearson = safe_cor("pearson"),
    rmse = rmse,
    reference_sd = reference_sd,
    sne = rmse / max(reference_sd, sd_floor)
  )
}

compare_projection_to_truth <- function(estimated_projection, truth_projection,
                                        direction_tolerance = 0.05) {
  required <- c("cell_type", "signed_projection")
  ctse_assert(is.data.frame(estimated_projection) && all(required %in% names(estimated_projection)),
              "ESTIMATED_PROJECTION_SCHEMA")
  ctse_assert(is.data.frame(truth_projection) && all(required %in% names(truth_projection)),
              "TRUTH_PROJECTION_SCHEMA")
  ctse_assert(!anyDuplicated(estimated_projection$cell_type), "ESTIMATED_CELLTYPE_DUPLICATED")
  ctse_assert(!anyDuplicated(truth_projection$cell_type), "TRUTH_CELLTYPE_DUPLICATED")
  ctse_assert(setequal(estimated_projection$cell_type, truth_projection$cell_type),
              "PROJECTION_CELLTYPE_MISMATCH")
  truth_projection <- truth_projection[
    match(estimated_projection$cell_type, truth_projection$cell_type), , drop = FALSE
  ]
  estimated_value <- estimated_projection$signed_projection
  truth_value <- truth_projection$signed_projection
  estimated_direction <- classify_projection_direction(estimated_value, direction_tolerance)
  truth_direction <- classify_projection_direction(truth_value, direction_tolerance)
  data.frame(
    cell_type = estimated_projection$cell_type,
    estimated_projection = estimated_value,
    truth_projection = truth_value,
    absolute_error = abs(estimated_value - truth_value),
    estimated_direction = estimated_direction,
    truth_direction = truth_direction,
    direction_correct = estimated_direction == truth_direction,
    stringsAsFactors = FALSE
  )
}

run_ctse_benchmark <- function(estimate, sample_metadata, pathway, group0, group1,
                               truth = NULL,
                               direction_tolerance = 0.05,
                               l2_tolerance = 1e-15) {
  tryCatch({
    projection <- signed_pathway_projection(
      estimate = estimate,
      sample_metadata = sample_metadata,
      pathway = pathway,
      group0 = group0,
      group1 = group1,
      l2_tolerance = l2_tolerance
    )
    projection$direction <- classify_projection_direction(
      projection$signed_projection, direction_tolerance
    )
    comparison <- NULL
    if (!is.null(truth)) {
      ctse_assert(is.array(truth) && identical(dim(truth), dim(estimate)) &&
                    identical(dimnames(truth), dimnames(estimate)),
                  "TRUTH_TENSOR_SCHEMA")
      ctse_assert(is.numeric(truth) && all(is.finite(truth)), "TRUTH_TENSOR_VALUES")
      truth_projection <- signed_pathway_projection(
        estimate = truth,
        sample_metadata = sample_metadata,
        pathway = pathway,
        group0 = group0,
        group1 = group1,
        l2_tolerance = l2_tolerance
      )
      comparison <- compare_projection_to_truth(
        projection, truth_projection, direction_tolerance
      )
    }
    list(
      schema_version = "CTSE_PATHWAY_BENCHMARK_V1",
      status = "SUCCESS",
      projection = projection,
      comparison = comparison,
      failure = NULL
    )
  }, ctse_benchmark_error = function(error) {
    list(
      schema_version = "CTSE_PATHWAY_BENCHMARK_V1",
      status = "STRUCTURED_FAILURE",
      projection = NULL,
      comparison = NULL,
      failure = data.frame(
        failure_class = "INPUT_OR_ENDPOINT_SCHEMA_FAILURE",
        failure_code = error$code,
        failure_detail = error$detail,
        stringsAsFactors = FALSE
      )
    )
  })
}

ctse_atomic_write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), ".tmp_"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(lines, temporary, useBytes = TRUE)
  ctse_assert(file.rename(temporary, path), "ATOMIC_PROMOTION_FAILED", basename(path))
  invisible(path)
}

ctse_atomic_write_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), ".tmp_"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3L, compress = TRUE)
  ctse_assert(file.rename(temporary, path), "ATOMIC_PROMOTION_FAILED", basename(path))
  invisible(path)
}

write_ctse_benchmark_result <- function(result, output_dir, task_id) {
  ctse_assert(is.list(result) && result$status %in% c("SUCCESS", "STRUCTURED_FAILURE"),
              "RESULT_SCHEMA_INVALID")
  ctse_assert(is.character(task_id) && length(task_id) == 1L &&
                grepl("^[A-Za-z0-9_.-]+$", task_id), "TASK_ID_INVALID")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(output_dir, paste0(task_id, ".rds"))
  checksum_path <- file.path(output_dir, paste0(task_id, ".sha256"))
  marker_path <- file.path(output_dir, paste0(task_id, ".", result$status))
  ctse_atomic_write_rds(result, result_path)
  checksum <- unname(tools::sha256sum(result_path))
  ctse_atomic_write_lines(paste(checksum, basename(result_path)), checksum_path)
  ctse_atomic_write_lines(c(paste0("status=", result$status), paste0("sha256=", checksum)), marker_path)
  list(result = result_path, checksum = checksum_path, marker = marker_path)
}
