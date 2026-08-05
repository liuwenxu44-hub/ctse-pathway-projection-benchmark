ctse_fraction_accuracy <- function(estimate, truth) {
  ctse_assert(is.matrix(estimate) && is.matrix(truth), "FRACTION_MATRIX_REQUIRED")
  ctse_assert(identical(dim(estimate), dim(truth)) &&
                identical(dimnames(estimate), dimnames(truth)),
              "FRACTION_SCHEMA_MISMATCH")
  estimate_value <- as.vector(t(estimate))
  truth_value <- as.vector(t(truth))
  data.frame(
    sample_id = rep(rownames(estimate), each = ncol(estimate)),
    cell_type = rep(colnames(estimate), times = nrow(estimate)),
    estimate = estimate_value, truth = truth_value,
    absolute_error = abs(estimate_value - truth_value),
    squared_error = (estimate_value - truth_value)^2,
    stringsAsFactors = FALSE
  )
}

ctse_gene_accuracy <- function(estimate, truth) {
  ctse_assert(is.array(estimate) && is.array(truth) &&
                identical(dim(estimate), dim(truth)) &&
                identical(dimnames(estimate), dimnames(truth)),
              "GENE_ACCURACY_SCHEMA")
  rows <- vector("list", dim(estimate)[[2L]] * dim(estimate)[[3L]])
  index <- 0L
  for (sample in dimnames(estimate)[[3L]]) {
    for (cell_type in dimnames(estimate)[[2L]]) {
      index <- index + 1L
      x <- estimate[, cell_type, sample]
      y <- truth[, cell_type, sample]
      metric <- ctse_vector_metrics(x, y)
      rows[[index]] <- data.frame(
        sample_id = sample, cell_type = cell_type,
        pearson = metric$pearson, spearman = metric$spearman,
        rmse = metric$rmse, bias = mean(x - y),
        negative_fraction = mean(x < 0), stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

ctse_task_summary <- function(fraction_accuracy, gene_accuracy,
                              projection_comparison, scenario = NA_character_) {
  ctse_assert(is.data.frame(fraction_accuracy) && is.data.frame(gene_accuracy),
              "TASK_SUMMARY_INPUT")
  data.frame(
    scenario = scenario,
    fraction_mae = mean(fraction_accuracy$absolute_error),
    fraction_rmse = sqrt(mean(fraction_accuracy$squared_error)),
    ctse_rmse_median = stats::median(gene_accuracy$rmse),
    ctse_pearson_median = stats::median(gene_accuracy$pearson, na.rm = TRUE),
    signal_projection_mae = mean(projection_comparison$absolute_error),
    signal_direction_accuracy = mean(projection_comparison$direction_correct),
    stringsAsFactors = FALSE
  )
}

ctse_bootstrap_mean_ci <- function(values, seed = 20260729L, resamples = 2000L,
                                   level = 0.95) {
  ctse_assert(is.numeric(values) && length(values) >= 2L && all(is.finite(values)),
              "BOOTSTRAP_INPUT")
  ctse_assert(resamples >= 2L && level > 0 && level < 1, "BOOTSTRAP_OPTIONS")
  ctse_set_rng(seed)
  draws <- replicate(resamples, mean(sample(values, length(values), replace = TRUE)))
  alpha <- (1 - level) / 2
  c(
    estimate = mean(values),
    low = unname(stats::quantile(draws, alpha, names = FALSE)),
    high = unname(stats::quantile(draws, 1 - alpha, names = FALSE)),
    n = length(values)
  )
}

ctse_paired_method_comparison <- function(data, method_a, method_b,
                                          metric = "signal_projection_mae",
                                          seed = 20260729L, resamples = 2000L) {
  required <- c("method", "scenario", "replicate_id", metric, "method_success")
  ctse_assert(is.data.frame(data) && all(required %in% names(data)),
              "PAIRED_COMPARISON_SCHEMA")
  a <- data[data$method == method_a & data$method_success,
            c("scenario", "replicate_id", metric), drop = FALSE]
  b <- data[data$method == method_b & data$method_success,
            c("scenario", "replicate_id", metric), drop = FALSE]
  names(a)[[3L]] <- "value_a"
  names(b)[[3L]] <- "value_b"
  paired <- merge(a, b, by = c("scenario", "replicate_id"))
  ctse_assert(nrow(paired) >= 2L, "PAIRED_COMPARISON_INSUFFICIENT")
  difference <- paired$value_a - paired$value_b
  interval <- ctse_bootstrap_mean_ci(difference, seed, resamples)
  data.frame(
    method_a = method_a, method_b = method_b, n_complete_pairs = nrow(paired),
    mean_difference = interval[["estimate"]], ci_low = interval[["low"]],
    ci_high = interval[["high"]],
    wilcox_p_value = suppressWarnings(stats::wilcox.test(difference, exact = FALSE)$p.value),
    stringsAsFactors = FALSE
  )
}

ctse_external_endpoint_levels <- function(estimate, proxy, sample_map,
                                          proxy_name = c("C1", "C2")) {
  proxy_name <- match.arg(proxy_name)
  ctse_assert(is.array(estimate) && length(dim(estimate)) == 3L,
              "EXTERNAL_ESTIMATE_SCHEMA")
  ctse_assert(is.data.frame(sample_map) &&
                all(c("sample_id", "mixture_id", "replicate_id") %in% names(sample_map)),
              "EXTERNAL_SAMPLE_MAP_SCHEMA")
  sample_map <- sample_map[match(dimnames(estimate)[[3L]], sample_map$sample_id), , drop = FALSE]
  ctse_assert(!anyNA(sample_map$sample_id), "EXTERNAL_SAMPLE_MAP_KEYS")
  genes <- dimnames(estimate)[[1L]]
  cell_types <- dimnames(estimate)[[2L]]
  mixtures <- unique(as.character(sample_map$mixture_id))
  expected <- if (proxy_name == "C1") c(length(genes), length(cell_types)) else
    c(length(genes), length(cell_types), length(mixtures))
  ctse_assert(identical(dim(proxy), expected), "EXTERNAL_PROXY_DIMENSION")
  ctse_assert(identical(dimnames(proxy)[[1L]], genes) &&
                identical(dimnames(proxy)[[2L]], cell_types), "EXTERNAL_PROXY_KEYS")
  if (proxy_name == "C2") {
    ctse_assert(identical(dimnames(proxy)[[3L]], mixtures), "EXTERNAL_PROXY_MIXTURE_KEYS")
  }

  level1 <- list()
  index <- 0L
  for (sample in dimnames(estimate)[[3L]]) {
    map <- sample_map[sample_map$sample_id == sample, , drop = FALSE]
    mixture_index <- match(map$mixture_id, mixtures)
    for (cell_type in cell_types) {
      reference <- if (proxy_name == "C1") proxy[, cell_type] else
        proxy[, cell_type, mixture_index]
      available <- all(is.finite(reference))
      metric <- if (available) ctse_vector_metrics(estimate[, cell_type, sample], reference) else
        data.frame(n = 0L, spearman = NA_real_, pearson = NA_real_, rmse = NA_real_,
                   reference_sd = NA_real_, sne = NA_real_)
      metric_success <- available && all(is.finite(unlist(
        metric[c("spearman", "pearson", "sne")], use.names = FALSE
      )))
      index <- index + 1L
      level1[[index]] <- data.frame(
        proxy = proxy_name, aggregation_level = "LEVEL1_GENE_VECTOR",
        sample_id = sample, mixture = map$mixture_id,
        replicate_id = map$replicate_id, cell_type = cell_type,
        paired_gene_count = metric$n, spearman = metric$spearman,
        pearson = metric$pearson, rmse = metric$rmse, sne = metric$sne,
        status = if (metric_success) "SUCCESS" else if (!available)
          "PROXY_UNAVAILABLE_PRESPECIFIED_ABSTENTION" else
          "ENDPOINT_UNAVAILABLE_CONSTANT_OR_NONFINITE_VECTOR",
        stringsAsFactors = FALSE
      )
    }
  }
  level1 <- do.call(rbind, level1)
  level2 <- list()
  index <- 0L
  for (mixture in mixtures) for (cell_type in cell_types) {
    x <- level1[level1$mixture == mixture & level1$cell_type == cell_type, , drop = FALSE]
    keep <- x$status == "SUCCESS"
    index <- index + 1L
    level2[[index]] <- data.frame(
      proxy = proxy_name, aggregation_level = "LEVEL2_CELLTYPE_MIXTURE_MEDIAN",
      sample_id = "", mixture = mixture, replicate_id = NA_integer_,
      cell_type = cell_type, paired_gene_count = if (any(keep)) length(genes) else 0L,
      spearman = if (any(keep)) stats::median(x$spearman[keep]) else NA_real_,
      pearson = if (any(keep)) stats::median(x$pearson[keep]) else NA_real_,
      rmse = NA_real_, sne = if (any(keep)) stats::median(x$sne[keep]) else NA_real_,
      status = if (any(keep)) "SUCCESS" else "PROXY_UNAVAILABLE_PRESPECIFIED_ABSTENTION",
      stringsAsFactors = FALSE
    )
  }
  level2 <- do.call(rbind, level2)
  level3 <- lapply(mixtures, function(mixture) {
    x <- level2[level2$mixture == mixture, , drop = FALSE]
    keep <- x$status == "SUCCESS"
    ctse_assert(any(keep), "EXTERNAL_MIXTURE_NO_AVAILABLE_CELLTYPE", mixture)
    data.frame(
      proxy = proxy_name, aggregation_level = "LEVEL3_MIXTURE_MEDIAN",
      sample_id = "", mixture = mixture, replicate_id = NA_integer_,
      cell_type = "ALL_AVAILABLE_CELL_TYPES", paired_gene_count = length(genes),
      spearman = stats::median(x$spearman[keep]),
      pearson = stats::median(x$pearson[keep]), rmse = NA_real_,
      sne = stats::median(x$sne[keep]), status = "SUCCESS",
      stringsAsFactors = FALSE
    )
  })
  rbind(level1, level2, do.call(rbind, level3))
}
