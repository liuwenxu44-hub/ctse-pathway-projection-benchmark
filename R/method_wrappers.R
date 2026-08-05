ctse_validate_joint_input <- function(bulk, fractions) {
  ctse_assert(is.matrix(bulk) && is.numeric(bulk), "METHOD_BULK_MATRIX_REQUIRED")
  ctse_assert(is.matrix(fractions) && is.numeric(fractions), "METHOD_FRACTION_MATRIX_REQUIRED")
  ctse_assert(!is.null(rownames(bulk)) && !is.null(colnames(bulk)), "METHOD_BULK_KEYS")
  ctse_assert(!is.null(rownames(fractions)) && !is.null(colnames(fractions)),
              "METHOD_FRACTION_KEYS")
  ctse_assert(identical(colnames(bulk), rownames(fractions)), "METHOD_SAMPLE_AXIS_MISMATCH")
  ctse_assert(all(is.finite(bulk)) && all(is.finite(fractions)), "METHOD_INPUT_NONFINITE")
  ctse_assert(all(fractions >= 0) && max(abs(rowSums(fractions) - 1)) < 1e-8,
              "METHOD_FRACTION_INVALID")
  invisible(TRUE)
}

ctse_require_package_version <- function(package, version) {
  ctse_assert(requireNamespace(package, quietly = TRUE),
              paste0(toupper(package), "_PACKAGE_REQUIRED"))
  observed <- as.character(utils::packageVersion(package))
  ctse_assert(identical(observed, version), paste0(toupper(package), "_VERSION_MISMATCH"),
              paste0("expected ", version, "; observed ", observed))
  invisible(TRUE)
}

ctse_method_failure <- function(method, error) {
  list(
    status = "STRUCTURED_FAILURE",
    method = method,
    canonical = NULL,
    native = NULL,
    failure = data.frame(
      failure_class = "METHOD_TERMINATION_OR_NUMERICAL_FAILURE",
      message = conditionMessage(error), stringsAsFactors = FALSE
    )
  )
}

run_unico_joint <- function(bulk, fractions, seed = 272000001L) {
  ctse_validate_joint_input(bulk, fractions)
  ctse_require_package_version("Unico", "0.1.0")
  tryCatch({
    ctse_set_rng(seed)
    fit <- Unico::Unico(
      bulk, fractions, C1 = NULL, C2 = NULL, fit_tau = FALSE,
      mean_penalty = 0, var_penalty = 0.01, covar_penalty = 0.01,
      mean_max_iterations = 2, var_max_iterations = 3,
      nloptr_opts_algorithm = "NLOPT_LN_COBYLA", max_stds = 2,
      init_weight = "default", max_u = 1, max_v = 1,
      parallel = FALSE, num_cores = 1L, log_file = NULL,
      verbose = FALSE, debug = FALSE
    )
    native <- Unico::tensor(
      bulk, fractions, C1 = NULL, C2 = NULL, Unico.mdl = fit,
      parallel = FALSE, num_cores = 1L, log_file = NULL,
      verbose = FALSE, debug = FALSE
    )
    canonical <- canonicalize_source_gene_sample(
      native, rownames(bulk), colnames(fractions), colnames(bulk)
    )
    list(status = "SUCCESS", method = "Unico", canonical = canonical,
         native = native, fit = fit, failure = NULL)
  }, error = function(error) ctse_method_failure("Unico", error))
}

run_tca_joint <- function(bulk, fractions, seed = 273000001L) {
  ctse_validate_joint_input(bulk, fractions)
  ctse_require_package_version("TCA", "1.2.1")
  tryCatch({
    ctse_set_rng(seed)
    fit <- TCA::tca(
      bulk, fractions, C1 = NULL, C1.map = NULL, C2 = NULL,
      refit_W = FALSE, refit_W.features = NULL, refit_W.sparsity = 500,
      refit_W.sd_threshold = 0.02, tau = NULL, vars.mle = FALSE,
      constrain_mu = FALSE, parallel = FALSE, num_cores = 1L,
      max_iters = 10L, log_file = NULL, debug = FALSE, verbose = FALSE
    )
    native <- TCA::tensor(
      bulk, tca.mdl = fit, scale = FALSE, parallel = FALSE,
      num_cores = 1L, log_file = NULL, debug = FALSE, verbose = FALSE
    )
    if (is.null(names(native))) names(native) <- colnames(fractions)
    if (setequal(names(native), colnames(fractions))) {
      native <- native[colnames(fractions)]
    }
    canonical <- canonicalize_celltype_matrices(
      native, rownames(bulk), colnames(fractions), colnames(bulk)
    )
    list(status = "SUCCESS", method = "TCA", canonical = canonical,
         native = native, fit = fit, failure = NULL)
  }, error = function(error) ctse_method_failure("TCA", error))
}

ctse_coerce_fraction <- function(value, sample_ids, cell_types) {
  value <- as.matrix(value)
  if (identical(rownames(value), sample_ids) && setequal(colnames(value), cell_types)) {
    value <- value[, cell_types, drop = FALSE]
  } else if (identical(colnames(value), sample_ids) &&
             setequal(rownames(value), cell_types)) {
    value <- t(value[cell_types, , drop = FALSE])
  } else {
    ctse_abort("BAYESPRISM_FRACTION_AXIS")
  }
  ctse_assert(all(is.finite(value)) && all(value >= 0) &&
                max(abs(rowSums(value) - 1)) < 1e-8,
              "BAYESPRISM_FRACTION_VALUE")
  value
}

ctse_coerce_expression <- function(value, sample_ids, cell_type) {
  value <- as.matrix(value)
  ctse_assert(!is.null(rownames(value)) && !is.null(colnames(value)),
              "BAYESPRISM_EXPRESSION_KEYS", cell_type)
  if (identical(rownames(value), sample_ids)) return(t(value))
  if (identical(colnames(value), sample_ids)) return(value)
  ctse_abort("BAYESPRISM_EXPRESSION_SAMPLE_AXIS", cell_type)
}

run_bayesprism <- function(reference_counts, cell_type_labels, cell_state_labels,
                           target_counts, seed, n_cores = 1L) {
  ctse_require_package_version("BayesPrism", "2.2.3")
  ctse_assert((is.matrix(reference_counts) || inherits(reference_counts, "Matrix")) &&
                (is.matrix(target_counts) || inherits(target_counts, "Matrix")),
              "BAYESPRISM_COUNT_MATRIX_REQUIRED")
  ctse_assert(identical(rownames(reference_counts), rownames(target_counts)),
              "BAYESPRISM_INPUT_GENE_AXIS")
  ctse_assert(length(cell_type_labels) == ncol(reference_counts) &&
                length(cell_state_labels) == ncol(reference_counts),
              "BAYESPRISM_REFERENCE_LABEL_LENGTH")
  tryCatch({
    ctse_set_rng(seed)
    prism <- BayesPrism::new.prism(
      reference = t(reference_counts), input.type = "count.matrix",
      cell.type.labels = cell_type_labels, cell.state.labels = cell_state_labels,
      key = NULL, mixture = t(target_counts), outlier.cut = 0.01,
      outlier.fraction = 0.1, pseudo.min = 1e-8
    )
    fit <- BayesPrism::run.prism(
      prism, n.cores = as.integer(n_cores), update.gibbs = TRUE,
      gibbs.control = list(), opt.control = list()
    )
    sample_ids <- colnames(target_counts)
    cell_types <- unique(as.character(cell_type_labels))
    theta0 <- ctse_coerce_fraction(
      BayesPrism::get.fraction(fit, "first", "type"), sample_ids, cell_types
    )
    thetaf <- ctse_coerce_fraction(
      BayesPrism::get.fraction(fit, "final", "type"), sample_ids, cell_types
    )
    expression <- setNames(lapply(cell_types, function(cell_type) {
      ctse_coerce_expression(
        BayesPrism::get.exp(fit, "type", cell_type), sample_ids, cell_type
      )
    }), cell_types)
    returned_genes <- Reduce(intersect, lapply(expression, rownames))
    ctse_assert(length(returned_genes) > 0L, "BAYESPRISM_EMPTY_RETURNED_UNIVERSE")
    canonical <- array(
      NA_real_, dim = c(length(returned_genes), length(cell_types), length(sample_ids)),
      dimnames = list(returned_genes, cell_types, sample_ids)
    )
    for (cell_type in cell_types) {
      canonical[, cell_type, ] <- expression[[cell_type]][returned_genes, sample_ids, drop = FALSE]
    }
    ctse_assert(all(is.finite(canonical)), "BAYESPRISM_CANONICAL_NONFINITE")
    list(
      status = "SUCCESS", method = "BayesPrism", canonical = canonical,
      theta0 = theta0, thetaf = thetaf, native_expression = expression,
      returned_genes = returned_genes,
      dropped_genes = setdiff(rownames(target_counts), returned_genes),
      fit = fit, failure = NULL
    )
  }, error = function(error) ctse_method_failure("BayesPrism", error))
}
