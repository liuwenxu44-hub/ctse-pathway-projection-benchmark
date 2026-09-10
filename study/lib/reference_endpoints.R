# Deterministic evaluation only. No package loaders, NNLS or model entrypoints.
# Sourceable library for artificial tests; real CLI requires both terminal ledgers.
re_assert <- function(ok, code, detail = "") {
  if (!isTRUE(ok)) stop(paste(code, detail), call. = FALSE)
}
re_read <- function(p) read.delim(p, sep = "\t", quote = "", comment.char = "",
                                check.names = FALSE, stringsAsFactors = FALSE,
                                colClasses = "character", na.strings = "NA")
re_write <- function(x, p, append = FALSE) write.table(x, p, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = !append, append = append, na = "NA")
re_sha <- function(p) unname(tools::sha256sum(p))
re_methods <- c("bayesprism", "unico", "tca", "epic", "reference", "bulk")
re_regimes <- c("ORIGINAL_POOLED_4800", "STATE0_ONLY_2400")
re_load_math <- function(path, expected) {
  re_assert(identical(re_sha(path), expected), "ENDPOINT_MATH_HASH_DRIFT")
  env <- new.env(parent = baseenv())
  # Only the original pure endpoint functions are loaded, not build_endpoints.R.
  # stats is explicitly inherited because safe_rank uses sd()/cor().
  parent.env(env) <- asNamespace("stats")
  sys.source(path, envir = env)
  re_assert(all(c("safe_rank", "projection", "classify", "strata_masks", "class_summary") %in%
                  ls(env)), "ORIGINAL_ENDPOINT_FUNCTION_MISSING")
  env
}
re_validate_tensor <- function(x, genes, cells, samples, name, allow_missing = TRUE) {
  if (is.null(x)) { re_assert(allow_missing, "REQUIRED_TENSOR_ABSENT", name); return(invisible(NULL)) }
  re_assert(is.array(x) && length(dim(x)) == 3L && all(is.finite(x)) &&
              length(rownames(x)) > 0L && !anyDuplicated(rownames(x)) &&
              all(rownames(x) %in% genes) &&
              identical(dimnames(x)[[2]], cells) && identical(dimnames(x)[[3]], samples),
            "TENSOR_AXIS_OR_VALUE_MISMATCH", name)
  invisible(NULL)
}
re_support <- function(regimes, genes) {
  re_assert(identical(names(regimes), re_regimes), "REFERENCE_REGIME_ROSTER_MISMATCH")
  supports <- lapply(regimes, function(z) {
    re_assert(identical(names(z$mats), re_methods), "METHOD_ROSTER_MISMATCH")
    available <- lapply(z$mats[!vapply(z$mats, is.null, logical(1))], rownames)
    re_assert(length(available) >= 2L, "BASELINES_MISSING")
    genes[genes %in% Reduce(intersect, available)]
  })
  paired <- genes[genes %in% Reduce(intersect, supports)]
  # Empty paired support blocks; it never triggers a union or support relaxation.
  re_assert(length(paired) > 0L, "NO_ACTUAL_PAIRED_SUPPORT")
  list(regime = supports, paired = paired)
}
re_unit_endpoints <- function(unit, genes, cells, samples, ev, regimes, math) {
  re_assert(length(genes) > 0L && !anyDuplicated(genes) && !anyDuplicated(cells) &&
              !anyDuplicated(samples), "INPUT_ROSTER_INVALID")
  re_validate_tensor(ev$truth_linear, genes, cells, samples, "truth_linear", FALSE)
  re_validate_tensor(ev$truth_native_log, genes, cells, samples, "truth_native_log", FALSE)
  re_assert(identical(rownames(ev$truth_linear), genes) &&
              identical(dimnames(ev$design_changed), list(genes, cells)) &&
              is.logical(ev$design_changed) && !anyNA(ev$design_changed) &&
              identical(as.character(ev$metadata$sample_id), samples), "EVALUATION_SCHEMA_INVALID")
  grp <- as.character(ev$metadata$group)
  re_assert(setequal(grp, c("group0", "group1")) &&
              length(unique(as.integer(table(grp)))) == 1L, "UNEQUAL_OR_MISSING_ORIGINAL_GROUPS")
  for (reg in names(regimes)) {
    z <- regimes[[reg]]
    for (m in re_methods) re_validate_tensor(z$mats[[m]], genes, cells, samples, paste(reg, m))
    for (m in c("reference", "bulk")) {
      re_assert(identical(rownames(z$mats[[m]]), genes), "BASELINE_FULL_SUPPORT_REQUIRED")
      re_validate_tensor(z$log_baselines[[m]], genes, cells, samples, paste(reg, "log", m), FALSE)
    }
  }
  re_assert(identical(regimes[[1]]$mats$bulk, regimes[[2]]$mats$bulk) &&
              identical(regimes[[1]]$log_baselines$bulk, regimes[[2]]$log_baselines$bulk),
            "UNCHANGED_BULK_BASELINE_DIFFERED")
  supports <- re_support(regimes, genes); common <- supports$paired
  complete <- vapply(regimes, function(z) all(!vapply(z$mats[re_methods[1:4]], is.null, logical(1))), logical(1))
  both_complete <- all(complete)
  base <- data.frame(task_id = as.character(unit$task_id), scenario = as.character(unit$scenario),
                     replicate_id = as.integer(unit$replicate_id), complete_four_both_regimes = both_complete)
  rows <- list(); emit <- function(name, x) {
    if (is.null(rows[[name]])) rows[[name]] <<- list()
    rows[[name]][[length(rows[[name]]) + 1L]] <<- cbind(base[rep(1L, nrow(x)), , drop = FALSE], x,
                                                       row.names = NULL)
  }
  emit("PAIRED_SUPPORT_LONG", data.frame(gene_id = genes,
    original_regime_common = genes %in% supports$regime[[1]],
    state0_regime_common = genes %in% supports$regime[[2]], paired_presence = genes %in% common))
  for (reg in re_regimes) {
    z <- regimes[[reg]]
    for (m in re_methods) {
      x <- z$mats[[m]]
      emit("METHOD_AVAILABILITY", data.frame(reference_regime = reg, method = m,
        original_input_genes = length(genes), returned_genes = if (is.null(x)) 0L else nrow(x),
        own_regime_common_genes = length(supports$regime[[reg]]), paired_common_genes = length(common),
        available = !is.null(x), failure = if (is.null(x)) z$failure[[m]] else "",
        complete_four_this_regime = complete[[reg]], celltypes = length(cells), samples = length(samples),
        denominator_status = "NEW_PAIRED_SUPPORT_SENSITIVITY_NOT_REPLACEMENT_OF_OLD_PRIMARY"))
      emit("RETURNED_SUPPORT_LONG", data.frame(reference_regime = reg, method = m,
        gene_id = genes, returned_presence = genes %in% rownames(x), paired_presence = genes %in% common))
    }
    for (k in seq_along(cells)) {
      masks <- math$strata_masks(ev, common, k)
      typechanged <- any(ev$design_changed[, k])
      emit("DESIGN_CHANGE_STRATA", data.frame(reference_regime = reg, celltype = cells[k],
        celltype_design_changed = typechanged, all_input_changed_genes = sum(ev$design_changed[, k]),
        own_regime_changed_genes = sum(ev$design_changed[supports$regime[[reg]], k]),
        paired_changed_genes = sum(masks$DESIGN_CHANGED), paired_unchanged_genes = sum(masks$DESIGN_UNCHANGED)))
      Yrank <- matrix(ev$truth_linear[common, k, , drop = FALSE], nrow = length(common))
      for (m in re_methods) {
        X <- z$mats[[m]]
        for (s in seq_along(samples)) for (st in names(masks)) {
          keep <- masks[[st]]; n <- sum(keep)
          xx <- if (is.null(X)) NULL else X[common, k, s]
          rho <- if (is.null(xx) || n < 3L) NA_real_ else math$safe_rank(xx[keep], Yrank[keep, s])
          reason <- if (is.null(xx)) "METHOD_UNAVAILABLE" else if (n < 3L)
            "INSUFFICIENT_DESIGN_STRATUM_GENES" else if (!is.finite(rho))
              "CONSTANT_VECTOR_RANK_UNDEFINED" else "EVALUABLE"
          emit("GENE_RANK_RECOVERY_LONG", data.frame(reference_regime = reg, method = m,
            celltype = cells[k], celltype_design_changed = typechanged, sample_id = samples[s], group = grp[s],
            gene_stratum = st, genes = n, spearman = rho, status = reason))
        }
      }
      contexts <- list(LINEAR_CPM_SECONDARY = list(unico = z$mats$unico, tca = z$mats$tca,
        reference = z$mats$reference, bulk = z$mats$bulk), EPIC_NATIVE_LOG_SECONDARY = list(
          epic = z$mats$epic, reference = z$log_baselines$reference, bulk = z$log_baselines$bulk))
      for (ctx in names(contexts)) {
        truth <- if (ctx == "LINEAR_CPM_SECONDARY") ev$truth_linear else ev$truth_native_log
        Y <- matrix(truth[common, k, , drop = FALSE], nrow = length(common))
        td <- rowMeans(Y[, grp == "group1", drop = FALSE]) - rowMeans(Y[, grp == "group0", drop = FALSE])
        for (m in names(contexts[[ctx]])) {
          X <- contexts[[ctx]][[m]]
          xx <- if (is.null(X)) NULL else matrix(X[common, k, , drop = FALSE], nrow = length(common))
          delta <- if (is.null(xx)) NULL else rowMeans(xx[, grp == "group1", drop = FALSE]) -
            rowMeans(xx[, grp == "group0", drop = FALSE])
          for (s in seq_along(samples)) for (st in names(masks)) {
            keep <- masks[[st]]; n <- sum(keep)
            err <- if (is.null(xx) || !n) NA_real_ else sqrt(mean((xx[keep, s] - Y[keep, s])^2))
            den <- if (n >= 2L) max(sd(Y[keep, s]), 1e-8) else NA_real_
            emit("SAME_SCALE_ERROR_LONG", data.frame(reference_regime = reg, context = ctx, method = m,
              celltype = cells[k], celltype_design_changed = typechanged, sample_id = samples[s], group = grp[s],
              gene_stratum = st, genes = n, rmse = err, sne = err / den,
              status = if (is.null(xx)) "METHOD_UNAVAILABLE" else if (!n) "EMPTY_DESIGN_STRATUM" else
                "SAME_UNIT_LATENT_TARGET_DIAGNOSTIC"))
          }
          for (st in names(masks)) {
            keep <- masks[[st]]; n <- sum(keep)
            emit("SCALE_SPECIFIC_CONTRAST_LONG", data.frame(reference_regime = reg, context = ctx, method = m,
              celltype = cells[k], celltype_design_changed = typechanged, gene_stratum = st, genes = n,
              truth_contrast_L2 = sqrt(sum(td[keep]^2)),
              estimated_contrast_L2 = if (is.null(delta) || !n) NA_real_ else sqrt(sum(delta[keep]^2)),
              contrast_error_L2 = if (is.null(delta) || !n) NA_real_ else sqrt(sum((delta[keep] - td[keep])^2))))
          }
          for (path in c("signal", "negative_control")) {
            tp <- math$projection(td, common, path)
            p <- if (is.null(delta)) c(P = NA_real_, norm = NA_real_, members = tp["members"][[1]]) else
              math$projection(delta, common, path)
            emit("SCALE_SPECIFIC_DIRECTION_LONG", data.frame(reference_regime = reg, context = ctx, method = m,
              celltype = cells[k], celltype_design_changed = typechanged, pathway = path,
              genes = length(common), pathway_members = unname(p["members"]),
              estimated_P = unname(p["P"]), truth_P = unname(tp["P"]),
              estimated_direction = unname(math$classify(p["P"])), truth_direction = unname(math$classify(tp["P"])),
              evaluable = is.finite(p["P"]), semantic_role = "SCALE_SPECIFIC_SECONDARY_NOT_FOUR_METHOD_RANKING"))
          }
        }
        for (path in c("signal", "negative_control")) {
          tp <- math$projection(td, common, path)
          emit("SCALE_SPECIFIC_DIRECTION_LONG", data.frame(reference_regime = reg, context = ctx, method = "always_zero",
            celltype = cells[k], celltype_design_changed = typechanged, pathway = path, genes = length(common),
            pathway_members = unname(tp["members"]), estimated_P = 0, truth_P = unname(tp["P"]),
            estimated_direction = "zero", truth_direction = unname(math$classify(tp["P"])), evaluable = TRUE,
            semantic_role = "ANALYTIC_ZERO_PREDICTION_NO_MODEL"))
        }
      }
    }
  }
  lapply(rows, function(x) do.call(rbind, x))
}

