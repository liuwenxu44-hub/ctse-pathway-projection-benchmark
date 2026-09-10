# Bounded EPIC chain instrumentation. Sourcing this file never runs a fit.
# Requires the frozen EPICunmix namespace, MCMCglmm, digest and posterior.

epic_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

epic_hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

epic_validate_pinned <- function(path) {
  epic_assert(requireNamespace("EPICunmix", quietly = TRUE), "EPICunmix unavailable")
  lines <- readLines(path, warn = FALSE)
  labels <- sub("^# ", "", lines[grepl("^# [A-Za-z_][A-Za-z0-9_]*$", lines)])
  expressions <- parse(path)
  epic_assert(length(labels) == length(expressions), "Pinned function inventory mismatch")
  for (i in seq_along(labels)) {
    expected <- eval(expressions[[i]], envir = baseenv())
    actual <- get(labels[i], envir = asNamespace("EPICunmix"), inherits = FALSE)
    epic_assert(identical(formals(expected), formals(actual)), paste("Formals drift", labels[i]))
    epic_assert(identical(paste(deparse(body(expected)), collapse = "\n"),
                          paste(deparse(body(actual)), collapse = "\n")),
                paste("Function body drift", labels[i]))
  }
  invisible(labels)
}

epic_validate_input <- function(a) {
  required <- c("L", "W", "prior", "genes", "full_genes", "sample_alias", "cell_alias")
  epic_assert(all(required %in% names(a)), "Missing fit-input fields")
  epic_assert(!any(grepl("truth|oracle|evaluation|metadata|phenotype|group", names(a),
                        ignore.case = TRUE)), "Evaluation fields in fit input")
  epic_assert(identical(dimnames(a$L), list(a$genes, a$sample_alias)), "L axes mismatch")
  epic_assert(identical(dimnames(a$W), list(a$sample_alias, a$cell_alias)), "W axes mismatch")
  epic_assert(identical(dimnames(a$prior$profile), list(a$genes, a$cell_alias)), "Prior axes mismatch")
  epic_assert(identical(dimnames(a$prior$covariance),
                        list(a$genes, a$cell_alias, a$cell_alias)), "Covariance axes mismatch")
  epic_assert(all(is.finite(a$L)) && all(is.finite(a$W)) &&
                all(is.finite(a$prior$profile)) && all(is.finite(a$prior$covariance)),
              "Nonfinite fit input")
  epic_assert(!anyDuplicated(a$genes) && !anyDuplicated(a$sample_alias) &&
                !anyDuplicated(a$cell_alias), "Duplicate fit axes")
  epic_assert(all(a$genes %in% a$full_genes), "Unrecognized gene IDs")
  epic_assert(all(grepl("^c[0-9]+$", a$cell_alias)) &&
                all(grepl("^s[0-9]+$", a$sample_alias)), "Unsafe official parser aliases")
  epic_assert(all(a$W >= 0) && max(abs(rowSums(a$W) - 1)) < 1e-8,
              "Invalid frozen fractions")
  invisible(TRUE)
}

# Obtain the actual official stage-2 arguments without sampling: only this
# private clone's bmind_de binding is intercepted. The package is never edited.
epic_prepare_target <- function(a, stage, first = NULL) {
  epic_validate_input(a)
  epic_assert(stage %in% c(1L, 2L), "Stage must be 1 or 2")
  if (stage == 1L) {
    X <- a$L
    profile <- a$prior$profile
    covariance <- a$prior$covariance
    nu <- 50
    upstream <- NULL
  } else {
    epic_assert(is.list(first) && length(dim(first$A)) == 3L, "Frozen stage-one A required")
    g <- rownames(first$A)
    epic_assert(identical(g, a$genes[a$genes %in% g]) && !anyDuplicated(g),
                "Frozen stage-one gene order mismatch")
    epic_assert(identical(dimnames(first$A)[[2]], a$cell_alias) &&
                  identical(dimnames(first$A)[[3]], a$sample_alias) &&
                  all(is.finite(first$A)), "Frozen stage-one A axes/values invalid")
    official <- get("run_epic_unmix", asNamespace("EPICunmix"))
    capture <- new.env(parent = emptyenv())
    shim <- new.env(parent = environment(official))
    shim$bmind_de <- function(bulk, frac, profile, covariance, noRE, nu, seed, ncore) {
      capture$args <- list(X = bulk, W = frac, profile = profile,
                           covariance = covariance, nu = nu, noRE = noRE)
      list(A = first$A[rownames(bulk), , , drop = FALSE])
    }
    environment(official) <- shim
    suppressMessages(official(bulk = a$L[g, , drop = FALSE], frac = a$W,
                              input_cts = first, outf = FALSE, nstop = 1,
                              delta = .1, nu0 = 50, nu1 = 50, seed = 123, ncore = 1))
    z <- capture$args
    epic_assert(is.list(z) && identical(z$W, a$W) && identical(z$noRE, FALSE),
                "Official stage-two argument capture failed")
    X <- z$X
    profile <- z$profile
    covariance <- z$covariance
    nu <- z$nu
    upstream <- epic_hash(first$A)
    epic_assert(nu == 100, "Unexpected second-stage nu")
  }
  epic_assert(identical(rownames(X), rownames(profile)) &&
                identical(rownames(profile), rownames(covariance)), "Captured prior alignment failure")
  # Official stage-two covariance array lacks its third dimnames in this build.
  # Do not repair/mutate that object; exact captured values enter bmind1.
  target <- list(schema = "EPIC_FIXED_POSTERIOR_TARGET_V1", stage = as.integer(stage),
                 X = X, W = a$W, profile = profile, covariance = covariance,
                 nu = nu, bounds = range(X), full_genes = a$full_genes,
                 denominator = a$full_denominator, sample_alias = a$sample_alias,
                 cell_alias = a$cell_alias, upstream_A_sha256 = upstream,
                 scale = a$scale, source_input_sha256 = a$source_input_sha256)
  target$target_sha256 <- epic_hash(target)
  target
}

epic_check_target <- function(target) {
  stored <- target$target_sha256
  copy <- target
  copy$target_sha256 <- NULL
  epic_assert(identical(stored, epic_hash(copy)), "Posterior target mutated")
  epic_assert(identical(target$bounds, range(target$X)), "Full-support clamp bounds changed")
  invisible(TRUE)
}

# Clone the exact original leaf function; intercept only MCMCglmm to retain
# its returned chain and the effective RNG state, without changing any call.
epic_instrument_leaf <- function(stage) {
  name <- if (stage == 1L) "lme_mc2" else "bmind1"
  fn <- get(name, asNamespace("EPICunmix"))
  capture <- new.env(parent = emptyenv())
  capture$n <- 0L
  env <- new.env(parent = environment(fn))
  env$MCMCglmm <- function(...) {
    capture$n <- capture$n + 1L
    capture$rng_before_sampler <- get(".Random.seed", envir = .GlobalEnv)
    capture$rng_kind <- RNGkind()
    callback <- getOption("comppath.on_sampler_start")
    if (is.function(callback)) callback()
    capture$fit <- MCMCglmm::MCMCglmm(...)
    capture$fit
  }
  environment(fn) <- env
  list(fn = fn, capture = capture, original_name = name)
}

epic_leaf_arguments <- function(target, gene, seed, nitt, burnin, thin) {
  epic_check_target(target)
  epic_assert(length(gene) == 1L && gene %in% rownames(target$X), "Gene unavailable in target")
  epic_assert(length(seed) == 1L && is.finite(seed) && seed == as.integer(seed) && seed > 0,
              "Invalid effective seed")
  epic_assert(nitt > burnin && burnin >= 0 && thin >= 1 &&
                (nitt - burnin) %% thin == 0, "Invalid MCMC schedule")
  common <- list(x = target$X[gene, ], W = target$W,
                 sample_id = rownames(target$W), mu = target$profile[gene, ],
                 nu = target$nu, nitt = nitt, burnin = burnin, thin = thin,
                 seed = as.integer(seed))
  if (target$stage == 1L) {
    common$V_fe <- diag(.5, ncol(target$W))
    common$V_re <- target$covariance[gene, , ]
  } else {
    common$var_fe <- target$covariance[gene, , ]
    common$noRE <- FALSE
    common$np <- FALSE
  }
  common
}

epic_cts_draws <- function(sol, samples, cells) {
  sol <- as.matrix(sol)
  output <- matrix(NA_real_, nrow(sol), length(samples) * length(cells))
  names <- character(ncol(output))
  at <- 0L
  for (cell in cells) for (sample in samples) {
    at <- at + 1L
    # Match exact MCMCglmm random-effect column by its actual cell/sample ID.
    re <- colnames(sol)[grepl(paste0("^", cell, "[.]sample_id[.]", sample, "$"), colnames(sol))]
    epic_assert(length(re) == 1L && cell %in% colnames(sol), "Ambiguous CTS chain column mapping")
    output[, at] <- sol[, cell] + sol[, re]
    names[at] <- paste(cell, sample, sep = "::")
  }
  colnames(output) <- names
  output
}

epic_fit_gene_chain <- function(target, gene, seed, nitt = 1300L,
                                burnin = 300L, thin = 1L) {
  args <- epic_leaf_arguments(target, gene, seed, nitt, burnin, thin)
  instrument <- epic_instrument_leaf(target$stage)
  warnings <- character()
  result <- withCallingHandlers(do.call(instrument$fn, args), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  cap <- instrument$capture
  epic_assert(cap$n == 1L, "Expected exactly one sampler invocation")
  sol <- as.matrix(cap$fit$Sol)
  vcv <- as.matrix(cap$fit$VCV)
  epic_assert(nrow(sol) == (nitt - burnin) / thin && nrow(vcv) == nrow(sol),
              "Unexpected saved-chain length")
  epic_assert(all(is.finite(sol)) && all(is.finite(vcv)), "Nonfinite sampler draws")
  cts <- epic_cts_draws(sol, target$sample_alias, target$cell_alias)
  clamped <- pmin(pmax(result$A, target$bounds[1]), target$bounds[2])
  epic_check_target(target)
  list(schema = "EPIC_GENE_CHAIN_V1", target_sha256 = target$target_sha256,
       stage = target$stage, gene = gene, effective_seed = as.integer(seed),
       schedule = c(nitt = nitt, burnin = burnin, thin = thin),
       rng_kind = cap$rng_kind, rng_before_sampler = cap$rng_before_sampler,
       original_function = instrument$original_name,
       original_function_body_sha256 = epic_hash(body(get(instrument$original_name,
                                                         asNamespace("EPICunmix")))),
       Sol = sol, VCV = vcv, CTS = cts,
       original_leaf_result = result, A_official_clamped = clamped,
       full_support_bounds = target$bounds, warnings = warnings,
       packages = sapply(c("EPICunmix", "MCMCglmm", "digest"),
                         function(p) as.character(utils::packageVersion(p))))
}

epic_diagnostic_matrix <- function(x, rhat_max = 1.01, ess_min = 400,
                                    relative_mcse_max = .05) {
  epic_assert(requireNamespace("posterior", quietly = TRUE),
              "posterior package required; do not substitute approximate diagnostics")
  epic_assert(is.matrix(x) && ncol(x) >= 4L && nrow(x) >= 20L,
              "At least four equal-length chains required")
  bad <- !all(is.finite(x))
  constant <- !bad && any(apply(x, 2, function(y) length(unique(y)) < 2L))
  values <- c(rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_,
              mcse_mean = NA_real_, posterior_sd = NA_real_, relative_mcse = NA_real_)
  if (!bad && !constant) {
    values["rhat"] <- posterior::rhat(x)
    values["ess_bulk"] <- posterior::ess_bulk(x)
    values["ess_tail"] <- posterior::ess_tail(x)
    values["mcse_mean"] <- posterior::mcse_mean(x)
    values["posterior_sd"] <- stats::sd(as.vector(x))
    values["relative_mcse"] <- values["mcse_mean"] / values["posterior_sd"]
  }
  pass <- !bad && !constant && all(is.finite(values)) &&
    values["rhat"] <= rhat_max && values["ess_bulk"] >= ess_min &&
    values["ess_tail"] >= ess_min && values["relative_mcse"] <= relative_mcse_max
  c(as.list(values), list(status = if (bad) "NONFINITE" else if (constant) "CONSTANT_CHAIN" else
                           if (pass) "PASS" else "DIAGNOSTICS_NOT_PASS"))
}

epic_diagnose_chains <- function(records) {
  epic_assert(length(records) == 4L, "Exactly four prespecified chains required")
  keys <- c("target_sha256", "stage", "gene", "schedule", "full_support_bounds")
  for (key in keys) epic_assert(all(vapply(records, function(z)
    identical(z[[key]], records[[1]][[key]]), logical(1))), paste("Different chain target:", key))
  epic_assert(length(unique(vapply(records, `[[`, integer(1), "effective_seed"))) == 4L,
              "Effective seeds are not distinct")
  states <- vapply(records, function(z) epic_hash(z$rng_before_sampler), character(1))
  epic_assert(length(unique(states)) == 4L, "Sampler RNG states are not distinct")
  answer <- list()
  for (quantity in c("Sol", "VCV", "CTS")) {
    first <- records[[1]][[quantity]]
    epic_assert(all(vapply(records, function(z) identical(dim(z[[quantity]]), dim(first)) &&
      identical(colnames(z[[quantity]]), colnames(first)), logical(1))), "Chain axes mismatch")
    for (j in seq_len(ncol(first))) {
      x <- do.call(cbind, lapply(records, function(z) z[[quantity]][, j]))
      row <- c(list(quantity = quantity, parameter = colnames(first)[j]), epic_diagnostic_matrix(x))
      answer[[length(answer) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, answer)
  if (any(lengths(lapply(records, `[[`, "warnings")) > 0L))
    out$status[out$status == "PASS"] <- "WARNING_REQUIRES_REVIEW"
  out
}

epic_run_gene <- function(target, gene, seeds, nitt = 1300L, burnin = 300L) {
  epic_assert(length(seeds) == 4L && !anyDuplicated(seeds), "Four unique effective seeds required")
  epic_assert(requireNamespace("posterior", quietly = TRUE), "posterior unavailable before fitting")
  records <- lapply(seeds, function(seed) epic_fit_gene_chain(target, gene, seed, nitt, burnin))
  list(target_sha256 = target$target_sha256, chains = records,
       diagnostics = epic_diagnose_chains(records),
       diagnostic_package = as.character(utils::packageVersion("posterior")))
}
