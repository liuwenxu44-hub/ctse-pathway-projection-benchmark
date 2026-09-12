epic_assert <- function(ok, message) { if (!isTRUE(ok)) stop(message, call. = FALSE) }
epic_hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

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

