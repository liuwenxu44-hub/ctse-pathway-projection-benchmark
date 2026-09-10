# Exact frozen CTS reconstruction from saved MCMCglmm coefficients. No sampling.
epic_cts_draws <- function(sol, samples, cells) {
  sol <- as.matrix(sol)
  output <- matrix(NA_real_, nrow(sol), length(samples) * length(cells))
  names <- character(ncol(output))
  at <- 0L
  for (cell in cells) for (sample in samples) {
    at <- at + 1L
    re <- colnames(sol)[grepl(paste0("^", cell, "[.]sample_id[.]", sample, "$"), colnames(sol))]
    epic_assert(length(re) == 1L && cell %in% colnames(sol), "Ambiguous CTS chain column mapping")
    output[, at] <- sol[, cell] + sol[, re]
    names[at] <- paste(cell, sample, sep = "::")
  }
  colnames(output) <- names
  output
}
