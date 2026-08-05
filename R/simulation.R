ctse_simulation_design <- function() {
  scenarios <- data.frame(
    scenario_order = 1:6,
    scenario = c(
      "composition_only", "intrinsic_only", "concordant_mixed",
      "discordant_mixed", "complete_cancellation", "null"
    ),
    truth_class = c(
      "composition_only", "intrinsic_only", "concordant_mixed",
      "discordant_mixed", "complete_cancellation", "null"
    ),
    p0_A = 0.25, p0_B = 0.25, p0_C = 0.25, p0_D = 0.25,
    p1_A = c(0.45, 0.25, 0.45, 0.45, 0.45, 0.25),
    p1_B = c(0.20, 0.25, 0.20, 0.20, 0.20, 0.25),
    p1_C = c(0.20, 0.25, 0.20, 0.20, 0.20, 0.25),
    p1_D = c(0.15, 0.25, 0.15, 0.15, 0.15, 0.25),
    q0_A = 0.08, q0_B = 0.04, q0_C = 0.02, q0_D = 0.01,
    q1_A = c(0.08, 0.12, 0.12, 0.07, 0.0533333333333333, 0.08),
    q1_B = 0.04, q1_C = 0.02, q1_D = 0.01,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  blocks <- data.frame(
    block_order = 1:7,
    block = c(
      "signal_pathway", "invariant_negative_control_pathway",
      "marker_A", "marker_B", "marker_C", "marker_D", "background"
    ),
    gene_start = c(1L, 51L, 101L, 201L, 301L, 401L, 501L),
    gene_end = c(50L, 100L, 200L, 300L, 400L, 500L, 2000L),
    stringsAsFactors = FALSE
  )
  blocks$gene_count <- blocks$gene_end - blocks$gene_start + 1L
  list(
    schema_version = "CTSE_CONTROLLED_SIMULATION_V1",
    genes = sprintf("sg%04d", 1:2000),
    cell_types = paste0("synthetic_cell_", LETTERS[1:4]),
    states = c("state0", "state1"),
    scenarios = scenarios,
    blocks = blocks,
    reference_donors = 6L,
    reference_cells_per_state = 100L,
    bulk_samples_per_group = 20L,
    reference_library = c(median = 5000, sdlog = 0.25, lower = 1000, upper = 20000),
    bulk_library = c(median = 2000000, sdlog = 0.20, lower = 500000, upper = 8000000),
    dirichlet_concentration = 2000
  )
}

ctse_simulation_scenario <- function(design, scenario) {
  row <- design$scenarios[design$scenarios$scenario == scenario, , drop = FALSE]
  ctse_assert(nrow(row) == 1L, "SIMULATION_UNKNOWN_SCENARIO", scenario)
  row
}

ctse_simulation_gene_metadata <- function(design = ctse_simulation_design()) {
  out <- data.frame(
    gene_id = design$genes,
    gene_block = NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(design$blocks))) {
    index <- design$blocks$gene_start[[i]]:design$blocks$gene_end[[i]]
    out$gene_block[index] <- design$blocks$block[[i]]
  }
  ctse_assert(!anyNA(out$gene_block), "SIMULATION_GENE_BLOCK_MISSING")
  out
}

ctse_simulation_seed <- function(stream, scenario, replicate_id = 0L,
                                 design = ctse_simulation_design()) {
  stream_code <- c(reference = 1L, bulk = 2L, controlled_input = 3L,
                   bayesprism = 4L, bootstrap = 5L)[[stream]]
  ctse_assert(!is.null(stream_code), "SIMULATION_UNKNOWN_STREAM", stream)
  if (identical(stream, "bootstrap")) return(250000001L)
  row <- ctse_simulation_scenario(design, scenario)
  seed <- 200000000 + 10000000 * stream_code +
    100000 * as.integer(row$scenario_order[[1L]]) + as.integer(replicate_id)
  ctse_assert(seed <= .Machine$integer.max, "SIMULATION_SEED_RANGE")
  as.integer(seed)
}

ctse_truncated_lognormal <- function(n, specification, max_attempts = 100000L) {
  result <- integer(n)
  for (i in seq_len(n)) {
    accepted <- FALSE
    for (attempt in seq_len(max_attempts)) {
      candidate <- stats::rlnorm(1L, log(specification[["median"]]),
                                 specification[["sdlog"]])
      if (is.finite(candidate) && candidate >= specification[["lower"]] &&
          candidate <= specification[["upper"]]) {
        result[[i]] <- as.integer(round(candidate))
        accepted <- TRUE
        break
      }
    }
    ctse_assert(accepted, "SIMULATION_LIBRARY_SIZE_REJECTION_LIMIT")
  }
  result
}

ctse_symmetric_compositions <- function(target, n_samples = 20L) {
  target <- as.numeric(target)
  ctse_assert(length(target) == 4L && n_samples %% 2L == 0L,
              "SIMULATION_COMPOSITION_DIMENSION")
  output <- matrix(NA_real_, nrow = n_samples, ncol = 4L)
  for (pair in seq_len(n_samples / 2L)) {
    accepted <- FALSE
    for (attempt in seq_len(100000L)) {
      delta <- c(stats::rnorm(3L, 0, 0.03), 0)
      delta[[4L]] <- -sum(delta[1:3])
      plus <- target + delta
      minus <- target - delta
      if (all(plus >= 0.05 & plus <= 0.85) &&
          all(minus >= 0.05 & minus <= 0.85)) {
        output[(pair - 1L) * 2L + 1L, ] <- plus
        output[(pair - 1L) * 2L + 2L, ] <- minus
        accepted <- TRUE
        break
      }
    }
    ctse_assert(accepted, "SIMULATION_COMPOSITION_REJECTION_LIMIT")
  }
  ctse_assert(max(abs(rowSums(output) - 1)) < 1e-12 && all(output >= 0),
              "SIMULATION_COMPOSITION_INVALID")
  output
}

ctse_population_profiles <- function(scenario, design = ctse_simulation_design()) {
  row <- ctse_simulation_scenario(design, scenario)
  metadata <- ctse_simulation_gene_metadata(design)
  profiles <- array(
    0,
    dim = c(length(design$genes), length(design$cell_types), length(design$states)),
    dimnames = list(design$genes, design$cell_types, design$states)
  )
  for (state_index in seq_along(design$states)) {
    state <- design$states[[state_index]]
    for (cell_index in seq_along(design$cell_types)) {
      letter <- LETTERS[[cell_index]]
      q <- as.numeric(row[[paste0(if (state == "state0") "q0_" else "q1_", letter)]])
      masses <- c(
        signal_pathway = q,
        invariant_negative_control_pathway = 0.02,
        marker_A = if (letter == "A") 0.25 else 0.02,
        marker_B = if (letter == "B") 0.25 else 0.02,
        marker_C = if (letter == "C") 0.25 else 0.02,
        marker_D = if (letter == "D") 0.25 else 0.02,
        background = 0.67 - q
      )
      ctse_assert(abs(sum(masses) - 1) < 1e-14 && all(masses > 0),
                  "SIMULATION_PROFILE_MASS_INVALID")
      for (block in names(masses)) {
        index <- which(metadata$gene_block == block)
        profiles[index, cell_index, state_index] <- masses[[block]] / length(index)
      }
    }
  }
  ctse_assert(max(abs(apply(profiles, c(2L, 3L), sum) - 1)) < 1e-14,
              "SIMULATION_PROFILE_NORMALIZATION")
  profiles
}

ctse_simulation_pathways <- function(design = ctse_simulation_design()) {
  list(
    signal = data.frame(gene_id = design$genes[1:50], weight = 1,
                        stringsAsFactors = FALSE),
    negative_control = data.frame(gene_id = design$genes[51:100], weight = 1,
                                  stringsAsFactors = FALSE)
  )
}

ctse_latent_truth <- function(scenario, profiles,
                              design = ctse_simulation_design()) {
  row <- ctse_simulation_scenario(design, scenario)
  p0 <- as.numeric(row[paste0("p0_", LETTERS[1:4])])
  p1 <- as.numeric(row[paste0("p1_", LETTERS[1:4])])
  state0 <- profiles[, , "state0", drop = FALSE][, , 1L]
  state1 <- profiles[, , "state1", drop = FALSE][, , 1L]
  states <- list(
    X00 = as.numeric(state0 %*% p0), X10 = as.numeric(state0 %*% p1),
    X01 = as.numeric(state1 %*% p0), X11 = as.numeric(state1 %*% p1)
  )
  for (name in names(states)) names(states[[name]]) <- design$genes
  states
}

ctse_dirichlet_multinomial <- function(profile, library_size, concentration = 2000) {
  ctse_assert(all(is.finite(profile)) && all(profile > 0) &&
                abs(sum(profile) - 1) < 1e-12,
              "SIMULATION_PROFILE_INVALID")
  gamma_value <- stats::rgamma(length(profile), shape = concentration * profile, rate = 1)
  ctse_assert(all(is.finite(gamma_value)) && all(gamma_value > 0),
              "SIMULATION_DIRICHLET_INVALID")
  probability <- gamma_value / sum(gamma_value)
  as.integer(stats::rmultinom(1L, as.integer(library_size), probability)[, 1L])
}

ctse_generate_reference <- function(scenario, design = ctse_simulation_design(),
                                    profiles = ctse_population_profiles(scenario, design),
                                    donors = design$reference_donors,
                                    cells_per_state = design$reference_cells_per_state) {
  ctse_set_rng(ctse_simulation_seed("reference", scenario, 0L, design))
  donor_ids <- sprintf("sim_ref_donor_%02d", seq_len(donors))
  total_cells <- donors * length(design$cell_types) * length(design$states) * cells_per_state
  counts <- matrix(0L, nrow = length(design$genes), ncol = total_cells,
                   dimnames = list(design$genes, NULL))
  metadata <- vector("list", total_cells)
  column <- 0L
  for (donor_index in seq_along(donor_ids)) {
    for (cell_index in seq_along(design$cell_types)) {
      donor_factor <- stats::rlnorm(length(design$genes), 0, 0.10)
      for (state_index in seq_along(design$states)) {
        profile <- profiles[, cell_index, state_index] * donor_factor
        profile <- profile / sum(profile)
        for (cell_number in seq_len(cells_per_state)) {
          column <- column + 1L
          size <- ctse_truncated_lognormal(1L, design$reference_library)[[1L]]
          counts[, column] <- as.integer(stats::rmultinom(1L, size, profile)[, 1L])
          cell_id <- sprintf("ref__%s__%s__%s__%s__c%03d", scenario,
                             donor_ids[[donor_index]], design$cell_types[[cell_index]],
                             design$states[[state_index]], cell_number)
          metadata[[column]] <- data.frame(
            cell_id = cell_id, donor_id = donor_ids[[donor_index]],
            cell_type = design$cell_types[[cell_index]],
            state = design$states[[state_index]], cell_number = cell_number,
            library_size = size, stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  metadata <- do.call(rbind, metadata)
  colnames(counts) <- metadata$cell_id
  list(counts = counts, metadata = metadata, population_profiles = profiles,
       scenario = scenario, seed = ctse_simulation_seed("reference", scenario, 0L, design))
}

ctse_generate_bulk <- function(scenario, replicate_id,
                               design = ctse_simulation_design(),
                               profiles = ctse_population_profiles(scenario, design),
                               samples_per_group = design$bulk_samples_per_group) {
  ctse_set_rng(ctse_simulation_seed("bulk", scenario, replicate_id, design))
  row <- ctse_simulation_scenario(design, scenario)
  p0 <- as.numeric(row[paste0("p0_", LETTERS[1:4])])
  p1 <- as.numeric(row[paste0("p1_", LETTERS[1:4])])
  composition <- rbind(
    ctse_symmetric_compositions(p0, samples_per_group),
    ctse_symmetric_compositions(p1, samples_per_group)
  )
  colnames(composition) <- design$cell_types
  sizes <- ctse_truncated_lognormal(samples_per_group, design$bulk_library)
  n_samples <- 2L * samples_per_group
  counts <- matrix(0L, nrow = length(design$genes), ncol = n_samples,
                   dimnames = list(design$genes, NULL))
  latent <- matrix(NA_real_, nrow = length(design$genes), ncol = n_samples,
                   dimnames = list(design$genes, NULL))
  metadata <- vector("list", n_samples)
  column <- 0L
  for (group_index in 1:2) {
    group <- paste0("group", group_index - 1L)
    state <- design$states[[group_index]]
    group_composition <- composition[((group_index - 1L) * samples_per_group + 1L):
                                       (group_index * samples_per_group), , drop = FALSE]
    for (sample_number in seq_len(samples_per_group)) {
      column <- column + 1L
      latent[, column] <- as.numeric(profiles[, , state] %*% group_composition[sample_number, ])
      counts[, column] <- ctse_dirichlet_multinomial(
        latent[, column], sizes[[sample_number]], design$dirichlet_concentration
      )
      sample_id <- sprintf("%s__r%03d__%s__s%02d", scenario, replicate_id,
                           group, sample_number)
      metadata[[column]] <- data.frame(
        sample_id = sample_id, group = group, sample_number = sample_number,
        library_size = sizes[[sample_number]], stringsAsFactors = FALSE
      )
    }
  }
  metadata <- do.call(rbind, metadata)
  colnames(counts) <- colnames(latent) <- metadata$sample_id
  rownames(composition) <- metadata$sample_id
  list(counts = counts, metadata = metadata, oracle_compositions = composition,
       latent_profiles = latent, scenario = scenario, replicate_id = replicate_id,
       seed = ctse_simulation_seed("bulk", scenario, replicate_id, design))
}

ctse_truth_tensor <- function(profiles, sample_metadata) {
  genes <- dimnames(profiles)[[1L]]
  cell_types <- dimnames(profiles)[[2L]]
  samples <- as.character(sample_metadata$sample_id)
  result <- array(NA_real_, dim = c(length(genes), length(cell_types), length(samples)),
                  dimnames = list(genes, cell_types, samples))
  for (i in seq_along(samples)) {
    state <- if (sample_metadata$group[[i]] == "group0") "state0" else "state1"
    result[, , i] <- profiles[, , state]
  }
  result
}

ctse_pooled_reference_profiles <- function(reference) {
  genes <- rownames(reference$counts)
  cell_types <- unique(as.character(reference$metadata$cell_type))
  result <- matrix(NA_real_, nrow = length(genes), ncol = length(cell_types),
                   dimnames = list(genes, cell_types))
  for (cell_type in cell_types) {
    index <- reference$metadata$cell_type == cell_type
    total <- rowSums(reference$counts[, index, drop = FALSE])
    result[, cell_type] <- total / sum(total)
  }
  result
}

ctse_shared_nnls_fraction <- function(reference_profiles, bulk_counts) {
  ctse_assert(requireNamespace("nnls", quietly = TRUE), "NNLS_PACKAGE_REQUIRED")
  ctse_assert(identical(rownames(reference_profiles), rownames(bulk_counts)),
              "NNLS_GENE_AXIS_MISMATCH")
  normalized_bulk <- sweep(bulk_counts, 2L, colSums(bulk_counts), "/")
  result <- matrix(NA_real_, nrow = ncol(bulk_counts), ncol = ncol(reference_profiles),
                   dimnames = list(colnames(bulk_counts), colnames(reference_profiles)))
  for (i in seq_len(ncol(bulk_counts))) {
    value <- as.numeric(nnls::nnls(reference_profiles, normalized_bulk[, i])$x)
    ctse_assert(all(is.finite(value)) && all(value >= 0) && sum(value) > 0,
                "NNLS_NUMERICAL_FAILURE")
    result[i, ] <- value / sum(value)
  }
  result
}
