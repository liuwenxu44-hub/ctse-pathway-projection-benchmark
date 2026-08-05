source("R/ctse_benchmark.R")
source("R/adapters.R")
source("R/simulation.R")
source("R/method_wrappers.R")
source("R/study_endpoints.R")

genes <- paste0("g", seq_len(6L))
cell_types <- c("A", "B")
samples <- paste0("s", seq_len(4L))
metadata <- data.frame(sample_id = samples, group = c("control", "control", "case", "case"))
pathway <- data.frame(gene_id = c("g1", "g2"), weight = c(1, 1))

x <- array(0, dim = c(6L, 2L, 4L), dimnames = list(genes, cell_types, samples))
x[1:2, "A", 3:4] <- 1
x[3, "A", 3:4] <- 1

canonical <- canonicalize_ctse_array(x, genes, cell_types, samples)
stopifnot(identical(canonical, x))

source_gene_sample <- aperm(x, c(2L, 1L, 3L))
stopifnot(identical(canonicalize_source_gene_sample(source_gene_sample, genes, cell_types, samples), x))

matrix_list <- setNames(lapply(seq_along(cell_types), function(k) x[, k, ]), cell_types)
stopifnot(identical(canonicalize_celltype_matrices(matrix_list, genes, cell_types, samples), x))

projection <- signed_pathway_projection(x, metadata, pathway, "control", "case")
expected <- 2 / (sqrt(2) * sqrt(3))
stopifnot(abs(projection$signed_projection[projection$cell_type == "A"] - expected) < 1e-12)
stopifnot(projection$signed_projection[projection$cell_type == "B"] == 0)

scaled <- signed_pathway_projection(x * 7, metadata, pathway, "control", "case")
stopifnot(isTRUE(all.equal(projection$signed_projection, scaled$signed_projection, tolerance = 1e-12)))

directions <- classify_projection_direction(c(-0.051, -0.05, 0, 0.05, 0.051, NA_real_))
stopifnot(identical(directions, c("negative", "zero", "zero", "zero", "positive", NA_character_)))

metrics <- ctse_vector_metrics(c(1, 2, 3), c(1, 2, 4))
stopifnot(metrics$n == 3L, is.finite(metrics$sne), metrics$sne >= 0)

success <- run_ctse_benchmark(x, metadata, pathway, "control", "case")
stopifnot(identical(success$status, "SUCCESS"), nrow(success$projection) == 2L)
stopifnot(is.null(success$comparison))

with_truth <- run_ctse_benchmark(x, metadata, pathway, "control", "case", truth = x)
stopifnot(identical(with_truth$status, "SUCCESS"))
stopifnot(all(with_truth$comparison$absolute_error == 0))
stopifnot(all(with_truth$comparison$direction_correct))

bad <- x
dimnames(bad)[[1L]][2L] <- dimnames(bad)[[1L]][1L]
failure <- run_ctse_benchmark(bad, metadata, pathway, "control", "case")
stopifnot(identical(failure$status, "STRUCTURED_FAILURE"))
stopifnot(identical(failure$failure$failure_code, "IDENTIFIER_DUPLICATED"))

output_dir <- file.path(tempdir(), "ctse_benchmark_tests")
unlink(output_dir, recursive = TRUE)
written <- write_ctse_benchmark_result(success, output_dir, "task-001")
stopifnot(all(file.exists(unlist(written))))
stopifnot(length(list.files(output_dir, pattern = "SUCCESS$")) == 1L)
stopifnot(length(list.files(output_dir, pattern = "STRUCTURED_FAILURE$")) == 0L)

failure_written <- write_ctse_benchmark_result(failure, output_dir, "task-002")
stopifnot(all(file.exists(unlist(failure_written))))
stopifnot(length(list.files(output_dir, pattern = "STRUCTURED_FAILURE$")) == 1L)

cat("Core CTSE benchmark tests passed.\n")

design <- ctse_simulation_design()
stopifnot(nrow(design$scenarios) == 6L)
stopifnot(length(design$genes) == 2000L)
profiles <- ctse_population_profiles("complete_cancellation", design)
stopifnot(identical(dim(profiles), c(2000L, 4L, 2L)))
stopifnot(max(abs(apply(profiles, c(2L, 3L), sum) - 1)) < 1e-14)

bulk <- ctse_generate_bulk(
  "complete_cancellation", 1L, design, profiles, samples_per_group = 4L
)
stopifnot(identical(dim(bulk$counts), c(2000L, 8L)))
stopifnot(all(colSums(bulk$counts) == bulk$metadata$library_size))
stopifnot(max(abs(rowSums(bulk$oracle_compositions) - 1)) < 1e-12)
truth_tensor <- ctse_truth_tensor(profiles, bulk$metadata)
stopifnot(identical(dim(truth_tensor), c(2000L, 4L, 8L)))

simulation_result <- run_ctse_benchmark(
  truth_tensor, bulk$metadata, ctse_simulation_pathways(design)$signal,
  "group0", "group1", truth = truth_tensor
)
stopifnot(simulation_result$status == "SUCCESS")
stopifnot(all(simulation_result$comparison$absolute_error == 0))

fraction_accuracy <- ctse_fraction_accuracy(
  bulk$oracle_compositions, bulk$oracle_compositions
)
stopifnot(all(fraction_accuracy$absolute_error == 0))
gene_accuracy <- ctse_gene_accuracy(truth_tensor, truth_tensor)
stopifnot(all(gene_accuracy$rmse == 0))
interval <- ctse_bootstrap_mean_ci(1:10, seed = 20260729L, resamples = 100L)
stopifnot(interval[["estimate"]] == 5.5)

sample_map <- data.frame(
  sample_id = dimnames(truth_tensor)[[3L]],
  mixture_id = rep(c("Mix1", "Mix2"), each = 4L),
  replicate_id = rep(1:4, times = 2L),
  stringsAsFactors = FALSE
)
c1_proxy <- profiles[, , "state0"]
c2_proxy <- array(
  NA_real_, dim = c(2000L, 4L, 2L),
  dimnames = list(design$genes, design$cell_types, c("Mix1", "Mix2"))
)
c2_proxy[, , 1L] <- profiles[, , "state0"]
c2_proxy[, , 2L] <- profiles[, , "state1"]
c1_endpoint <- ctse_external_endpoint_levels(truth_tensor, c1_proxy, sample_map, "C1")
c2_endpoint <- ctse_external_endpoint_levels(truth_tensor, c2_proxy, sample_map, "C2")
stopifnot(sum(c1_endpoint$aggregation_level == "LEVEL3_MIXTURE_MEDIAN") == 2L)
stopifnot(sum(c2_endpoint$aggregation_level == "LEVEL3_MIXTURE_MEDIAN") == 2L)
stopifnot(all(c1_endpoint$status == "SUCCESS"), all(c2_endpoint$status == "SUCCESS"))

mock_fraction <- matrix(
  c(0.7, 0.3, 0.2, 0.8), nrow = 2L, byrow = TRUE,
  dimnames = list(c("sample_1", "sample_2"), c("A", "B"))
)
stopifnot(identical(
  ctse_coerce_fraction(mock_fraction, c("sample_1", "sample_2"), c("A", "B")),
  mock_fraction
))
mock_expression <- matrix(
  1:6, nrow = 2L,
  dimnames = list(c("sample_1", "sample_2"), c("gene_1", "gene_2", "gene_3"))
)
stopifnot(identical(
  ctse_coerce_expression(mock_expression, c("sample_1", "sample_2"), "A"),
  t(mock_expression)
))

cat("All study companion tests passed.\n")
