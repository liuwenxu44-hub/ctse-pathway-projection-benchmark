source("R/ctse_benchmark.R")
source("R/adapters.R")

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

cat("All CTSE benchmark tests passed.\n")
