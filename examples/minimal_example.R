source("R/ctse_benchmark.R")
source("R/adapters.R")

set.seed(20260806)
genes <- sprintf("gene_%03d", seq_len(100L))
cell_types <- c("cell_A", "cell_B")
samples <- sprintf("sample_%02d", seq_len(8L))
groups <- rep(c("control", "case"), each = 4L)

estimate <- array(
  rexp(length(genes) * length(cell_types) * length(samples), rate = 2),
  dim = c(length(genes), length(cell_types), length(samples)),
  dimnames = list(genes, cell_types, samples)
)
estimate[1:10, "cell_A", groups == "case"] <-
  estimate[1:10, "cell_A", groups == "case"] + 1

sample_metadata <- data.frame(sample_id = samples, group = groups)
pathway <- data.frame(gene_id = genes[1:10], weight = 1)

result <- run_ctse_benchmark(
  estimate = estimate,
  sample_metadata = sample_metadata,
  pathway = pathway,
  group0 = "control",
  group1 = "case"
)

stopifnot(identical(result$status, "SUCCESS"))
print(result$projection)

output_dir <- file.path(tempdir(), "ctse_benchmark_example")
paths <- write_ctse_benchmark_result(result, output_dir, "synthetic_example")
stopifnot(all(file.exists(unlist(paths))))
