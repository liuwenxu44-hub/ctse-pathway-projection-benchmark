# Artificial tensors only. Does not read any scientific RDS or model package.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
source(args[1])
math <- re_load_math(args[2], re_sha(args[2]))
out <- args[3]
stopifnot(!file.exists(out)); dir.create(out)
tests <- list()
check <- function(name, ok) {
  stopifnot(isTRUE(ok)); tests[[length(tests) + 1L]] <<- data.frame(test = name, status = "PASS")
}
fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
genes <- sprintf("sg%04d", c(1:6, 51:56, 101:106, 501:506))
cells <- c("synthetic_A", "synthetic_B"); samples <- paste0("s", 1:4)
group <- c("group0", "group0", "group1", "group1")
y <- array(seq_len(24L * 2L * 4L) / 10, c(24L, 2L, 4L), dimnames = list(genes, cells, samples))
y[, , 3:4] <- y[, , 1:2]
y[1:6, 1, 3:4] <- y[1:6, 1, 3:4] + .8
changed <- matrix(FALSE, length(genes), length(cells), dimnames = list(genes, cells)); changed[1:6, 1] <- TRUE
ev <- list(truth_linear = y, truth_native_log = log2(y + 1),
           metadata = data.frame(sample_id = samples, group = group), design_changed = changed)
ref <- y
for (s in 2:4) ref[, , s] <- ref[, , 1]
bulk <- y * 1.1
make_regime <- function() list(mats = list(bayesprism = y / 100, unico = y * .9, tca = y * 1.2,
  epic = log2(y + 1), reference = ref, bulk = bulk),
  failure = setNames(as.list(rep("", 6)), re_methods),
  log_baselines = list(reference = log2(ref + 1), bulk = log2(bulk + 1)))
old <- make_regime(); new <- make_regime()
old$mats["tca"] <- list(NULL); old$failure$tca <- "PRESERVED_ORIGINAL_TCA_FAILURE"
old$mats$unico <- old$mats$unico[genes != "sg0501", , , drop = FALSE]
old$mats$epic <- old$mats$epic[genes != "sg0051", , , drop = FALSE]
new$mats$epic <- new$mats$epic[genes != "sg0001", , , drop = FALSE]
new$mats$reference <- ref * .8
new$log_baselines$reference <- log2(new$mats$reference + 1)
regimes <- setNames(list(old, new), re_regimes)
support <- re_support(regimes, genes)
check("intersection_across_both_regimes_all_available_methods",
      identical(support$paired, setdiff(genes, c("sg0501", "sg0051", "sg0001"))))
check("old_support_recorded_separately_from_new_paired_denominator",
      length(support$regime[[1]]) == 22L && length(support$regime[[2]]) == 23L && length(support$paired) == 21L)
unit <- data.frame(task_id = "ARTIFICIAL__r001", scenario = "ARTIFICIAL", replicate_id = 1L)
r <- re_unit_endpoints(unit, genes, cells, samples, ev, regimes, math)
av <- r$METHOD_AVAILABILITY
check("all_six_method_slots_in_both_regimes", nrow(av) == 12L && all(table(av$reference_regime) == 6L))
check("original_tca_failure_not_imputed", sum(!av$available) == 1L &&
      av$failure[!av$available] == "PRESERVED_ORIGINAL_TCA_FAILURE" && !any(av$complete_four_both_regimes))
rank <- r$GENE_RANK_RECOVERY_LONG
check("rank_all_methods_both_groups_and_design_strata_retained", nrow(rank) == 2 * 6 * 2 * 4 * 6 &&
      setequal(rank$gene_stratum, names(math$strata_masks(ev, support$paired, 1))) && setequal(rank$group, group))
check("rank_missing_method_stays_NA", all(is.na(rank$spearman[rank$reference_regime == re_regimes[1] & rank$method == "tca"])))
check("unchanged_celltype_empty_changed_stratum_retained", all(rank$genes[rank$celltype == cells[2] &
      rank$gene_stratum == "DESIGN_CHANGED"] == 0L))
g <- support$paired
row <- rank[rank$reference_regime == re_regimes[1] & rank$method == "unico" & rank$celltype == cells[1] &
      rank$sample_id == samples[1] & rank$gene_stratum == "ALL", ]
check("exact_original_gene_vector_rank", identical(row$spearman, math$safe_rank(old$mats$unico[g, 1, 1], y[g, 1, 1])))
er <- r$SAME_SCALE_ERROR_LONG
check("legal_scale_contexts_only", !any(er$method == "bayesprism") &&
      setequal(er$method[er$context == "LINEAR_CPM_SECONDARY"], c("unico", "tca", "reference", "bulk")) &&
      setequal(er$method[er$context == "EPIC_NATIVE_LOG_SECONDARY"], c("epic", "reference", "bulk")))
row <- er[er$reference_regime == re_regimes[1] & er$context == "LINEAR_CPM_SECONDARY" &
      er$method == "unico" & er$celltype == cells[1] & er$sample_id == samples[1] & er$gene_stratum == "ALL", ]
rmse <- sqrt(mean((old$mats$unico[g, 1, 1] - y[g, 1, 1])^2))
check("exact_original_rmse_sne", identical(row$rmse, rmse) && identical(row$sne, rmse / max(sd(y[g, 1, 1]), 1e-8)))
dr <- r$SCALE_SPECIFIC_DIRECTION_LONG
row <- dr[dr$reference_regime == re_regimes[1] & dr$context == "LINEAR_CPM_SECONDARY" & dr$method == "unico" &
      dr$celltype == cells[1] & dr$pathway == "signal", ]
delta <- rowMeans(old$mats$unico[g, 1, 3:4]) - rowMeans(old$mats$unico[g, 1, 1:2])
check("exact_original_projection_uses_entire_paired_delta_norm", identical(row$estimated_P,
      unname(math$projection(delta, g, "signal")["P"])))
check("original_point05_direction_threshold_unchanged", identical(unname(math$classify(c(-.051, -.05, 0, .05, .051))),
      c("negative", "zero", "zero", "zero", "positive")))
check("reference_and_analytic_zero_contrasts_retained", all(dr$estimated_P[dr$method %in% c("reference", "always_zero")] == 0))
check("missing_method_not_predicted_zero", all(is.na(dr$estimated_direction[dr$reference_regime == re_regimes[1] & dr$method == "tca"])))
bad <- regimes; bad[[2]]$mats$bulk[1] <- -1
check("unchanged_bulk_drift_blocks", fails(re_unit_endpoints(unit, genes, cells, samples, ev, bad, math)))
bad <- regimes; dimnames(bad[[2]]$mats$unico)[[2]] <- rev(cells)
check("axis_permutation_is_not_silently_relabelled", fails(re_unit_endpoints(unit, genes, cells, samples, ev, bad, math)))
bad <- regimes; bad[[2]]$mats$epic[1] <- NA_real_
check("nonfinite_canonical_blocks_instead_of_gene_dropping", fails(re_unit_endpoints(unit, genes, cells, samples, ev, bad, math)))
bad <- regimes; bad[[1]]$mats$unico <- bad[[1]]$mats$unico[1:2, , , drop = FALSE]
bad[[2]]$mats$unico <- bad[[2]]$mats$unico[10:11, , , drop = FALSE]
check("disjoint_support_blocks_without_relaxation", fails(re_support(bad, genes)))
check("constant_rank_is_undefined_not_zero", is.na(math$safe_rank(1:4, rep(1, 4))))
second <- re_unit_endpoints(unit, genes, cells, samples, ev, regimes, math)
check("deterministic_rebuild_all_tables_identical", identical(r, second))
for (n in names(r)) re_write(r[[n]], file.path(out, paste0(n, ".tsv")))
re_write(do.call(rbind, tests), file.path(out, "ARTIFICIAL_R_ENDPOINT_TESTS.tsv"))
sources <- c(args[1:2], normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
re_write(data.frame(path = sources, sha256 = vapply(sources, re_sha, character(1))),
         file.path(out, "ARTIFICIAL_TEST_SOURCE_HASHES.tsv"))
cat("ARTIFICIAL_REFERENCE_ENDPOINT_PASS", length(tests), "REAL_FITS_0_REAL_ENDPOINT_UNITS_0\n")
