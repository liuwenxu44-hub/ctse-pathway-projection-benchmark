source("R/ctse_benchmark.R")
source("R/adapters.R")
source("R/simulation.R")
source("R/study_endpoints.R")

design <- ctse_simulation_design()
scenario <- "complete_cancellation"
profiles <- ctse_population_profiles(scenario, design)
bulk <- ctse_generate_bulk(
  scenario = scenario,
  replicate_id = 1L,
  design = design,
  profiles = profiles,
  samples_per_group = 4L
)
truth <- ctse_truth_tensor(profiles, bulk$metadata)
pathway <- ctse_simulation_pathways(design)$signal

result <- run_ctse_benchmark(
  estimate = truth,
  sample_metadata = bulk$metadata,
  pathway = pathway,
  group0 = "group0",
  group1 = "group1",
  truth = truth,
  comparison_contract = ctse_comparison_contract("probability", "probability",
                                                  "conditional_profile", "conditional_profile")
)

stopifnot(result$status == "SUCCESS")
print(result$projection)
