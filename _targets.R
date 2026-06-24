# targets pipeline for the parasite-host niche overlap workflow.
#
# Randomization is enabled by default here. Increase test_repetitions for
# higher p-value resolution when running final inferential analyses.

library(targets)
library(tarchetypes)

tar_option_set(
  packages = c(
    "ade4", "colorspace", "dplyr", "ecospat", "geodata", "ggplot2",
    "ggrepel", "kableExtra", "knitr", "patchwork", "purrr", "readr",
    "readxl", "rnaturalearth", "scales", "sessioninfo", "sf", "stringr",
    "terra", "tibble", "tidyr", "viridisLite", "future", "future.apply",
    "furrr"
  )
)

source("Scripts/functions/niche_overlap_functions.R")
source("Scripts/functions/cache_helpers.R")
source("Scripts/analysis/run_workflow.R")

list(
  tar_target(
    workflow_settings,
    list(
      minimum_occurrences = 5,
      buffer_degrees = 2,
      worldclim_resolution = 10,
      grid_resolution = 100,
      run_randomization_tests = TRUE,
      test_repetitions = 100,
      random_seed = 42,
      future_workers = 4,
      future_strategy = "multisession"
    )
  ),
  tar_target(
    workflow_results,
    run_niche_overlap_workflow(
      project_dir = getwd(),
      settings = workflow_settings,
      use_cache = TRUE
    )
  ),
  tar_render(
    supplementary_pdf,
    "Rmarkdown/parasite_host_bees.Rmd",
    output_file = "parasite_host_niche_overlap_supplement.pdf",
    output_dir = "Results"
  )
)
