# Executable workflow for parasite-host niche overlap analyses.
#
# The RMarkdown file calls run_niche_overlap_workflow() so the same code can be
# rendered interactively, sourced from R, or scheduled through targets. Expensive
# stages are cached in Results/cache with explicit keys.

run_niche_overlap_workflow <- function(project_dir, settings = list(), use_cache = TRUE) {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)

  function_file <- file.path(project_dir, "Scripts", "functions", "niche_overlap_functions.R")
  cache_file <- file.path(project_dir, "Scripts", "functions", "cache_helpers.R")
  source(function_file)
  source(cache_file)

  defaults <- list(
    occurrence_file = file.path(project_dir, "Data", "occurrences_bees_parasite_host.csv"),
    pair_file = file.path(project_dir, "Data", "parasite_host_pairs.csv"),
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
  settings <- utils::modifyList(defaults, as.list(settings))

  create_output_dirs(project_dir)

  cache_dir <- file.path(project_dir, "Results", "cache")
  tables_dir <- file.path(project_dir, "Results", "tables")
  figure_dir <- file.path(project_dir, "Figures")
  validation_dir <- file.path(project_dir, "Results", "validation")

  cache_status <- list()
  cache_version <- 5L

  validation_cache <- cache_eval(
    cache_dir = cache_dir,
    name = "occurrence_cleaning",
    key = list(
      cache_version = cache_version,
      occurrence = input_fingerprint(settings$occurrence_file),
      minimum_occurrences = settings$minimum_occurrences
    ),
    code = {
      raw_occurrences <- read_occurrence_file(settings$occurrence_file)
      validate_occurrence_data(
        raw_occurrences,
        min_occ = settings$minimum_occurrences,
        validation_dir = validation_dir
      )
    },
    use_cache = use_cache
  )
  validation <- validation_cache$value
  cache_status[["occurrence_cleaning"]] <- validation_cache[c("status", "path")]

  pair_cache <- cache_eval(
    cache_dir = cache_dir,
    name = "pair_definition",
    key = list(
      cache_version = cache_version,
      occurrence = input_fingerprint(settings$occurrence_file),
      pair_file = input_fingerprint(settings$pair_file),
      minimum_occurrences = settings$minimum_occurrences
    ),
    code = {
      build_pair_table(
        validation$clean,
        validation$species_summary,
        pair_file = settings$pair_file,
        min_occ = settings$minimum_occurrences
      )
    },
    use_cache = use_cache
  )
  pairs <- pair_cache$value
  cache_status[["pair_definition"]] <- pair_cache[c("status", "path")]

  variables <- load_worldclim_bioclim(project_dir, resolution = settings$worldclim_resolution)

  species_to_analyze <- unique(c(pairs$valid$parasite, pairs$valid$host))
  env_cache <- cache_eval(
    cache_dir = cache_dir,
    name = "environmental_extraction_accessible_areas",
    key = list(
      cache_version = cache_version,
      occurrence = input_fingerprint(settings$occurrence_file),
      species = sort(species_to_analyze),
      buffer_degrees = settings$buffer_degrees,
      worldclim_resolution = settings$worldclim_resolution,
      minimum_occurrences = settings$minimum_occurrences
    ),
    code = {
      prepare_environment_space(
        clean_occ = validation$clean,
        variables = variables,
        species_to_analyze = species_to_analyze,
        buffer_size = settings$buffer_degrees,
        min_occ = settings$minimum_occurrences
      )
    },
    use_cache = use_cache
  )
  env_space <- env_cache$value
  cache_status[["environmental_extraction_accessible_areas"]] <- env_cache[c("status", "path")]

  pairs_environment <- filter_pairs_by_environment(pairs$valid, env_space$species_summary)
  pairs$valid <- pairs_environment |>
    dplyr::filter(.data$environment_valid) |>
    dplyr::select(-environment_valid, -environment_skip_reason)
  pairs$skipped <- dplyr::bind_rows(
    pairs$skipped,
    pairs_environment |>
      dplyr::filter(!.data$environment_valid) |>
      dplyr::mutate(valid_pair = FALSE, skip_reason = .data$environment_skip_reason) |>
      dplyr::select(-environment_valid, -environment_skip_reason)
  )
  if (nrow(pairs$valid) == 0) {
    stop("No valid parasite-host pairs remain after occurrence and environmental filtering.", call. = FALSE)
  }

  pca_cache <- cache_eval(
    cache_dir = cache_dir,
    name = "environmental_pca",
    key = list(
      cache_version = cache_version,
      species = sort(names(env_space$spec_env)),
      background_cells = vapply(env_space$back_env, nrow, integer(1)),
      occurrence_cells = vapply(env_space$spec_env, nrow, integer(1))
    ),
    code = {
      run_environment_pca(
        spec_env = env_space$spec_env,
        back_env = env_space$back_env
      )
    },
    use_cache = use_cache
  )
  pca_results <- pca_cache$value
  cache_status[["environmental_pca"]] <- pca_cache[c("status", "path")]

  grid_cache <- cache_eval(
    cache_dir = cache_dir,
    name = "niche_density_grids",
    key = list(
      cache_version = cache_version,
      species = sort(names(pca_results$scores_spec)),
      grid_resolution = settings$grid_resolution,
      pca_eigenvalues = round(pca_results$eigenvalues$eigenvalue[seq_len(2)], 8)
    ),
    code = {
      build_density_grids(pca_results, grid_resolution = settings$grid_resolution)
    },
    use_cache = use_cache
  )
  grids <- grid_cache$value
  cache_status[["niche_density_grids"]] <- grid_cache[c("status", "path")]

  test_cache <- cache_eval(
    cache_dir = cache_dir,
    name = "overlap_dynamics_randomization",
    key = list(
      cache_version = cache_version,
      pairs = pairs$valid$pair_id,
      grid_resolution = settings$grid_resolution,
      repetitions = settings$test_repetitions,
      random_seed = settings$random_seed,
      randomization = settings$run_randomization_tests,
      future_workers = settings$future_workers
    ),
    code = {
      run_pairwise_niche_tests(
        valid_pairs = pairs$valid,
        grids = grids,
        repetitions = settings$test_repetitions,
        seed = settings$random_seed,
        ncores = settings$future_workers,
        run_randomization_tests = settings$run_randomization_tests,
        future_strategy = settings$future_strategy
      )
    },
    use_cache = use_cache
  )
  test_results <- test_cache$value
  cache_status[["overlap_dynamics_randomization"]] <- test_cache[c("status", "path")]

  run_settings <- tibble::tibble(
    minimum_occurrences = settings$minimum_occurrences,
    buffer_degrees = settings$buffer_degrees,
    worldclim_resolution_arcmin = settings$worldclim_resolution,
    grid_resolution = settings$grid_resolution,
    run_randomization_tests = settings$run_randomization_tests,
    test_repetitions = settings$test_repetitions,
    random_seed = settings$random_seed,
    future_workers = settings$future_workers,
    future_strategy = settings$future_strategy,
    cache_enabled = use_cache
  )

  save_analysis_outputs(
    project_dir = project_dir,
    validation = validation,
    pairs = pairs,
    env_space = env_space,
    pca_results = pca_results,
    grids = grids,
    test_results = test_results,
    run_settings = run_settings
  )
  modification_log <- write_modification_log(project_dir, run_settings, validation, pairs)

  expected_niche_files <- file.path(
    figure_dir,
    "Niche_Plots",
    paste(test_results$metrics$parasite, test_results$metrics$host, "niche_space.pdf", sep = "_")
  )
  figure_cache_allowed <- isTRUE(use_cache) && all(file.exists(expected_niche_files))

  figure_cache <- cache_eval(
    cache_dir = cache_dir,
    name = "figure_generation",
    key = list(
      cache_version = cache_version,
      pairs = test_results$metrics$pair_id,
      metrics = round(test_results$metrics$schoener_d, 8),
      figure_version = 11L
    ),
    code = {
      fig_occurrences <- plot_occurrence_map(validation$clean, env_space$world)
      fig_backgrounds <- plot_background_map(env_space$background_polygons, validation$clean, env_space$world)
      fig_overlap <- plot_overlap_heatmap(test_results$metrics)
      fig_dynamic_metrics <- plot_dynamic_metrics_stacked_bar(test_results$metrics)
      fig_loadings <- plot_pca_loadings(pca_results$loadings, pca_results$eigenvalues)
      top_pair <- test_results$metrics |>
        dplyr::slice_max(.data$schoener_d, n = 1, with_ties = FALSE)
      fig_top_dynamics <- plot_niche_dynamics_pair(top_pair, grids)

      save_ggplot_dual(fig_occurrences, file.path(figure_dir, "Figure1_cleaned_occurrences"), width = 10, height = 8)
      save_ggplot_dual(fig_backgrounds, file.path(figure_dir, "Figure2_accessible_area_backgrounds"), width = 10, height = 8)
      save_ggplot_dual(fig_overlap, file.path(figure_dir, "Figure3_schoener_d_heatmap"), width = 8.5, height = 5.5)
      save_ggplot_dual(fig_dynamic_metrics, file.path(figure_dir, "Figure4_niche_dynamics_metrics"), width = 8.6, height = 9.0)
      save_ggplot_dual(fig_top_dynamics, file.path(figure_dir, "Figure5_niche_dynamics_top_pair"), width = 7.5, height = 6.5)
      save_ggplot_dual(fig_loadings, file.path(figure_dir, "Figure6_pca_loadings"), width = 7, height = 6)
      niche_dynamic <- save_niche_dynamic_plots(test_results$metrics, grids, figure_dir, tables_dir = tables_dir)

      list(
        fig_occurrences = fig_occurrences,
        fig_backgrounds = fig_backgrounds,
        fig_overlap = fig_overlap,
        fig_dynamic_metrics = fig_dynamic_metrics,
        fig_top_dynamics = fig_top_dynamics,
        fig_loadings = fig_loadings,
        niche_dynamic = niche_dynamic
      )
    },
    use_cache = figure_cache_allowed
  )
  figures <- figure_cache$value
  cache_status[["figure_generation"]] <- figure_cache[c("status", "path")]

  cache_status_table <- dplyr::bind_rows(lapply(names(cache_status), function(stage) {
    absolute_path <- normalizePath(cache_status[[stage]]$path, mustWork = FALSE)
    project_prefix <- paste0(normalizePath(project_dir, mustWork = TRUE), .Platform$file.sep)
    relative_path <- if (startsWith(absolute_path, project_prefix)) {
      substring(absolute_path, nchar(project_prefix) + 1L)
    } else {
      absolute_path
    }
    tibble::tibble(
      stage = stage,
      status = cache_status[[stage]]$status,
      path = relative_path
    )
  }))
  readr::write_csv(cache_status_table, file.path(tables_dir, "cache_status.csv"))

  list(
    validation = validation,
    pairs = pairs,
    variables = variables,
    env_space = env_space,
    pca_results = pca_results,
    grids = grids,
    test_results = test_results,
    run_settings = run_settings,
    figures = figures,
    cache_status = cache_status_table,
    modification_log = modification_log
  )
}
