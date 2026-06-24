# Helper functions for the parasite-host niche overlap workflow.
#
# The RMarkdown report sources this file so the document stays readable and the
# computational steps can be tested independently. The functions below are kept
# in pipeline order: project setup, data validation, pair construction,
# environmental-space preparation, niche metrics, plotting, and output writing.

# Create all folders required by the requested project layout.
create_output_dirs <- function(project_dir) {
  dirs <- file.path(
    project_dir,
    c(
      "Figures",
      "Figures/Niche_Plots",
      "Results",
      "Results/cache",
      "Results/cache/worldclim",
      "Results/objects",
      "Results/tables",
      "Results/validation",
      "Scripts",
      "Scripts/functions",
      "Scripts/analysis",
      "Scripts/visualization"
    )
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

# Convert underscore-separated species names to labels suitable for tables and
# figure facets. The occurrence file stores species as Genus_species.
species_label <- function(x) {
  gsub("_", " ", x, fixed = TRUE)
}

# Build filesystem-safe file-name fragments from species names or pair IDs.
file_stub <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  x <- gsub("^_+|_+$", "", x)
  x
}

# Parse numeric fields that may be stored as text, including decimal commas.
parse_numeric_column <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "na", "Na")] <- NA_character_
  suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
}

# Read CSV/TXT/XLS/XLSX occurrence-like files. CSV delimiters are detected from
# the header because the current occurrence CSV is semicolon-delimited.
read_occurrence_file <- function(path) {
  if (!file.exists(path)) {
    stop("Occurrence file not found: ", path, call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "txt")) {
    first_line <- readLines(path, n = 1, warn = FALSE)
    delim <- if (grepl(";", first_line, fixed = TRUE)) ";" else ","
    readr::read_delim(
      path,
      delim = delim,
      col_types = readr::cols(.default = readr::col_character()),
      trim_ws = TRUE,
      show_col_types = FALSE
    )
  } else if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(path, col_types = "text")
  } else {
    stop("Unsupported occurrence file extension: ", ext, call. = FALSE)
  }
}

# Validate and clean occurrence data before any environmental extraction.
# Corrections are written to disk so every removed or changed record is auditable.
validate_occurrence_data <- function(raw_occ, min_occ = 5, validation_dir) {
  required <- c("ID", "Species", "Interaction", "Longitude", "Latitude", "Year_collect")
  missing_cols <- setdiff(required, names(raw_occ))
  if (length(missing_cols) > 0) {
    stop(
      "Occurrence data are missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Standardize classes without trusting spreadsheet/CSV type guesses.
  raw <- tibble::as_tibble(raw_occ)
  prepared <- raw |>
    dplyr::mutate(
      row_number = dplyr::row_number(),
      ID = as.character(.data$ID),
      Species = stringr::str_squish(as.character(.data$Species)),
      Interaction = tolower(stringr::str_squish(as.character(.data$Interaction))),
      Longitude_original = as.character(.data$Longitude),
      Latitude_original = as.character(.data$Latitude),
      Longitude = parse_numeric_column(.data$Longitude),
      Latitude = parse_numeric_column(.data$Latitude),
      Year_collect = parse_numeric_column(.data$Year_collect)
    )

  # Keep issue tables separate so the report can state exactly what was removed.
  missing_species <- prepared |>
    dplyr::filter(is.na(.data$Species) | .data$Species == "") |>
    dplyr::mutate(issue = "missing_species")

  invalid_interaction <- prepared |>
    dplyr::filter(!.data$Interaction %in% c("host", "parasite")) |>
    dplyr::mutate(issue = "invalid_interaction")

  missing_coordinates <- prepared |>
    dplyr::filter(is.na(.data$Longitude) | is.na(.data$Latitude)) |>
    dplyr::mutate(issue = "missing_or_non_numeric_coordinates")

  invalid_coordinates <- prepared |>
    dplyr::filter(
      !is.na(.data$Longitude),
      !is.na(.data$Latitude),
      .data$Longitude < -180 | .data$Longitude > 180 |
        .data$Latitude < -90 | .data$Latitude > 90
    ) |>
    dplyr::mutate(issue = "coordinates_outside_valid_range")

  # Only valid geographic coordinates are allowed into the analytical dataset.
  valid_before_duplicates <- prepared |>
    dplyr::filter(
      !is.na(.data$Species),
      .data$Species != "",
      .data$Interaction %in% c("host", "parasite"),
      !is.na(.data$Longitude),
      !is.na(.data$Latitude),
      dplyr::between(.data$Longitude, -180, 180),
      dplyr::between(.data$Latitude, -90, 90)
    ) |>
    dplyr::mutate(
      coordinate_key = paste(
        .data$Species,
        .data$Interaction,
        round(.data$Longitude, 6),
        round(.data$Latitude, 6),
        sep = "|"
      )
    )

  # Duplicate removal is species- and interaction-specific. The first occurrence
  # is retained and later records at the same rounded coordinate are documented.
  duplicate_rows <- valid_before_duplicates |>
    dplyr::filter(duplicated(.data$coordinate_key)) |>
    dplyr::mutate(issue = "duplicate_species_interaction_coordinate")

  # The cleaned table is the single source used by all later workflow steps.
  clean <- valid_before_duplicates |>
    dplyr::filter(!duplicated(.data$coordinate_key)) |>
    dplyr::mutate(
      Species_display = species_label(.data$Species),
      Interaction = factor(.data$Interaction, levels = c("parasite", "host"))
    ) |>
    dplyr::select(
      row_number,
      ID,
      Species,
      Species_display,
      Interaction,
      Longitude,
      Latitude,
      Year_collect
    )

  issue_records <- dplyr::bind_rows(
    missing_species,
    invalid_interaction,
    missing_coordinates,
    invalid_coordinates,
    duplicate_rows
  ) |>
    dplyr::select(
      issue,
      row_number,
      ID,
      Species,
      Interaction,
      Longitude_original,
      Latitude_original,
      Longitude,
      Latitude,
      Year_collect
    )

  raw_counts <- prepared |>
    dplyr::count(.data$Interaction, .data$Species, name = "raw_records")

  clean_counts <- clean |>
    dplyr::mutate(Interaction = as.character(.data$Interaction)) |>
    dplyr::count(.data$Interaction, .data$Species, name = "clean_records")

  duplicate_counts <- duplicate_rows |>
    dplyr::count(.data$Interaction, .data$Species, name = "duplicates_removed")

  # Summaries are used both for reporting and for minimum-sample filtering.
  species_summary <- raw_counts |>
    dplyr::full_join(clean_counts, by = c("Interaction", "Species")) |>
    dplyr::full_join(duplicate_counts, by = c("Interaction", "Species")) |>
    dplyr::mutate(
      dplyr::across(
        c("raw_records", "clean_records", "duplicates_removed"),
        ~ tidyr::replace_na(.x, 0L)
      ),
      invalid_or_missing_removed = .data$raw_records -
        .data$clean_records - .data$duplicates_removed,
      minimum_records_required = min_occ,
      analysis_status = dplyr::if_else(
        .data$clean_records >= min_occ,
        "eligible",
        "below_minimum_records"
      )
    ) |>
    dplyr::arrange(.data$Interaction, .data$Species)

  corrections <- tibble::tibble(
    check = c(
      "delimiter_and_types",
      "missing_or_non_numeric_coordinates",
      "coordinates_outside_valid_range",
      "invalid_interaction_labels",
      "duplicate_species_interaction_coordinate",
      "minimum_occurrence_threshold"
    ),
    records_affected = c(
      nrow(prepared),
      nrow(missing_coordinates),
      nrow(invalid_coordinates),
      nrow(invalid_interaction),
      nrow(duplicate_rows),
      sum(species_summary$analysis_status == "below_minimum_records")
    ),
    correction_applied = c(
      "Read the semicolon-delimited file and parsed coordinates/year as numeric values.",
      "Removed records with coordinates that could not be parsed as numeric values.",
      "Removed records with longitude outside [-180, 180] or latitude outside [-90, 90].",
      "Removed records whose Interaction value was not host or parasite.",
      "Kept the first record per species, interaction, and rounded coordinate pair; removed later duplicates.",
      "Excluded species with fewer cleaned unique coordinates than the configured minimum from pairwise tests."
    )
  )

  validation_report <- issue_records |>
    dplyr::transmute(
      record_identifier = dplyr::if_else(
        is.na(.data$ID) | .data$ID == "",
        paste0("row_", .data$row_number),
        .data$ID
      ),
      species = .data$Species,
      interaction = .data$Interaction,
      validation_issue = .data$issue,
      action_performed = dplyr::case_when(
        .data$issue == "missing_species" ~ "Removed from analysis because species identity is missing.",
        .data$issue == "invalid_interaction" ~ "Removed from analysis because the interaction label is not host or parasite.",
        .data$issue == "missing_or_non_numeric_coordinates" ~ "Removed from analysis because longitude or latitude could not be parsed.",
        .data$issue == "coordinates_outside_valid_range" ~ "Removed from analysis because coordinates are outside valid longitude or latitude ranges.",
        .data$issue == "duplicate_species_interaction_coordinate" ~ "Removed as a duplicate; the first record for the same species, interaction, and rounded coordinate was retained.",
        TRUE ~ "Removed from analysis during occurrence validation."
      )
    )

  below_minimum_records <- clean |>
    dplyr::mutate(Interaction = as.character(.data$Interaction)) |>
    dplyr::inner_join(
      species_summary |>
        dplyr::filter(.data$analysis_status == "below_minimum_records") |>
        dplyr::select(Interaction, Species),
      by = c("Interaction", "Species")
    ) |>
    dplyr::transmute(
      record_identifier = dplyr::if_else(
        is.na(.data$ID) | .data$ID == "",
        paste0("row_", .data$row_number),
        .data$ID
      ),
      species = .data$Species,
      interaction = .data$Interaction,
      validation_issue = "species_below_minimum_cleaned_occurrences",
      action_performed = "Retained in cleaned occurrence data but excluded from pairwise niche analyses."
    )

  validation_report <- dplyr::bind_rows(validation_report, below_minimum_records) |>
    dplyr::arrange(.data$species, .data$interaction, .data$record_identifier)

  readr::write_csv(clean, file.path(validation_dir, "occurrences_cleaned.csv"))
  readr::write_csv(issue_records, file.path(validation_dir, "occurrence_issues.csv"))
  readr::write_csv(corrections, file.path(validation_dir, "validation_corrections.csv"))
  readr::write_csv(species_summary, file.path(validation_dir, "species_occurrence_summary.csv"))
  readr::write_csv(validation_report, file.path(validation_dir, "validation_report.csv"))

  list(
    raw = raw,
    prepared = prepared,
    clean = clean,
    issues = issue_records,
    corrections = corrections,
    validation_report = validation_report,
    species_summary = species_summary
  )
}

# Construct parasite-host candidate pairs. When a curated pair table exists, it
# is used; otherwise every parasite-host combination is considered a candidate.
build_pair_table <- function(clean_occ, species_summary, pair_file = NULL, min_occ = 5) {
  count_lookup <- species_summary |>
    dplyr::select(Interaction, Species, clean_records)

  if (!is.null(pair_file) && file.exists(pair_file)) {
    raw_pairs <- read_occurrence_file(pair_file)
    names(raw_pairs) <- tolower(names(raw_pairs))
    if (!all(c("parasite", "host") %in% names(raw_pairs))) {
      stop(
        "Pair file must contain at least 'parasite' and 'host' columns: ",
        pair_file,
        call. = FALSE
      )
    }

    candidate_pairs <- raw_pairs |>
      dplyr::transmute(
        parasite = gsub(" ", "_", stringr::str_squish(.data$parasite)),
        host = gsub(" ", "_", stringr::str_squish(.data$host)),
        pairing_rule = "user_supplied_pair_file"
      ) |>
      dplyr::distinct()
  } else {
    # No explicit interaction matrix is available in the provided data, so this
    # fallback preserves all candidate combinations for transparent screening.
    parasites <- count_lookup |>
      dplyr::filter(.data$Interaction == "parasite") |>
      dplyr::pull(.data$Species)
    hosts <- count_lookup |>
      dplyr::filter(.data$Interaction == "host") |>
      dplyr::pull(.data$Species)

    candidate_pairs <- tidyr::expand_grid(parasite = parasites, host = hosts) |>
      dplyr::mutate(pairing_rule = "all_candidate_parasite_host_combinations")
  }

  pair_counts <- candidate_pairs |>
    dplyr::left_join(
      count_lookup |>
        dplyr::filter(.data$Interaction == "parasite") |>
        dplyr::select(parasite = Species, n_parasite = clean_records),
      by = "parasite"
    ) |>
    dplyr::left_join(
      count_lookup |>
        dplyr::filter(.data$Interaction == "host") |>
        dplyr::select(host = Species, n_host = clean_records),
      by = "host"
    ) |>
    dplyr::mutate(
      n_parasite = tidyr::replace_na(.data$n_parasite, 0L),
      n_host = tidyr::replace_na(.data$n_host, 0L),
      pair_id = paste(file_stub(.data$parasite), file_stub(.data$host), sep = "__"),
      parasite_label = species_label(.data$parasite),
      host_label = species_label(.data$host),
      valid_pair = .data$n_parasite >= min_occ & .data$n_host >= min_occ,
      skip_reason = dplyr::case_when(
        .data$n_parasite < min_occ & .data$n_host < min_occ ~ "parasite_and_host_below_minimum_records",
        .data$n_parasite < min_occ ~ "parasite_below_minimum_records",
        .data$n_host < min_occ ~ "host_below_minimum_records",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::arrange(.data$parasite, .data$host)

  list(
    valid = pair_counts |> dplyr::filter(.data$valid_pair),
    skipped = pair_counts |> dplyr::filter(!.data$valid_pair),
    all = pair_counts
  )
}

# Download or reuse cached WorldClim bioclimatic variables. Files are cached
# inside Results/ so reruns do not repeatedly download the raster archive.
load_worldclim_bioclim <- function(project_dir, resolution = 10) {
  cache_dir <- file.path(project_dir, "Results", "cache", "worldclim")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  variables <- geodata::worldclim_global(var = "bio", res = resolution, path = cache_dir)
  names(variables) <- paste0("bio", seq_len(terra::nlyr(variables)))
  variables
}

# Create a species-specific accessible-area background from a buffered minimum
# convex polygon and clip it to land. This mirrors the original workflow logic.
make_mcp_background <- function(occ, buffer_size, land_union) {
  pts <- sf::st_as_sf(
    occ,
    coords = c("Longitude", "Latitude"),
    crs = 4326,
    remove = FALSE
  )
  # The original workflow used degree-based buffers; keep that behavior for
  # compatibility and document the limitation in the RMarkdown.
  hull <- suppressWarnings(sf::st_convex_hull(sf::st_union(sf::st_geometry(pts))))
  buffered <- suppressWarnings(sf::st_buffer(hull, dist = buffer_size))
  buffered <- sf::st_make_valid(sf::st_as_sf(data.frame(id = 1), geometry = buffered))
  buffered <- sf::st_set_crs(buffered, 4326)
  bg <- suppressWarnings(sf::st_intersection(buffered, land_union))
  bg <- sf::st_make_valid(bg)
  if (nrow(bg) == 0 || all(sf::st_is_empty(bg))) {
    stop("Background polygon is empty for ", unique(occ$Species), call. = FALSE)
  }
  bg
}

# Extract environmental values for each eligible species and its background.
# Species that lose too many records after raster extraction are flagged.
prepare_environment_space <- function(clean_occ, variables, species_to_analyze,
                                      buffer_size = 2, min_occ = 5) {
  # Planar operations are used because the original MCP/buffer logic is in
  # decimal degrees rather than projected meters.
  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  world <- sf::st_transform(world, 4326)
  land_union <- sf::st_as_sf(sf::st_union(sf::st_make_valid(world)))

  species_to_analyze <- unique(species_to_analyze)
  background_polygons <- list()
  back_env <- list()
  spec_env <- list()
  spec_occ_used <- list()

  for (sp_name in species_to_analyze) {
    occ_sp <- clean_occ |>
      dplyr::filter(.data$Species == sp_name) |>
      dplyr::arrange(.data$row_number)

    if (nrow(occ_sp) < min_occ) {
      next
    }

    bg <- make_mcp_background(occ_sp, buffer_size = buffer_size, land_union = land_union)
    bg$Species <- sp_name
    bg$Species_display <- species_label(sp_name)
    bg$Interaction <- as.character(unique(occ_sp$Interaction))
    background_polygons[[sp_name]] <- bg

    # terra::extract returns an ID column first; remove it and keep complete
    # bioclimatic rows only.
    extracted_bg <- terra::extract(variables, terra::vect(bg))[, -1, drop = FALSE]
    extracted_bg <- as.data.frame(extracted_bg)
    extracted_bg <- extracted_bg[stats::complete.cases(extracted_bg), , drop = FALSE]

    extracted_occ <- terra::extract(variables, occ_sp[, c("Longitude", "Latitude")])[, -1, drop = FALSE]
    extracted_occ <- as.data.frame(extracted_occ)
    complete_occ <- stats::complete.cases(extracted_occ)
    extracted_occ <- extracted_occ[complete_occ, , drop = FALSE]
    occ_sp <- occ_sp[complete_occ, , drop = FALSE]

    back_env[[sp_name]] <- extracted_bg
    spec_env[[sp_name]] <- extracted_occ
    spec_occ_used[[sp_name]] <- occ_sp
  }

  env_species_summary <- tibble::tibble(
    Species = names(spec_env),
    Species_display = species_label(names(spec_env)),
    occurrence_records_for_environment = vapply(spec_env, nrow, integer(1)),
    background_cells = vapply(back_env, nrow, integer(1)),
    environmental_status = dplyr::case_when(
      .data$occurrence_records_for_environment < min_occ ~ "below_minimum_after_environment_extraction",
      .data$background_cells == 0 ~ "no_background_environment_cells",
      TRUE ~ "eligible"
    )
  )

  eligible_species <- env_species_summary |>
    dplyr::filter(.data$environmental_status == "eligible") |>
    dplyr::pull(.data$Species)

  list(
    world = world,
    background_polygons = background_polygons[eligible_species],
    back_env = back_env[eligible_species],
    spec_env = spec_env[eligible_species],
    spec_occ_used = spec_occ_used[eligible_species],
    species_summary = env_species_summary
  )
}

# Fit the Broennimann-style global PCA. Background cells receive weight 1 and
# occurrence rows receive weight 0, matching the original implementation.
run_environment_pca <- function(spec_env, back_env) {
  species_order <- names(spec_env)
  all_spec_env <- dplyr::bind_rows(spec_env)
  all_back_env <- dplyr::bind_rows(back_env)
  data_env <- dplyr::bind_rows(all_spec_env, all_back_env)

  weights <- c(rep(0, nrow(all_spec_env)), rep(1, nrow(all_back_env)))
  pca <- ade4::dudi.pca(
    data_env,
    row.w = weights,
    center = TRUE,
    scale = TRUE,
    scannf = FALSE,
    nf = 2
  )

  spec_counts <- vapply(spec_env, nrow, integer(1))
  back_counts <- vapply(back_env, nrow, integer(1))
  spec_offsets <- cumsum(c(0L, spec_counts))
  back_offsets <- cumsum(c(0L, back_counts))
  first_back_row <- nrow(all_spec_env)

  scores_spec <- list()
  scores_back <- list()
  # Split the global PCA scores back into per-species occurrence and background
  # tables required by ecospat.grid.clim.dyn().
  for (i in seq_along(species_order)) {
    sp_name <- species_order[[i]]
    spec_idx <- seq.int(spec_offsets[[i]] + 1L, spec_offsets[[i + 1L]])
    back_idx <- seq.int(
      first_back_row + back_offsets[[i]] + 1L,
      first_back_row + back_offsets[[i + 1L]]
    )
    scores_spec[[sp_name]] <- as.data.frame(pca$li[spec_idx, 1:2, drop = FALSE])
    scores_back[[sp_name]] <- as.data.frame(pca$li[back_idx, 1:2, drop = FALSE])
  }

  spec_scores_table <- dplyr::bind_rows(
    lapply(names(scores_spec), function(sp_name) {
      dplyr::bind_cols(
        tibble::tibble(Species = sp_name, Species_display = species_label(sp_name), score_type = "occurrence"),
        scores_spec[[sp_name]]
      )
    })
  )

  back_scores_table <- dplyr::bind_rows(
    lapply(names(scores_back), function(sp_name) {
      dplyr::bind_cols(
        tibble::tibble(Species = sp_name, Species_display = species_label(sp_name), score_type = "background"),
        scores_back[[sp_name]]
      )
    })
  )

  loadings <- as.data.frame(pca$co[, 1:2, drop = FALSE])
  loadings$variable <- rownames(loadings)
  names(loadings)[1:2] <- c("axis1", "axis2")

  eigenvalues <- tibble::tibble(
    axis = seq_along(pca$eig),
    eigenvalue = pca$eig,
    variance_percent = 100 * pca$eig / sum(pca$eig)
  )

  list(
    pca = pca,
    all_spec_env = all_spec_env,
    all_back_env = all_back_env,
    data_env = data_env,
    scores_spec = scores_spec,
    scores_back = scores_back,
    total_scores_back = dplyr::bind_rows(scores_back),
    spec_scores_table = spec_scores_table,
    back_scores_table = back_scores_table,
    loadings = loadings,
    eigenvalues = eigenvalues
  )
}

# Build ecospat density grids in the two-dimensional PCA environmental space.
build_density_grids <- function(pca_results, grid_resolution = 100) {
  species_order <- names(pca_results$scores_spec)
  grids <- list()
  for (sp_name in species_order) {
    grids[[sp_name]] <- ecospat::ecospat.grid.clim.dyn(
      glob = pca_results$total_scores_back,
      glob1 = pca_results$scores_back[[sp_name]],
      sp = pca_results$scores_spec[[sp_name]],
      R = grid_resolution
    )
  }
  grids
}

# Remove pairs where either species failed environmental extraction.
filter_pairs_by_environment <- function(valid_pairs, env_species_summary) {
  eligible_species <- env_species_summary |>
    dplyr::filter(.data$environmental_status == "eligible") |>
    dplyr::pull(.data$Species)

  valid_pairs |>
    dplyr::mutate(
      environment_valid = .data$parasite %in% eligible_species & .data$host %in% eligible_species,
      environment_skip_reason = dplyr::case_when(
        !.data$parasite %in% eligible_species & !.data$host %in% eligible_species ~ "parasite_and_host_failed_environment_extraction",
        !.data$parasite %in% eligible_species ~ "parasite_failed_environment_extraction",
        !.data$host %in% eligible_species ~ "host_failed_environment_extraction",
        TRUE ~ NA_character_
      )
    )
}

# Calculate pairwise overlap metrics. Randomization tests are optional because
# equivalency and similarity tests are computationally expensive for many pairs.
run_pairwise_niche_tests <- function(valid_pairs, grids, repetitions = 100,
                                     seed = 42, ncores = 1,
                                     run_randomization_tests = TRUE,
                                     future_strategy = "multisession") {
  ncores <- max(1L, as.integer(ncores))
  randomization_enabled <- isTRUE(run_randomization_tests) && repetitions > 0
  if (nrow(valid_pairs) == 0) {
    return(list(
      metrics = tibble::tibble(),
      equivalency_null = tibble::tibble(),
      similarity_null = tibble::tibble(),
      errors = tibble::tibble()
    ))
  }

  analyze_pair <- function(i) {
    pair <- valid_pairs[i, , drop = FALSE]
    set.seed(seed + i)

    tryCatch({
      z_parasite <- grids[[pair$parasite]]
      z_host <- grids[[pair$host]]

      overlap <- ecospat::ecospat.niche.overlap(z_parasite, z_host, cor = TRUE)
      dyn_ph <- ecospat::ecospat.niche.dyn.index(z_parasite, z_host)$dynamic.index.w
      dyn_hp <- ecospat::ecospat.niche.dyn.index(z_host, z_parasite)$dynamic.index.w

      if (randomization_enabled) {
        # These three tests generate null distributions by randomization.
        ecospat_cores <- ncores
        equivalency <- ecospat::ecospat.niche.equivalency.test(
          z_parasite,
          z_host,
          rep = repetitions,
          ncores = ecospat_cores
        )
        similarity_ph <- ecospat::ecospat.niche.similarity.test(
          z_parasite,
          z_host,
          rep = repetitions,
          ncores = ecospat_cores
        )
        similarity_hp <- ecospat::ecospat.niche.similarity.test(
          z_host,
          z_parasite,
          rep = repetitions,
          ncores = ecospat_cores
        )

        equivalency_p_d <- equivalency$p.D
        equivalency_p_i <- equivalency$p.I
        similarity_p_d_ph <- similarity_ph$p.D
        similarity_p_d_hp <- similarity_hp$p.D
        similarity_p_i_ph <- similarity_ph$p.I
        similarity_p_i_hp <- similarity_hp$p.I
      } else {
        # Keep the output schema stable when randomization is disabled.
        equivalency <- similarity_ph <- similarity_hp <- NULL
        equivalency_p_d <- equivalency_p_i <- NA_real_
        similarity_p_d_ph <- similarity_p_d_hp <- NA_real_
        similarity_p_i_ph <- similarity_p_i_hp <- NA_real_
      }

      metrics <- tibble::tibble(
        pair_id = pair$pair_id,
        parasite = pair$parasite,
        host = pair$host,
        parasite_label = pair$parasite_label,
        host_label = pair$host_label,
        n_parasite = pair$n_parasite,
        n_host = pair$n_host,
        schoener_d = unname(overlap[["D"]]),
        warren_i = unname(overlap[["I"]]),
        randomization_tests_run = randomization_enabled,
        equivalency_p_d = equivalency_p_d,
        equivalency_p_i = equivalency_p_i,
        similarity_p_d_parasite_to_host = similarity_p_d_ph,
        similarity_p_d_host_to_parasite = similarity_p_d_hp,
        similarity_p_i_parasite_to_host = similarity_p_i_ph,
        similarity_p_i_host_to_parasite = similarity_p_i_hp,
        # Parasite-host terminology. The host is the reference niche and the
        # parasite is the focal niche, matching ecospat's z1/reference and
        # z2/focal convention for dynamic indices.
        parasite_host_shared_stability = unname(dyn_hp[["stability"]]),
        parasite_exclusive_environment_use = unname(dyn_hp[["expansion"]]),
        host_environment_unfilled_by_parasite = unname(dyn_hp[["unfilling"]]),
        # Reverse-direction values are retained for auditability.
        host_parasite_shared_stability = unname(dyn_ph[["stability"]]),
        host_exclusive_environment_use = unname(dyn_ph[["expansion"]]),
        parasite_environment_unfilled_by_host = unname(dyn_ph[["unfilling"]]),
        # Backward-compatible columns retained for older downstream outputs.
        expansion_parasite_vs_host = unname(dyn_ph[["expansion"]]),
        stability_parasite_vs_host = unname(dyn_ph[["stability"]]),
        unfilling_parasite_vs_host = unname(dyn_ph[["unfilling"]]),
        expansion_host_vs_parasite = unname(dyn_hp[["expansion"]]),
        stability_host_vs_parasite = unname(dyn_hp[["stability"]]),
        unfilling_host_vs_parasite = unname(dyn_hp[["unfilling"]])
      )

      if (randomization_enabled) {
        equivalency_null <- equivalency$sim |>
          tibble::as_tibble() |>
          dplyr::mutate(
            pair_id = pair$pair_id,
            parasite = pair$parasite,
            host = pair$host,
            replicate = dplyr::row_number(),
            test = "equivalency"
          ) |>
          dplyr::relocate(pair_id, parasite, host, test, replicate)

        similarity_null <- dplyr::bind_rows(
          similarity_ph$sim |>
            tibble::as_tibble() |>
            dplyr::mutate(
              pair_id = pair$pair_id,
              parasite = pair$parasite,
              host = pair$host,
              direction = "parasite_to_host",
              replicate = dplyr::row_number()
            ),
          similarity_hp$sim |>
            tibble::as_tibble() |>
            dplyr::mutate(
              pair_id = pair$pair_id,
              parasite = pair$parasite,
              host = pair$host,
              direction = "host_to_parasite",
              replicate = dplyr::row_number()
            )
        ) |>
          dplyr::relocate(pair_id, parasite, host, direction, replicate)
      } else {
        equivalency_null <- tibble::tibble()
        similarity_null <- tibble::tibble()
      }

      list(
        metrics = metrics,
        equivalency_null = equivalency_null,
        similarity_null = similarity_null,
        errors = tibble::tibble()
      )
    }, error = function(e) {
      list(
        metrics = tibble::tibble(),
        equivalency_null = tibble::tibble(),
        similarity_null = tibble::tibble(),
        errors = tibble::tibble(
          pair_id = pair$pair_id,
          parasite = pair$parasite,
          host = pair$host,
          error_message = conditionMessage(e)
        )
      )
    })
  }

  pair_indices <- seq_len(nrow(valid_pairs))
  parallelize_pairs <- !randomization_enabled &&
    ncores > 1L &&
    requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)

  if (parallelize_pairs) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future_strategy, workers = ncores)
    pair_results <- future.apply::future_lapply(
      pair_indices,
      analyze_pair,
      future.seed = seed
    )
  } else {
    pair_results <- lapply(pair_indices, analyze_pair)
  }

  metrics <- dplyr::bind_rows(lapply(pair_results, `[[`, "metrics"))
  equivalency_null <- dplyr::bind_rows(lapply(pair_results, `[[`, "equivalency_null"))
  similarity_null <- dplyr::bind_rows(lapply(pair_results, `[[`, "similarity_null"))
  errors <- dplyr::bind_rows(lapply(pair_results, `[[`, "errors"))

  if (nrow(metrics) == 0 && nrow(errors) > 0) {
    stop(
      "All pairwise niche tests failed. First error: ",
      errors$error_message[[1]],
      call. = FALSE
    )
  }

  list(
    metrics = metrics,
    equivalency_null = equivalency_null,
    similarity_null = similarity_null,
    errors = errors
  )
}

# Shared ggplot theme for publication-oriented figures.
theme_publication <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "grey88"),
      strip.background = ggplot2::element_rect(fill = "grey94", colour = "grey70"),
      strip.text = ggplot2::element_text(face = "italic"),
      axis.title = ggplot2::element_text(colour = "grey15"),
      axis.text = ggplot2::element_text(colour = "grey20"),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(colour = "grey25")
    )
}

# Save each ggplot as PDF for vector output and TIFF for journal-style raster
# submission. TIFFs use LZW compression to keep file sizes manageable.
save_ggplot_dual <- function(plot, filename_base, width, height, dpi = 600) {
  ggplot2::ggsave(
    paste0(filename_base, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = "pdf"
  )
  ggplot2::ggsave(
    paste0(filename_base, ".tiff"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    compression = "lzw"
  )
}

# Map all cleaned occurrences after validation and duplicate filtering.
plot_occurrence_map <- function(clean_occ, world) {
  xlim <- range(clean_occ$Longitude, na.rm = TRUE) + c(-4, 4)
  ylim <- range(clean_occ$Latitude, na.rm = TRUE) + c(-4, 4)

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = world, fill = "grey96", colour = "grey78", linewidth = 0.15) +
    ggplot2::geom_point(
      data = clean_occ,
      ggplot2::aes(
        x = .data$Longitude,
        y = .data$Latitude,
        colour = .data$Interaction,
        shape = .data$Interaction
      ),
      size = 1.6,
      alpha = 0.82,
      stroke = 0.2
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$Species_display), ncol = 4) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::scale_colour_manual(
      values = c(parasite = "#D55E00", host = "#0072B2"),
      labels = c(parasite = "Parasite", host = "Host")
    ) +
    ggplot2::scale_shape_manual(
      values = c(parasite = 17, host = 16),
      labels = c(parasite = "Parasite", host = "Host")
    ) +
    ggplot2::labs(
      x = "Longitude",
      y = "Latitude",
      colour = "Interaction",
      shape = "Interaction",
      title = "Validated occurrence records"
    ) +
    theme_publication(base_size = 9.5) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.text = ggplot2::element_text(face = "italic", size = 8)
    )
}

# Map the MCP + buffer accessible-area polygons against cleaned occurrences.
plot_background_map <- function(background_polygons, clean_occ, world) {
  bg <- do.call(rbind, background_polygons)
  xlim <- range(clean_occ$Longitude, na.rm = TRUE) + c(-4, 4)
  ylim <- range(clean_occ$Latitude, na.rm = TRUE) + c(-4, 4)

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = world, fill = "grey97", colour = "grey82", linewidth = 0.15) +
    ggplot2::geom_sf(
      data = bg,
      ggplot2::aes(fill = .data$Interaction),
      alpha = 0.22,
      colour = "grey25",
      linewidth = 0.18
    ) +
    ggplot2::geom_point(
      data = clean_occ |> dplyr::filter(.data$Species %in% unique(bg$Species)),
      ggplot2::aes(
        x = .data$Longitude,
        y = .data$Latitude,
        colour = .data$Interaction,
        shape = .data$Interaction
      ),
      size = 0.8,
      alpha = 0.85
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$Species_display), ncol = 4) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::scale_fill_manual(
      values = c(parasite = "#D55E00", host = "#0072B2"),
      labels = c(parasite = "Parasite", host = "Host")
    ) +
    ggplot2::scale_colour_manual(
      values = c(parasite = "#D55E00", host = "#0072B2"),
      labels = c(parasite = "Parasite", host = "Host")
    ) +
    ggplot2::scale_shape_manual(
      values = c(parasite = 17, host = 16),
      labels = c(parasite = "Parasite", host = "Host")
    ) +
    ggplot2::labs(
      x = "Longitude",
      y = "Latitude",
      fill = "Interaction",
      colour = "Interaction",
      shape = "Interaction",
      title = "Accessible-area backgrounds",
      subtitle = "Minimum convex polygons buffered in degrees and clipped to land"
    ) +
    theme_publication(base_size = 9.5) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.text = ggplot2::element_text(face = "italic", size = 8)
    )
}

# Heatmap of pairwise Schoener's D values.
plot_overlap_heatmap <- function(metrics) {
  ggplot2::ggplot(
    metrics,
    ggplot2::aes(x = .data$host_label, y = .data$parasite_label, fill = .data$schoener_d)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$schoener_d)),
      size = 3,
      colour = "grey10"
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      limits = c(0, 1),
      name = "Schoener's D"
    ) +
    ggplot2::labs(
      x = "Host bee species",
      y = "Parasite species",
      title = "Environmental niche overlap"
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "italic"),
      axis.text.y = ggplot2::element_text(face = "italic"),
      panel.grid = ggplot2::element_blank()
    )
}

# Correlation-circle style PCA loading plot for the 19 bioclimatic variables.
plot_pca_loadings <- function(loadings, eigenvalues) {
  axis1 <- eigenvalues$variance_percent[[1]]
  axis2 <- eigenvalues$variance_percent[[2]]

  ggplot2::ggplot(loadings, ggplot2::aes(x = .data$axis1, y = .data$axis2)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey72", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey72", linewidth = 0.25) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, y = 0, xend = .data$axis1, yend = .data$axis2),
      arrow = ggplot2::arrow(length = grid::unit(0.14, "cm")),
      linewidth = 0.35,
      colour = "#333333"
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$variable),
      size = 3,
      min.segment.length = 0,
      box.padding = 0.25,
      seed = 42
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", axis1),
      y = sprintf("PC2 (%.1f%%)", axis2),
      title = "Bioclimatic variable loadings"
    ) +
    theme_publication(base_size = 10)
}

# Convert a SpatRaster with x/y cells into a data frame for ggplot2.
raster_to_plot_df <- function(raster, value_name = "value") {
  df <- terra::as.data.frame(raster, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- value_name
  df
}

# Prepare a long-format table of the parasite-host dynamic metrics used in the
# report. The labels avoid invasive/native wording while preserving the original
# mathematical definitions from Broennimann-style dynamic indices.
make_dynamic_metrics_long <- function(metrics) {
  metrics |>
    dplyr::select(
      pair_id,
      parasite_label,
      host_label,
      shared = parasite_host_shared_stability,
      parasite_exclusive = parasite_exclusive_environment_use,
      host_unfilled = host_environment_unfilled_by_parasite
    ) |>
    tidyr::pivot_longer(
      cols = c("shared", "parasite_exclusive", "host_unfilled"),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::mutate(
      metric_label = dplyr::case_when(
        .data$metric == "shared" ~ "Shared parasite-host environment",
        .data$metric == "parasite_exclusive" ~ "Parasite-exclusive environment",
        .data$metric == "host_unfilled" ~ "Host environment unfilled by parasite",
        TRUE ~ .data$metric
      ),
      pair_label = paste0(.data$parasite_label, " x ", .data$host_label)
    )
}

# Prepare a normalized visual partition for Figure 4. The exact dynamic metrics
# remain in the tables; this derived table rescales the three displayed
# components so each stacked bar sums to one.
make_dynamic_visual_partition <- function(metrics) {
  metrics |>
    dplyr::transmute(
      pair_id = .data$pair_id,
      pair_label = paste0(.data$parasite_label, " x ", .data$host_label),
      schoener_d = .data$schoener_d,
      shared = .data$parasite_host_shared_stability,
      parasite_exclusive = .data$parasite_exclusive_environment_use,
      host_unfilled = .data$host_environment_unfilled_by_parasite
    ) |>
    tidyr::pivot_longer(
      cols = c("shared", "parasite_exclusive", "host_unfilled"),
      names_to = "component",
      values_to = "raw_value"
    ) |>
    dplyr::group_by(.data$pair_id) |>
    dplyr::mutate(
      component_sum = sum(.data$raw_value, na.rm = TRUE),
      visual_value = dplyr::if_else(.data$component_sum > 0, .data$raw_value / .data$component_sum, NA_real_)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      component_label = factor(
        dplyr::case_when(
          .data$component == "shared" ~ "Shared parasite-host environment",
          .data$component == "parasite_exclusive" ~ "Parasite-exclusive environment",
          .data$component == "host_unfilled" ~ "Host environment unfilled by parasite",
          TRUE ~ .data$component
        ),
        levels = c(
          "Shared parasite-host environment",
          "Parasite-exclusive environment",
          "Host environment unfilled by parasite"
        )
      )
    )
}

# Figure-level summary of the three parasite-host dynamic components. The plot is
# a stacked-bar replacement for the earlier heatmap and is sorted by Schoener's D.
plot_dynamic_metrics_stacked_bar <- function(metrics) {
  dynamic_long <- make_dynamic_visual_partition(metrics)
  pair_order <- metrics |>
    dplyr::arrange(dplyr::desc(.data$schoener_d)) |>
    dplyr::transmute(pair_label = paste0(.data$parasite_label, " x ", .data$host_label)) |>
    dplyr::pull(.data$pair_label)

  ggplot2::ggplot(
    dynamic_long,
    ggplot2::aes(
      y = factor(.data$pair_label, levels = rev(pair_order)),
      x = .data$visual_value,
      fill = .data$component_label
    )
  ) +
    ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.18) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Shared parasite-host environment" = "#1B9E77",
        "Parasite-exclusive environment" = "#D95F02",
        "Host environment unfilled by parasite" = "#1F78B4"
      ),
      name = "Component"
    ) +
    ggplot2::labs(
      x = "Normalized contribution to displayed dynamics components",
      y = "Parasite-host pair",
      title = "Parasite-host niche dynamics partition",
      subtitle = "Pairs are sorted by Schoener's D; exact metric values are reported in tables"
    ) +
    theme_publication(base_size = 8) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "italic", size = 6.8),
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank()
    )
}

# Backward-compatible alias for older reports that used the previous function
# name before Figure 4 was changed from a heatmap to stacked bars.
plot_dynamic_metrics_heatmap <- plot_dynamic_metrics_stacked_bar

# ggplot2 reconstruction of the ecospat dynamic-category map. The host is the
# reference niche and the parasite is the focal niche.
plot_niche_dynamics_pair <- function(row, grids) {
  z_host <- grids[[row$host]]
  z_parasite <- grids[[row$parasite]]
  dynamic_raster <- ecospat::ecospat.niche.dyn.index(z_host, z_parasite)$dyn

  dynamic_df <- raster_to_plot_df(dynamic_raster, "category") |>
    dplyr::mutate(
      category = as.integer(.data$category),
      category_label = factor(
        dplyr::case_when(
          .data$category == 2L ~ "Host unfilled by parasite",
          .data$category == 3L ~ "Shared environment",
          .data$category == 4L ~ "Parasite-exclusive environment",
          .data$category == 1L ~ "Host-only non-analog",
          .data$category == 5L ~ "Parasite-only non-analog",
          .data$category == 6L ~ "Other available environment",
          TRUE ~ NA_character_
        ),
        levels = c(
          "Shared environment",
          "Parasite-exclusive environment",
          "Host unfilled by parasite",
          "Host-only non-analog",
          "Parasite-only non-analog",
          "Other available environment"
        )
      )
    ) |>
    dplyr::filter(!is.na(.data$category_label))

  host_contours <- niche_contour_df(z_host$z.uncor)
  parasite_contours <- niche_contour_df(z_parasite$z.uncor)
  host_half_density <- half_density_break(host_contours$density)
  parasite_half_density <- half_density_break(parasite_contours$density)

  ggplot2::ggplot(dynamic_df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$category_label), alpha = 0.95) +
    ggplot2::geom_contour(
      data = host_contours,
      ggplot2::aes(z = .data$available, colour = "Host", linetype = "Full available niche"),
      breaks = 0.5,
      linewidth = 0.48
    ) +
    ggplot2::geom_contour(
      data = host_contours,
      ggplot2::aes(z = .data$density, colour = "Host", linetype = "50% density"),
      breaks = host_half_density,
      linewidth = 0.48
    ) +
    ggplot2::geom_contour(
      data = parasite_contours,
      ggplot2::aes(z = .data$available, colour = "Parasite", linetype = "Full available niche"),
      breaks = 0.5,
      linewidth = 0.48
    ) +
    ggplot2::geom_contour(
      data = parasite_contours,
      ggplot2::aes(z = .data$density, colour = "Parasite", linetype = "50% density"),
      breaks = parasite_half_density,
      linewidth = 0.48
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Shared environment" = "#1B9E77",
        "Parasite-exclusive environment" = "#D95F02",
        "Host unfilled by parasite" = "#1F78B4",
        "Host-only non-analog" = "#D6EAF8",
        "Parasite-only non-analog" = "#FDE0C5",
        "Other available environment" = "#F2F2F2"
      ),
      drop = TRUE,
      name = "Dynamic category"
    ) +
    ggplot2::scale_colour_manual(
      values = c("Host" = "#08519C", "Parasite" = "#E6550D"),
      name = "Niche outline"
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Full available niche" = "solid", "50% density" = "22"),
      name = "Contour"
    ) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(
      x = "PC1",
      y = "PC2",
      title = paste(row$parasite_label, "x", row$host_label),
      subtitle = sprintf(
        "D = %.3f; shared = %.3f; parasite-exclusive = %.3f; host-unfilled = %.3f",
        row$schoener_d,
        row$parasite_host_shared_stability,
        row$parasite_exclusive_environment_use,
        row$host_environment_unfilled_by_parasite
      )
    ) +
    theme_publication(base_size = 9) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      panel.grid = ggplot2::element_blank()
    )
}

# Prepare a binary full-available-niche surface and a density surface for the
# two contour lines requested in pair-specific niche dynamics figures.
niche_contour_df <- function(raster) {
  raster_to_plot_df(raster, "density") |>
    dplyr::mutate(
      density = tidyr::replace_na(.data$density, 0),
      available = as.integer(.data$density > 0)
    )
}

# Use one dashed contour at 50% of the maximum density. Empty or flat density
# rasters return no breaks, allowing the dynamic map to render without warnings.
half_density_break <- function(x) {
  max_density <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(max_density) || max_density <= 0) {
    return(numeric(0))
  }
  max_density * 0.5
}

# Save ggplot2 niche dynamic plots for every analyzed pair.
save_niche_dynamic_plots <- function(metrics, grids, figure_dir, tables_dir = NULL) {
  niche_dir <- file.path(figure_dir, "Niche_Plots")
  dir.create(niche_dir, recursive = TRUE, showWarnings = FALSE)

  multipage_pdf <- file.path(figure_dir, "Figure5_niche_dynamics_all_pairs.pdf")
  grDevices::pdf(multipage_pdf, width = 7, height = 6.5, onefile = TRUE)
  for (i in seq_len(nrow(metrics))) {
    print(plot_niche_dynamics_pair(metrics[i, ], grids))
  }
  grDevices::dev.off()

  for (i in seq_len(nrow(metrics))) {
    row <- metrics[i, ]
    plot <- plot_niche_dynamics_pair(row, grids)
    filename_base <- paste(row$parasite, row$host, "niche_space", sep = "_")
    ggplot2::ggsave(
      file.path(niche_dir, paste0(filename_base, ".pdf")),
      plot = plot,
      width = 7,
      height = 6.5,
      device = "pdf"
    )
    ggplot2::ggsave(
      file.path(niche_dir, paste0(filename_base, ".png")),
      plot = plot,
      width = 7,
      height = 6.5,
      dpi = 600
    )
  }

  figure_summary <- metrics |>
    dplyr::transmute(
      Parasite = .data$parasite_label,
      Host = .data$host_label,
      D = .data$schoener_d,
      I = .data$warren_i,
      Figure_filename = file.path(
        "Figures",
        "Niche_Plots",
        paste(.data$parasite, .data$host, "niche_space.pdf", sep = "_")
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$D))

  if (!is.null(tables_dir)) {
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(figure_summary, file.path(tables_dir, "niche_space_figures.csv"))
  }

  list(multipage_pdf = multipage_pdf, summary = figure_summary)
}

# Write every reproducibility artifact requested by the project specification.
save_analysis_outputs <- function(project_dir, validation, pairs, env_space,
                                  pca_results, grids, test_results,
                                  run_settings) {
  tables_dir <- file.path(project_dir, "Results", "tables")
  objects_dir <- file.path(project_dir, "Results", "objects")
  validation_dir <- file.path(project_dir, "Results", "validation")

  readr::write_csv(validation$clean, file.path(tables_dir, "occurrences_cleaned.csv"))
  readr::write_csv(validation$species_summary, file.path(tables_dir, "species_occurrence_summary.csv"))
  readr::write_csv(validation$corrections, file.path(tables_dir, "validation_corrections.csv"))
  readr::write_csv(validation$issues, file.path(validation_dir, "occurrence_issues.csv"))
  readr::write_csv(validation$validation_report, file.path(validation_dir, "validation_report.csv"))
  readr::write_csv(pairs$all, file.path(tables_dir, "parasite_host_pairs_all.csv"))
  readr::write_csv(pairs$valid, file.path(tables_dir, "parasite_host_pairs_used.csv"))
  readr::write_csv(pairs$skipped, file.path(tables_dir, "parasite_host_pairs_skipped.csv"))
  readr::write_csv(env_space$species_summary, file.path(tables_dir, "environment_extraction_summary.csv"))
  readr::write_csv(pca_results$loadings, file.path(tables_dir, "pca_loadings.csv"))
  readr::write_csv(pca_results$eigenvalues, file.path(tables_dir, "pca_eigenvalues.csv"))
  readr::write_csv(pca_results$spec_scores_table, file.path(tables_dir, "pca_scores_occurrences.csv"))
  readr::write_csv(pca_results$back_scores_table, file.path(tables_dir, "pca_scores_background.csv"))
  readr::write_csv(test_results$metrics, file.path(tables_dir, "niche_overlap_metrics.csv"))
  readr::write_csv(
    make_dynamic_metrics_long(test_results$metrics),
    file.path(tables_dir, "niche_dynamics_metrics_long.csv")
  )
  readr::write_csv(test_results$equivalency_null, file.path(tables_dir, "niche_equivalency_null_models.csv"))
  readr::write_csv(test_results$similarity_null, file.path(tables_dir, "niche_similarity_null_models.csv"))
  readr::write_csv(test_results$errors, file.path(tables_dir, "analysis_errors.csv"))
  readr::write_csv(run_settings, file.path(tables_dir, "run_settings.csv"))

  analysis_object <- list(
    validation = validation,
    pairs = pairs,
    env_space = env_space,
    pca_results = pca_results,
    grids = grids,
    test_results = test_results,
    run_settings = run_settings
  )
  if (exists("pack_cache_value", mode = "function")) {
    analysis_object <- pack_cache_value(analysis_object)
  }

  saveRDS(analysis_object, file.path(objects_dir, "niche_overlap_analysis.rds"))
}

# Save a compact log of implementation changes and assumptions.
write_modification_log <- function(project_dir, run_settings, validation, pairs) {
  path <- file.path(project_dir, "Results", "workflow_modifications.md")
  lines <- c(
    "# Workflow modifications",
    "",
    "- Replaced the previous three-group Oxytrigona workflow with an automated parasite-host workflow based on the `Interaction` column.",
    "- Reads `Data/occurrences_bees_parasite_host.csv` as a semicolon-delimited file and parses coordinate fields robustly.",
    "- Validates missing, non-numeric, out-of-range, and duplicate species-coordinate records before analysis.",
    "- Uses all candidate parasite-host combinations that pass the minimum cleaned occurrence threshold unless `Data/parasite_host_pairs.csv` is supplied.",
    "- Keeps the Broennimann et al. environmental PCA and ecospat density-grid approach, but automates it for all valid pairs.",
    "- Adds parasite-host terminology for stability, parasite-exclusive environmental use, and host environmental space unfilled by the parasite.",
    "- Enables randomization-based niche equivalency and directional niche similarity tests with fixed seeds and ecospat-level parallel workers.",
    "- Revises pair-specific niche dynamics figures to use solid full-available-niche outlines and dashed 50% density contours rather than multiple density contour levels.",
    "- Revises the supplementary PDF to expose reproducibility-relevant analytical code while hiding PDF table styling and other presentational infrastructure.",
    "- Displays reported data objects before formatted tables so readers can identify the summarized objects without seeing table-formatting helper calls.",
    "- Saves figures under `Figures/`, pair-specific niche-space plots under `Figures/Niche_Plots/`, analytical tables and null-model outputs under `Results/tables/`, validation outputs under `Results/validation/`, and R objects under `Results/objects/`.",
    "- Produces publication-oriented occurrence maps, background maps, PCA loading plots, Schoener's D heatmaps, niche dynamics stacked bars, and pair-specific ggplot2 niche dynamics figures.",
    "",
    "## Key settings",
    "",
    paste0("- Minimum cleaned unique occurrences per species: ", run_settings$minimum_occurrences),
    paste0("- MCP buffer size in decimal degrees: ", run_settings$buffer_degrees),
    paste0("- Environmental grid resolution: ", run_settings$grid_resolution),
    paste0("- Test repetitions: ", run_settings$test_repetitions),
    paste0("- Randomization tests run: ", run_settings$run_randomization_tests),
    paste0("- Candidate pairs: ", nrow(pairs$all)),
    paste0("- Valid analyzed pairs: ", nrow(pairs$valid)),
    paste0("- Occurrence issue records written: ", nrow(validation$issues))
  )
  writeLines(lines, path)
  path
}
