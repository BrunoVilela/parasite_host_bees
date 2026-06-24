# Parasite-host bee niche overlap

This project estimates environmental niche overlap between cleptoparasitic bee species and host bee species using the framework of Broennimann et al. (2012), as implemented in `ecospat`.

## RStudio Project use

Open `Capitulo 2.Rproj` in RStudio before running analyses. The workflow assumes that the working directory is the project root defined by the `.Rproj` file, so all paths are relative to that root.

## Repository layout

- `Data/`: occurrence data and optional curated parasite-host pair definitions.
- `Rmarkdown/`: supplementary report source.
- `Scripts/functions/`: reusable validation, analysis, plotting, and cache helpers.
- `Scripts/analysis/`: executable workflow wrapper.
- `Scripts/visualization/`: reserved for future standalone plotting scripts.
- `Figures/`: publication figures.
- `Figures/Niche_Plots/`: pair-specific niche-space plots in PDF and PNG.
- `Results/tables/`: analytical tables, overlap metrics, PCA outputs, and test summaries.
- `Results/validation/`: occurrence validation reports and correction logs.
- `Results/objects/`: serialized analysis objects.
- `Results/cache/`: local workflow cache ignored by git.

Generated figures, rendered PDFs, caches, and serialized intermediate objects are intentionally ignored by git so the repository remains lightweight and reproducible.

## Dependencies

Install the required R packages before rendering the supplement. The RMarkdown document checks for required packages and stops if any are missing; it does not install packages automatically during rendering.

## Main workflow

After opening the `.Rproj`, render the supplementary PDF without running randomization tests:

```r
rmarkdown::render(
  "Rmarkdown/parasite_host_bees.Rmd",
  output_dir = "Results",
  output_file = "parasite_host_niche_overlap_supplement.pdf",
  envir = new.env(parent = globalenv())
)
```

Run the modular workflow directly from the project root:

```r
source("Scripts/analysis/run_workflow.R")

results <- run_niche_overlap_workflow(
  project_dir = getwd(),
  settings = list(
    run_randomization_tests = FALSE,
    test_repetitions = 0,
    future_workers = 1
  ),
  use_cache = TRUE
)
```

The `_targets.R` file provides a reproducible pipeline scaffold. It uses the same defaults and keeps randomization disabled until explicitly changed.

## Pair definitions

If `Data/parasite_host_pairs.csv` exists, the analysis is restricted to the curated parasite-host pairs in that file. The file must contain at least `parasite` and `host` columns. If the file is absent, all parasite-host combinations passing occurrence and environmental filters are analyzed and flagged as exploratory in the report.

## Randomization tests

Niche equivalency and similarity tests are implemented but disabled by default because they are computationally expensive. For a final inferential run, set `run_randomization_tests = TRUE`, choose a defensible number of repetitions, and configure `future_workers` and `future_strategy` according to available hardware. Run these analyses from the RStudio Project root so relative paths and caches resolve correctly.

## GitHub remote

The intended remote repository is:

```sh
https://github.com/BrunoVilela/parasite_host_bees.git
```

No credentials are stored in the project files.
