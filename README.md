# Parasite-host bee niche overlap

This project estimates environmental niche overlap between cleptoparasitic bee species and host bee species using the framework of Broennimann et al. (2012), as implemented in `ecospat`.

## RStudio Project use

Open `bee_parasites.Rproj` in RStudio before running analyses. The workflow assumes that the working directory is the project root defined by the `.Rproj` file, so all paths are relative to that root.

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

After opening the `.Rproj`, render the supplementary PDF with randomization tests enabled:

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
    run_randomization_tests = TRUE,
    test_repetitions = 100,
    future_workers = 4,
    future_strategy = "multisession"
  ),
  use_cache = TRUE
)
```

The `_targets.R` file provides a reproducible pipeline scaffold using the same randomization defaults.

## Pair definitions

If `Data/parasite_host_pairs.csv` exists, the analysis is restricted to the curated parasite-host pairs in that file. The file must contain at least `parasite` and `host` columns. If the file is absent, all parasite-host combinations passing occurrence and environmental filters are analyzed and flagged as exploratory in the report.

## Randomization tests

Niche equivalency and similarity tests are enabled by default with 100 repetitions, `random_seed = 42`, and four workers passed to each `ecospat` randomization test. The workflow runs one niche equivalency test and two directional niche similarity tests for each valid parasite-host pair. Pair-level `future` parallelization is used only for non-randomization runs because serialized `ecospat` grid objects can fail in separate R sessions during randomization.

The 100-repetition setting provides a reproducible inferential run and p-value resolution of approximately 0.01. For final publication inference, consider increasing `test_repetitions` to 999 or more if computing time permits, then rerun the workflow from the RStudio Project root.

Randomization outputs are written to `Results/tables/niche_equivalency_null_models.csv`, `Results/tables/niche_similarity_null_models.csv`, and `Results/tables/niche_overlap_metrics.csv`. These generated result files are ignored by git to keep the repository lightweight; they can be regenerated from the tracked code and input data.

## Niche dynamics figures

The all-pairs niche dynamics metric figure displays stability, parasite-exclusive environment, and host-unfilled environment as separate 0-1 metrics; the values are not normalized to make every pair sum to 100%. Pair-specific niche dynamics figures use the host niche as the reference and the parasite niche as the focal niche. Filled cells show dynamic occupancy categories only; unclassified available background cells are left blank. Solid contours mark species-specific background availability and dashed contours mark the 50% occurrence-density contour for each species. `Figures/Figure5_niche_dynamics_all_pairs.pdf` is a single faceted all-pairs niche-space figure with shared legends, and `Figures/Niche_Plots/` contains the individual pair plots. Multiple density contour levels are intentionally not drawn.

## GitHub remote

The intended remote repository is:

```sh
https://github.com/BrunoVilela/parasite_host_bees.git
```

No credentials are stored in the project files.
