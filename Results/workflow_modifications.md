# Workflow modifications

- Replaced the previous three-group Oxytrigona workflow with an automated parasite-host workflow based on the `Interaction` column.
- Reads `Data/occurrences_bees_parasite_host.csv` as a semicolon-delimited file and parses coordinate fields robustly.
- Validates missing, non-numeric, out-of-range, and duplicate species-coordinate records before analysis.
- Uses all candidate parasite-host combinations that pass the minimum cleaned occurrence threshold unless `Data/parasite_host_pairs.csv` is supplied.
- Keeps the Broennimann et al. environmental PCA and ecospat density-grid approach, but automates it for all valid pairs.
- Adds parasite-host terminology for stability, parasite-exclusive environmental use, and host environmental space unfilled by the parasite.
- Saves figures under `Figures/`, pair-specific niche-space plots under `Figures/Niche_Plots/`, analytical tables and null-model outputs under `Results/tables/`, validation outputs under `Results/validation/`, and R objects under `Results/objects/`.
- Produces publication-oriented occurrence maps, background maps, PCA loading plots, Schoener's D heatmaps, niche dynamics stacked bars, and pair-specific ggplot2 niche dynamics figures.

## Key settings

- Minimum cleaned unique occurrences per species: 5
- MCP buffer size in decimal degrees: 2
- Environmental grid resolution: 100
- Test repetitions: 0
- Randomization tests run: FALSE
- Candidate pairs: 49
- Valid analyzed pairs: 30
- Occurrence issue records written: 42
