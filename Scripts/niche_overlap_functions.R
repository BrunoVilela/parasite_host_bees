# Compatibility wrapper. The reusable functions now live under Scripts/functions/.
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (is.na(this_file)) {
  this_file <- file.path(getwd(), "Scripts", "niche_overlap_functions.R")
}
source(file.path(dirname(this_file), "functions", "niche_overlap_functions.R"))
