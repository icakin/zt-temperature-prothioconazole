# =============================================================================
# run_all.R — Master pipeline runner
# =============================================================================
# Sources all scripts in order. Each script sources 00_config.R internally.
#
# Pipeline order:
#   01  Wide -> long oxygen data
#   02  Spline trimming + diagnostics
#  (03  Interactive trim/exclusion selector — MANUAL, not run here)
#   04  Oxygen model fits (nlsLM) -> carbon units + descriptive plots
#   05  Descriptive result figures (review before MCMC)
#   06  Bayesian thermal models (Arrhenius, Sharpe-Schoolfield, hierarchical)
#   07  All plots, effect sizes, synergy analysis, publication figures
#   08  Dose-response Hill + log-linear models (reads fig1_extra CSV from 07)
#  (09  Optional prior-sensitivity refit — run manually)
#
# You can also run any script individually — each one loads its own
# dependencies via source("00_config.R").
#
# NOTE ON 03_trim_selector.R (interactive Shiny app):
#   It is NOT sourced here because it needs manual clicking. To use it, run
#   01 and 02, then open 03_trim_selector.R and "Run App" to set fit windows
#   and choose excluded samples. It writes manual_fit_windows.csv and
#   plot_exclude_points.csv, which 04 picks up automatically. Then run
#   run_all.R (or 04 onward). If you skip 03, run_all.R uses whatever choices
#   are already saved in those CSVs (or none).
# =============================================================================

script_dir <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
                   rstudioapi::isAvailable() &&
                   nzchar(rstudioapi::getActiveDocumentContext()$path)) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}

run_script <- function(name) {
  path <- file.path(script_dir, name)
  message("\n", strrep("=", 70))
  message("Running: ", name)
  message(strrep("=", 70), "\n")
  source(path, local = FALSE)
}

run_script("01_longdata.R")
run_script("02_trimming.R")
run_script("04_oxygen_fits.R")
run_script("05_result_figures.R")   # descriptive figures (review before MCMC)
run_script("06_bayesian_models.R")
run_script("07_plots_and_effects.R")
run_script("08_dose_response_brms.R")

message("\n", strrep("=", 70))
message("Pipeline complete.")
message(strrep("=", 70))
