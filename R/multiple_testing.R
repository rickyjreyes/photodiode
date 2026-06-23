# ==============================================================================
# multiple_testing.R
# Declared analysis registry + multiplicity corrections. The canonical
# scan-global p (look-elsewhere over k only) is NOT a correction for capture,
# window, detrending, peak-threshold, ratio, or control selection.
# ==============================================================================

pd_multiple_testing <- function(pvalues, labels = NULL,
                                 scan_global_p = NA, spec_grid_size = NA) {
  if (is.null(labels)) labels <- paste0("test_", seq_along(pvalues))
  m <- length(pvalues)
  holm <- stats::p.adjust(pvalues, method = "holm")
  bh   <- stats::p.adjust(pvalues, method = "BH")
  bonf <- stats::p.adjust(pvalues, method = "bonferroni")
  data.frame(
    test = labels,
    raw_p = pvalues,
    holm_p = holm,
    bh_fdr_p = bh,
    bonferroni_p = bonf
  )
}

# Family-wise maximum-statistic correction: repeats the full selection
# (here: max delta_chi2 over the whole specification grid) per simulated draw.
pd_familywise_maxstat <- function(observed_max, null_max_draws) {
  (1 + sum(null_max_draws >= observed_max)) / (1 + length(null_max_draws))
}

pd_analysis_registry <- function() {
  data.frame(
    family = c("captures","frequency_windows","spectral_estimators",
               "peak_thresholds","k_ranges","detrending_models",
               "harmonic_hypotheses","ratio_hypotheses","controls",
               "persistence_tests","ringdown_models"),
    description = c(
      "raw waveforms analyzed",
      "10-500,20-500,20-1000,40-1000 Hz",
      "Hann FFT, Welch, multitaper",
      "peak-height fraction of max",
      "0.5-40, 0.5-80, 0.5-120",
      "linear/quadratic/cubic/spline",
      "f0=60 fixed vs free, M0..M3",
      "2/3, 1/2, 3/4 ratio searches",
      "baseline/sham/disconnected/dark/terminated/channel-swap/battery",
      "within- and across-capture persistence",
      "M0..M5 decay models"),
    stringsAsFactors = FALSE
  )
}
