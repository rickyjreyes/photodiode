# ==============================================================================
# window_sensitivity.R
# Prespecified specification grid over frequency window, k range, detrending,
# and spectral estimator. Produces a specification curve.
# ==============================================================================

pd_spec_grid <- function() {
  expand.grid(
    fwin = c("10-500","20-500","20-1000","40-1000"),
    krange = c("0.5-40","0.5-80","0.5-120"),
    detrend = c("linear","quadratic","cubic"),
    estimator = c("hann","welch"),
    stringsAsFactors = FALSE
  )
}

pd_parse_range <- function(s) as.numeric(strsplit(s, "-")[[1]])

# Run the log-cos scan for one specification on a raw waveform.
pd_run_spec <- function(wf, fwin, krange, detrend, estimator, nk = 1000) {
  fr <- pd_parse_range(fwin); kr <- pd_parse_range(krange)
  sr <- wf$sample_rate
  if (estimator == "welch") {
    w <- pd_welch(wf$y, sr); sp <- data.frame(frequency_Hz = w$spectrum$frequency_Hz,
                                              amplitude = sqrt(pmax(0, w$spectrum$psd)))
  } else {
    sp <- pd_fft_canonical(wf$y, sr)
  }
  sel <- sp$frequency_Hz >= fr[1] & sp$frequency_Hz <= fr[2] & sp$frequency_Hz > 0
  spc <- sp[sel, ]
  if (nrow(spc) < 8) return(NULL)
  ld <- pd_build_log_domain(spc, detrend)
  sc <- pd_scan_logcos(ld$ell_ln_f, ld$y_detrended, kr[1], kr[2], nk)
  data.frame(fwin = fwin, krange = krange, detrend = detrend, estimator = estimator,
             best_k = sc$best$k, delta_chi2 = sc$best$delta_chi2, A = sc$best$A,
             at_k_boundary = abs(sc$best$k - kr[1]) < 1e-6 | abs(sc$best$k - kr[2]) < 1e-6,
             in_peak_region = sc$best$k > 30 & sc$best$k < 60)
}

pd_window_sensitivity <- function(wf, grid = pd_spec_grid(), nk = 1000) {
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    pd_run_spec(wf, g$fwin, g$krange, g$detrend, g$estimator, nk)
  })
  out <- do.call(rbind, rows)
  out[order(-out$delta_chi2), ]
}
