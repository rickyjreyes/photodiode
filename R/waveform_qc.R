# ==============================================================================
# waveform_qc.R
# Time-domain quality control, calibration audit, and stationarity diagnostics.
# ==============================================================================

pd_skewness <- function(x) { x <- x - mean(x); mean(x^3) / (mean(x^2)^1.5) }
pd_kurtosis <- function(x) { x <- x - mean(x); mean(x^4) / (mean(x^2)^2) }

# Estimate quantization step from the smallest positive difference of sorted
# unique values (robust to a uniformly sampled ADC grid).
pd_quantization_step <- function(y, raw = NULL) {
  u <- sort(unique(round(y, 12)))
  d <- diff(u)
  d <- d[d > 0]
  if (!length(d)) return(NA_real_)
  stats::median(d)
}

pd_clip_fraction <- function(y) {
  rng <- range(y)
  ext <- sum(y == rng[1]) + sum(y == rng[2])
  ext / length(y)
}

# Full QC record for one DC-removed waveform `y` (and optional raw `y_raw`).
pd_waveform_qc <- function(wf) {
  y <- wf$y; yr <- wf$y_raw; t <- wf$t; n <- length(y)
  # linear trend / drift
  lt <- stats::lm(y ~ t); slope <- stats::coef(lt)[2]
  # change point: largest jump in cumulative mean halves
  half <- floor(n/2)
  drift <- abs(mean(yr[1:half]) - mean(yr[(half+1):n]))
  # autocorrelation lag-1 & effective sample size
  ac <- stats::acf(y, lag.max = 50, plot = FALSE)$acf[-1]
  ess <- n / (1 + 2 * sum(ac[ac > 0]))
  # Ljung-Box
  lb <- suppressWarnings(stats::Box.test(y, lag = 20, type = "Ljung-Box"))
  data.frame(
    n = n,
    mean = mean(yr),
    median = stats::median(yr),
    rms = sqrt(mean(yr^2)),
    sd = stats::sd(yr),
    peak_to_peak = diff(range(yr)),
    min = min(yr), max = max(yr),
    skewness = pd_skewness(y),
    kurtosis = pd_kurtosis(y),
    clip_fraction = pd_clip_fraction(yr),
    quantization_step = pd_quantization_step(yr),
    linear_trend_slope = unname(slope),
    half_split_drift = drift,
    stationary_mean = drift < 0.05 * (diff(range(yr)) + 1e-12),
    acf_lag1 = ac[1],
    effective_sample_size = ess,
    ljung_box_p = unname(lb$p.value)
  )
}

# Calibration audit against the documented Siglent convention.
# voltage = raw_byte * vscale / 25 - offset  ->  quantization step = vscale / 25.
pd_calibration_audit <- function(wf, vscale = NA, offset = NA) {
  q <- pd_quantization_step(wf$y_raw)
  expected_q <- if (is.finite(vscale)) vscale / 25 else NA_real_
  list(
    measured_quantization_step = q,
    expected_step_vscale_div25 = expected_q,
    step_ratio = if (is.finite(expected_q) && expected_q > 0) q / expected_q else NA_real_,
    consistent_with_8bit = if (is.finite(expected_q))
      isTRUE(abs(q / expected_q - 1) < 0.05) else NA,
    byte_encoding = "stored as <f8 voltage in .npy (already calibrated)",
    formula = "voltage = raw_byte * vscale / 25 - offset",
    vscale = vscale, offset = offset,
    note = paste("The committed .npy contains float voltages, not raw bytes;",
                 "the calibration formula is verified indirectly via the",
                 "quantization step vscale/25.")
  )
}
