# ==============================================================================
# peak_estimation.R
# Peak detection (matching scipy.find_peaks height rule), sub-bin frequency
# estimation, and harmonic-ladder model comparison M0-M3.
# ==============================================================================

# Strict-interior local maxima with height >= frac * max(amp); top_n by amplitude.
pd_detect_peaks <- function(spec, peak_height_frac = 0.05, top_n = 30) {
  f <- spec$frequency_Hz; a <- spec$amplitude
  n <- length(a)
  if (n < 3) return(data.frame())
  h <- max(a) * peak_height_frac
  idx <- which(a[2:(n-1)] > a[1:(n-2)] & a[2:(n-1)] > a[3:n]) + 1
  idx <- idx[a[idx] >= h]
  if (!length(idx)) return(data.frame())
  d <- data.frame(idx = idx, frequency_Hz = f[idx], amplitude = a[idx])
  d <- d[order(-d$amplitude), ]
  d <- utils::head(d, top_n)
  d$rank_by_amplitude <- seq_len(nrow(d))
  d <- d[order(d$frequency_Hz), ]
  d$rank_by_frequency <- seq_len(nrow(d))
  rownames(d) <- NULL
  d
}

# Quadratic (parabolic) sub-bin interpolation around a spectral bin.
pd_subbin_quadratic <- function(spec, idx) {
  f <- spec$frequency_Hz; a <- spec$amplitude
  if (idx <= 1 || idx >= length(a)) return(list(frequency = f[idx], amplitude = a[idx]))
  ym1 <- a[idx-1]; y0 <- a[idx]; yp1 <- a[idx+1]
  denom <- (ym1 - 2*y0 + yp1)
  delta <- if (denom != 0) 0.5 * (ym1 - yp1) / denom else 0
  df <- f[2] - f[1]
  list(frequency = f[idx] + delta * df,
       amplitude = y0 - 0.25 * (ym1 - yp1) * delta)
}

# Time-domain single-sinusoid NLS refinement near a seed frequency (CI included).
pd_sine_nls <- function(t, y, f_seed) {
  fit <- try(stats::nls(
    y ~ A * cos(2*pi*f*t + phi),
    start = list(A = stats::sd(y) * sqrt(2), f = f_seed, phi = 0),
    control = stats::nls.control(maxiter = 200, warnOnly = TRUE)
  ), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  co <- summary(fit)$coefficients
  fhat <- co["f", "Estimate"]; fse <- co["f", "Std. Error"]
  list(frequency = fhat, freq_se = fse,
       ci_low = fhat - 1.96 * fse, ci_high = fhat + 1.96 * fse,
       amplitude = abs(co["A", "Estimate"]))
}

# Build a peak-estimate table with sub-bin frequencies and 60 Hz residuals.
pd_peak_estimates <- function(spec, peaks, t = NULL, y = NULL,
                              mains_candidates = c(50, 59.94, 60)) {
  if (!nrow(peaks)) return(data.frame())
  df <- spec$frequency_Hz[2] - spec$frequency_Hz[1]
  out <- lapply(seq_len(nrow(peaks)), function(i) {
    idx <- peaks$idx[i]
    sb <- pd_subbin_quadratic(spec, idx)
    nls_est <- NULL
    if (!is.null(t) && !is.null(y)) nls_est <- pd_sine_nls(t, y, sb$frequency)
    fsb <- sb$frequency
    nearest_h <- round(fsb / 60)
    # local SNR
    lo <- max(1, idx - 10); hi <- min(nrow(spec), idx + 10)
    noise <- stats::median(spec$amplitude[lo:hi])
    data.frame(
      bin_frequency_Hz = peaks$frequency_Hz[i],
      subbin_frequency_Hz = fsb,
      nls_frequency_Hz = if (!is.null(nls_est)) nls_est$frequency else NA_real_,
      nls_freq_ci_low  = if (!is.null(nls_est)) nls_est$ci_low else NA_real_,
      nls_freq_ci_high = if (!is.null(nls_est)) nls_est$ci_high else NA_real_,
      amplitude = sb$amplitude,
      local_snr = sb$amplitude / noise,
      peak_width_Hz = df,
      nearest_60Hz_harmonic = nearest_h,
      residual_from_n60 = fsb - nearest_h * 60,
      ci_includes_60_multiple = if (!is.null(nls_est))
        (nls_est$ci_low <= nearest_h*60 & nls_est$ci_high >= nearest_h*60) else NA
    )
  })
  do.call(rbind, out)
}

# Harmonic-ladder model comparison: f_n = f_offset + n * f_0.
# M0: f0=60, offset=0 ; M1: f0 free, offset=0 ; M2: f0,offset free ; M3: free peaks.
pd_harmonic_models <- function(peak_freqs, harmonic_numbers = NULL) {
  f <- sort(peak_freqs)
  if (is.null(harmonic_numbers)) harmonic_numbers <- round(f / 60)
  n <- harmonic_numbers
  rss_aic_bic <- function(rss, npar, m) {
    k <- npar; ll <- -m/2 * (log(2*pi*rss/m) + 1)
    list(rss = rss, aic = -2*ll + 2*k, bic = -2*ll + k*log(m))
  }
  m <- length(f)
  # M0
  r0 <- f - n * 60; s0 <- rss_aic_bic(sum(r0^2), 1, m)
  # M1: minimize over f0 -> f0 = sum(n f)/sum(n^2)
  f0_1 <- sum(n * f) / sum(n^2); r1 <- f - n * f0_1; s1 <- rss_aic_bic(sum(r1^2), 2, m)
  # M2: lm f ~ n
  fit2 <- stats::lm(f ~ n); r2 <- stats::resid(fit2); s2 <- rss_aic_bic(sum(r2^2), 3, m)
  off2 <- stats::coef(fit2)[1]; f0_2 <- stats::coef(fit2)[2]
  # M3: free -> rss 0
  s3 <- list(rss = 0, aic = 2*m, bic = m*log(m))
  data.frame(
    model = c("M0_f0=60,off=0","M1_f0free,off=0","M2_f0,off free","M3_free"),
    f0_Hz = c(60, f0_1, unname(f0_2), NA),
    offset_Hz = c(0, 0, unname(off2), NA),
    rss = c(s0$rss, s1$rss, s2$rss, s3$rss),
    rms_resid = c(sqrt(s0$rss/m), sqrt(s1$rss/m), sqrt(s2$rss/m), 0),
    aic = c(s0$aic, s1$aic, s2$aic, s3$aic),
    bic = c(s0$bic, s1$bic, s2$bic, s3$bic)
  )
}
