# ==============================================================================
# null_calibration.R
# Verify that the log-cos scan-max + shuffle-null decision rule has its intended
# false-positive rate. Generate complete null waveforms, run the identical
# pipeline, and measure the empirical alpha at nominal levels.
# ==============================================================================

# One calibration trial: simulate colored-noise waveform, run full pipeline,
# build a small shuffle null, return the shuffle p-value.
pd_calibration_trial <- function(n, sr, fmin, fmax, ar_coef, ar_sd,
                                 kgrid, inner_B, detrend = "quadratic") {
  sim <- if (length(ar_coef))
    as.numeric(stats::arima.sim(n = n, model = list(ar = ar_coef), sd = ar_sd))
    else stats::rnorm(n, sd = ar_sd)
  sim <- sim - mean(sim)
  sp <- pd_fft_canonical(sim, sr)
  sel <- sp$frequency_Hz >= fmin & sp$frequency_Hz <= fmax & sp$frequency_Hz > 0
  ld <- pd_build_log_domain(sp[sel, ], detrend)
  ell <- ld$ell_ln_f; y <- ld$y_detrended
  obs <- pd_scanmax_logcos(ell, y, kgrid)$delta_chi2
  null_best <- numeric(inner_B)
  for (b in seq_len(inner_B)) null_best[b] <- pd_scanmax_logcos(ell, sample(y), kgrid)$delta_chi2
  (1 + sum(null_best >= obs)) / (1 + inner_B)
}

pd_null_calibration <- function(wf, fmin, fmax, kgrid, outer_n = 1000,
                                inner_B = 200, alphas = c(0.10, 0.05, 0.01, 0.001),
                                detrend = "quadratic") {
  ar_fit <- try(stats::ar(wf$y, order.max = 8, aic = TRUE), silent = TRUE)
  ar_coef <- if (!inherits(ar_fit, "try-error")) ar_fit$ar else numeric(0)
  ar_sd <- if (!inherits(ar_fit, "try-error")) sqrt(ar_fit$var.pred) else stats::sd(wf$y)
  ps <- numeric(outer_n)
  for (i in seq_len(outer_n))
    ps[i] <- pd_calibration_trial(wf$n, wf$sample_rate, fmin, fmax,
                                  ar_coef, ar_sd, kgrid, inner_B, detrend)
  rows <- lapply(alphas, function(a) {
    obs <- mean(ps <= a); k <- sum(ps <= a)
    ci <- stats::binom.test(k, outer_n, a)$conf.int
    data.frame(nominal_alpha = a, observed_fpr = obs,
               ci_low = ci[1], ci_high = ci[2],
               expected_count = a * outer_n, observed_count = k,
               calibration_ratio = obs / a, n_outer = outer_n, inner_B = inner_B)
  })
  do.call(rbind, rows)
}
