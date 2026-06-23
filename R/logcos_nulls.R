# ==============================================================================
# logcos_nulls.R
# Canonical shuffle null (PARITY) plus stronger audit nulls:
#   A. block-permutation        D. colored-noise time-domain (AR)
#   B. stationary bootstrap     E. mains-harmonic null
#   C. smooth-envelope residual F. baseline/control empirical null
#
# Every null repeats the FULL log-cos scan-max pipeline and reports the
# (r+1)/(B+1) p-value (never 0).
# ==============================================================================

# Generic null runner: `gen()` returns a detrended/standardized residual vector
# `y` aligned to the fixed `ell`; we scan-max each and compare to observed.
pd_run_null <- function(ell, observed_delta, gen, B, kgrid, label) {
  null_best <- numeric(B); null_k <- numeric(B)
  for (i in seq_len(B)) {
    y <- gen()
    sm <- pd_scanmax_logcos(ell, y, kgrid)
    null_best[i] <- sm$delta_chi2; null_k[i] <- sm$k
  }
  p <- (1 + sum(null_best >= observed_delta)) / (1 + B)
  list(
    label = label, B = B,
    exceedances = sum(null_best >= observed_delta),
    p_value = p,
    resolution_floor = 1 / (B + 1),
    null_mean = mean(null_best), null_sd = stats::sd(null_best),
    null_median = stats::median(null_best),
    p95 = stats::quantile(null_best, 0.95, names = FALSE),
    p99 = stats::quantile(null_best, 0.99, names = FALSE),
    null_max = max(null_best),
    null_best = null_best, null_k = null_k
  )
}

# A. Canonical shuffle null (permute standardized residual).
pd_null_shuffle <- function(ell, y, observed_delta, B, kgrid) {
  pd_run_null(ell, observed_delta, function() sample(y), B, kgrid, "shuffle_canonical")
}

# B. Block-permutation null (permute contiguous blocks of length L).
pd_null_block <- function(ell, y, observed_delta, B, kgrid, block_len = 8) {
  n <- length(y)
  pd_run_null(ell, observed_delta, function() {
    nb <- ceiling(n / block_len)
    blocks <- split(seq_len(n), ceiling(seq_len(n) / block_len))
    perm <- sample(seq_along(blocks))
    unlist(blocks[perm])[seq_len(n)] -> idx
    y[idx]
  }, B, kgrid, paste0("block_perm_L", block_len))
}

# C. Stationary bootstrap (Politis-Romano), preserves local autocorrelation.
pd_null_stationary_boot <- function(ell, y, observed_delta, B, kgrid, p_geom = 1/8) {
  n <- length(y)
  pd_run_null(ell, observed_delta, function() {
    out <- numeric(n); i <- 1
    cur <- sample.int(n, 1)
    while (i <= n) {
      out[i] <- y[cur]; i <- i + 1
      if (stats::runif(1) < p_geom) cur <- sample.int(n, 1)
      else cur <- if (cur == n) 1 else cur + 1
    }
    # restandardize
    (out - mean(out)) / (sqrt(mean((out - mean(out))^2)) + 1e-12)
  }, B, kgrid, "stationary_bootstrap")
}

# D. Colored-noise time-domain null: fit AR to the waveform, simulate, rerun the
#    full FFT + detrend pipeline. Requires the raw waveform + spectral fns.
pd_null_colored_ar <- function(wf, fmin, fmax, observed_delta, B, kgrid,
                               detrend = "quadratic", order_max = 12) {
  ar_fit <- stats::ar(wf$y, order.max = order_max, aic = TRUE)
  sr <- wf$sample_rate; n <- wf$n
  # The FFT frequency grid (hence ell) is identical for every simulation of the
  # same length and sample rate, so compute it once.
  sp0 <- pd_fft_canonical(wf$y, sr)
  sel0 <- sp0$frequency_Hz >= fmin & sp0$frequency_Hz <= fmax & sp0$frequency_Hz > 0
  ell <- log(sp0$frequency_Hz[sel0])
  gen <- function() {
    sim <- as.numeric(stats::arima.sim(
      n = n, model = list(ar = ar_fit$ar),
      sd = sqrt(ar_fit$var.pred)))
    sp <- pd_fft_canonical(sim, sr)
    sel <- sp$frequency_Hz >= fmin & sp$frequency_Hz <= fmax & sp$frequency_Hz > 0
    pd_build_log_domain(sp[sel, ], detrend)$y_detrended
  }
  pd_run_null(ell, observed_delta, gen, B, kgrid, "colored_noise_AR")
}

# E. Mains-harmonic null: synthesize a 60 Hz ladder (fitted amplitudes, random
#    phases) + colored broadband noise + quantization, rerun full pipeline.
pd_null_mains <- function(wf, harmonic_amps, fmin, fmax, observed_delta, B, kgrid,
                          f0 = 60, detrend = "quadratic", quant = 0.12) {
  sr <- wf$sample_rate; n <- wf$n
  tt <- seq(0, by = 1/sr, length.out = n)
  noise_sd <- stats::sd(wf$y)
  ar_fit <- try(stats::ar(wf$y, order.max = 8, aic = TRUE), silent = TRUE)
  null_best <- numeric(B)
  for (i in seq_len(B)) {
    sig <- numeric(n)
    for (h in seq_along(harmonic_amps)) {
      ph <- stats::runif(1, 0, 2*pi)
      sig <- sig + harmonic_amps[h] * cos(2*pi*h*f0*tt + ph)
    }
    noise <- if (!inherits(ar_fit, "try-error") && length(ar_fit$ar))
      as.numeric(stats::arima.sim(n = n, model = list(ar = ar_fit$ar),
                                  sd = sqrt(ar_fit$var.pred)))
      else stats::rnorm(n, sd = noise_sd)
    z <- sig + noise
    if (quant > 0) z <- round(z / quant) * quant
    z <- z - mean(z)
    sp <- pd_fft_canonical(z, sr)
    sel <- sp$frequency_Hz >= fmin & sp$frequency_Hz <= fmax & sp$frequency_Hz > 0
    ld <- pd_build_log_domain(sp[sel, ], detrend)
    null_best[i] <- pd_scanmax_logcos(ld$ell_ln_f, ld$y_detrended, kgrid)$delta_chi2
  }
  p <- (1 + sum(null_best >= observed_delta)) / (1 + B)
  list(label = "mains_harmonic", B = B,
       exceedances = sum(null_best >= observed_delta), p_value = p,
       resolution_floor = 1/(B+1), null_mean = mean(null_best),
       null_sd = stats::sd(null_best), null_median = stats::median(null_best),
       p95 = stats::quantile(null_best, .95, names = FALSE),
       p99 = stats::quantile(null_best, .99, names = FALSE),
       null_max = max(null_best), null_best = null_best,
       null_k = rep(NA_real_, B))
}

# Helper: tidy a list of null results into one comparison data.frame.
pd_null_comparison_table <- function(null_list, observed_delta) {
  rows <- lapply(null_list, function(r) {
    if (is.null(r)) return(NULL)
    data.frame(null_model = r$label, B = r$B,
               observed_delta_chi2 = observed_delta,
               exceedances = r$exceedances, p_value = r$p_value,
               resolution_floor = r$resolution_floor,
               null_mean = r$null_mean, null_sd = r$null_sd,
               null_median = r$null_median, p95 = r$p95, p99 = r$p99,
               null_max = r$null_max)
  })
  do.call(rbind, rows)
}
