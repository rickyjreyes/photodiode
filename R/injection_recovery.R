# ==============================================================================
# injection_recovery.R
# Synthesize waveforms with known structure and measure recovery / power.
# Critically distinguishes (a) detecting a 60 Hz ladder, (b) detecting exact
# 2:3:4 arithmetic, (c) detecting log-cos modulation beyond the ladder.
# ==============================================================================

# Build a synthetic waveform of a chosen type.
pd_synth_waveform <- function(type, n, sr, amp = 1, noise_sd = 0.05,
                              f0 = 60, n_harm = 6, quant = 0,
                              ar_coef = numeric(0)) {
  tt <- seq(0, by = 1/sr, length.out = n)
  noise <- if (length(ar_coef))
    as.numeric(stats::arima.sim(n = n, model = list(ar = ar_coef), sd = noise_sd))
    else stats::rnorm(n, sd = noise_sd)
  sig <- switch(type,
    white       = numeric(n),
    colored     = numeric(n),
    mains_ladder= {
      s <- numeric(n)
      for (h in 1:n_harm) s <- s + (amp / h) * cos(2*pi*h*f0*tt + stats::runif(1,0,2*pi))
      s
    },
    ladder_234  = amp*cos(2*pi*120*tt) + 1.5*amp*cos(2*pi*180*tt) + 2*amp*cos(2*pi*240*tt),
    logcos_env  = {
      s <- numeric(n)
      for (h in 1:n_harm) {
        env <- 1 + 0.5*cos(44.5*log(h*f0))
        s <- s + (amp*env/h) * cos(2*pi*h*f0*tt + stats::runif(1,0,2*pi))
      }
      s
    },
    decaying    = {
      s <- numeric(n); dec <- exp(-tt/0.2)
      for (h in 1:n_harm) s <- s + dec*(amp/h)*cos(2*pi*h*f0*tt)
      s
    },
    numeric(n))
  z <- sig + noise
  if (quant > 0) z <- round(z / quant) * quant
  z - mean(z)
}

# One injection-recovery trial.
pd_injection_trial <- function(type, n, sr, amp, fmin, fmax, kgrid, params,
                               null_thresh) {
  z <- pd_synth_waveform(type, n, sr, amp = amp)
  sp <- pd_fft_canonical(z, sr)
  sel <- sp$frequency_Hz >= fmin & sp$frequency_Hz <= fmax & sp$frequency_Hz > 0
  spc <- sp[sel, ]
  peaks <- pd_detect_peaks(spc, params$peak_height_frac, params$top_peaks)
  ladder_detected <- nrow(peaks) >= 3
  f0_recovered <- if (nrow(peaks)) {
    hm <- pd_harmonic_models(peaks$frequency_Hz)
    abs(hm$f0_Hz[hm$model == "M1_f0free,off=0"] - 60) < 1
  } else FALSE
  ld <- pd_build_log_domain(spc, params$detrend)
  dchi <- pd_scanmax_logcos(ld$ell_ln_f, ld$y_detrended, kgrid)$delta_chi2
  data.frame(type = type, amp = amp,
             ladder_detected = ladder_detected,
             f0_recovered = f0_recovered,
             n_peaks = nrow(peaks),
             logcos_delta_chi2 = dchi,
             logcos_significant = dchi >= null_thresh)
}

pd_injection_recovery <- function(wf, fmin, fmax, kgrid, params,
                                  types = c("white","mains_ladder","ladder_234","logcos_env"),
                                  amps = c(0, 0.25, 0.5, 1, 2, 4),
                                  trials = 20, null_thresh = 21.14) {
  rows <- list(); idx <- 1
  for (ty in types) for (a in amps) {
    res <- lapply(seq_len(trials), function(i)
      pd_injection_trial(ty, wf$n, wf$sample_rate, a, fmin, fmax, kgrid, params, null_thresh))
    df <- do.call(rbind, res)
    rows[[idx]] <- data.frame(
      type = ty, amplitude = a, trials = trials,
      p_ladder_detected = mean(df$ladder_detected),
      p_f0_recovered = mean(df$f0_recovered),
      p_logcos_significant = mean(df$logcos_significant),
      mean_delta_chi2 = mean(df$logcos_delta_chi2))
    idx <- idx + 1
  }
  do.call(rbind, rows)
}
