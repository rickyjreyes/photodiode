# ==============================================================================
# python_r_parity.R
# PARITY_MODE: reproduce the canonical scan_npy_logcos.py calculation in pure R
# and compare deterministic outputs against the committed Python summary.json.
#
# Stochastic outputs (shuffle null) are NOT required to match draw-for-draw;
# they are compared by distribution summaries. R never calls the Python script.
# ==============================================================================

# Canonical deterministic targets from outputs_npy_logcos/summary.json.
PD_PARITY_TARGETS <- list(
  best_k        = 44.53413353338335,
  delta_chi2    = 21.13908290898513,
  chi2_null     = 491.00000000000006,
  chi2_logcos   = 469.8609170910149,
  A             = 0.2926851234475271,
  phi           = 0.8552744569803362,
  active_delta_ell = 3.912023005428146,
  active_winding   = 27.72774419215855
)

# Run the full canonical pipeline in R and write parity outputs.
pd_run_parity <- function(wf, dirs, fmin = 20, fmax = 1000,
                          kmin = 0.5, kmax = 80, nk = 4000,
                          null_n = 1000, seed = 12345) {
  sr <- wf$sample_rate
  spec_all <- pd_fft_canonical(wf$y, sr)
  sel <- spec_all$frequency_Hz >= fmin & spec_all$frequency_Hz <= fmax &
         spec_all$frequency_Hz > 0
  spec <- spec_all[sel, ]
  peaks <- pd_detect_peaks(spec, 0.05, 30)
  ratios <- pd_ratio_pairs(peaks)
  triplets <- pd_koide_triplets(peaks)
  ld <- pd_build_log_domain(spec, "quadratic")
  sc <- pd_scan_logcos(ld$ell_ln_f, ld$y_detrended, kmin, kmax, nk)
  best <- sc$best
  best$active_delta_ell <- pd_active_delta_ell(fmin, fmax)
  best$active_winding_n <- pd_active_winding(best$k, fmin, fmax)

  # shuffle null
  kgrid <- seq(kmin, kmax, length.out = nk)
  set.seed(seed)
  null_best <- numeric(null_n); null_k <- numeric(null_n)
  for (i in seq_len(null_n)) {
    sm <- pd_scanmax_logcos(ld$ell_ln_f, sample(ld$y_detrended), kgrid)
    null_best[i] <- sm$delta_chi2; null_k[i] <- sm$k
  }
  p <- pd_pvalue_geq(best$delta_chi2, null_best)

  # write parity outputs
  utils::write.csv(spec, file.path(dirs$parity, "spectrum.csv"), row.names = FALSE)
  utils::write.csv(peaks, file.path(dirs$parity, "top_peaks.csv"), row.names = FALSE)
  utils::write.csv(ratios, file.path(dirs$parity, "ratio_pairs.csv"), row.names = FALSE)
  utils::write.csv(triplets, file.path(dirs$parity, "koide_triplets.csv"), row.names = FALSE)
  utils::write.csv(sc$scan, file.path(dirs$parity, "logcos_scan.csv"), row.names = FALSE)
  utils::write.csv(data.frame(null_best_delta_chi2 = null_best, null_best_k = null_k),
                   file.path(dirs$parity, "logcos_null.csv"), row.names = FALSE)
  summ <- list(
    input = wf$source, samples = wf$n, sample_rate_Hz = sr,
    duration_s = wf$n / sr, frequency_window_Hz = c(fmin, fmax),
    k_scan = c(kmin, kmax, nk),
    best_logcos = c(best, list(p_shuffle_scanmax = p, null_n = null_n)),
    null_summary = list(mean = mean(null_best), sd = stats::sd(null_best),
                        median = stats::median(null_best),
                        p95 = unname(stats::quantile(null_best, .95)),
                        p99 = unname(stats::quantile(null_best, .99)),
                        max = max(null_best),
                        exceedances = sum(null_best >= best$delta_chi2)))
  pd_write_json(summ, file.path(dirs$parity, "summary.json"))

  list(spec = spec, peaks = peaks, ratios = ratios, triplets = triplets,
       scan = sc, best = best, p = p, null_best = null_best, summary = summ)
}

# Compare R deterministic results against canonical Python targets.
pd_parity_compare <- function(best, tol = 1e-6) {
  got <- list(
    best_k = best$k, delta_chi2 = best$delta_chi2,
    chi2_null = best$chi2_null, chi2_logcos = best$chi2_logcos,
    A = best$A, phi = best$phi,
    active_delta_ell = best$active_delta_ell,
    active_winding = best$active_winding_n
  )
  rows <- lapply(names(PD_PARITY_TARGETS), function(nm) {
    py <- PD_PARITY_TARGETS[[nm]]; r <- got[[nm]]
    reldiff <- abs(r - py) / (abs(py) + 1e-30)
    data.frame(quantity = nm, python = py, r = r,
               abs_diff = abs(r - py), rel_diff = reldiff,
               within_tol = reldiff < tol)
  })
  do.call(rbind, rows)
}
