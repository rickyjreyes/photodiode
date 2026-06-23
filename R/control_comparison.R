# ==============================================================================
# control_comparison.R
# ONE standardized analysis function applied unchanged to every capture, plus
# the control matrix and baseline/post harmonic effect sizes.
#
# The same thresholds are used for primary and control runs (no post-hoc
# threshold changes).
# ==============================================================================

# Standardized per-capture analysis. Returns a list with a one-row `summary`
# and the underlying objects. Identical parameters for all captures.
pd_analyze_capture <- function(wf, params = pd_default_params()) {
  sr <- wf$sample_rate
  spec <- pd_fft_canonical(wf$y, sr)
  sel <- spec$frequency_Hz >= params$fmin & spec$frequency_Hz <= params$fmax &
         spec$frequency_Hz > 0
  spc <- spec[sel, ]
  peaks <- pd_detect_peaks(spc, params$peak_height_frac, params$top_peaks)
  ladder <- pd_ladder_statistic(spc, f0 = 60, n_harmonics = params$n_harmonics)
  triplets <- if (nrow(peaks) >= 3) pd_koide_triplets(peaks) else data.frame()
  q_exact <- if (nrow(triplets)) sum(triplets$koide_error <= 1e-6) else 0
  min_koide <- if (nrow(triplets)) min(triplets$koide_error) else NA_real_
  ld <- pd_build_log_domain(spc, params$detrend)
  sc <- pd_scan_logcos(ld$ell_ln_f, ld$y_detrended, params$kmin, params$kmax, params$nk)
  amp60 <- {
    i60 <- which.min(abs(spc$frequency_Hz - 60)); spc$amplitude[i60]
  }
  qc <- pd_waveform_qc(wf)
  summary <- data.frame(
    rms_voltage = qc$rms,
    amp_60Hz = amp60,
    total_harmonic_power = ladder$total_ladder_power,
    n_harmonics = ladder$n_detected,
    ladder_statistic = ladder$statistic,
    q_triplet_count = q_exact,
    min_koide_error = min_koide,
    logcos_best_k = sc$best$k,
    logcos_delta_chi2 = sc$best$delta_chi2,
    logcos_A = sc$best$A
  )
  list(summary = summary, spectrum = spc, peaks = peaks, ladder = ladder,
       triplets = triplets, logdomain = ld, scan = sc, qc = qc, full_spectrum = spec)
}

pd_default_params <- function() {
  list(fmin = 20, fmax = 1000, peak_height_frac = 0.05, top_peaks = 30,
       n_harmonics = 16, detrend = "quadratic", kmin = 0.5, kmax = 80, nk = 4000)
}

# Baseline-vs-post harmonic effect sizes from peak amplitude tables.
# NOTE: amplitudes must be on a comparable normalization; mismatches are flagged.
pd_harmonic_effect_sizes <- function(baseline_peaks, post_peaks,
                                     harmonics = c(1,2,3,4,5,6,8,10),
                                     f0 = 60, comparable = TRUE) {
  rows <- lapply(harmonics, function(h) {
    ft <- h * f0
    bi <- which.min(abs(baseline_peaks$frequency_hz - ft))
    pi_ <- which.min(abs(post_peaks$frequency_Hz - ft))
    ba <- if (length(bi) && abs(baseline_peaks$frequency_hz[bi] - ft) < 5)
            baseline_peaks$amplitude[bi] else NA_real_
    pa <- if (length(pi_) && abs(post_peaks$frequency_Hz[pi_] - ft) < 5)
            post_peaks$amplitude[pi_] else NA_real_
    data.frame(harmonic = h, target_Hz = ft,
               baseline_amplitude = ba, post_amplitude = pa,
               absolute_change = pa - ba,
               fold_change = if (!is.na(ba) && ba != 0) pa / ba else NA_real_)
  })
  out <- do.call(rbind, rows)
  attr(out, "comparable_normalization") <- comparable
  out
}

# Assemble the control matrix from a named list of per-capture analyses + meta.
pd_control_matrix <- function(analyses, meta) {
  rows <- lapply(names(analyses), function(id) {
    s <- analyses[[id]]$summary
    m <- meta[[id]]
    data.frame(analysis_id = id,
               condition = m$phase %||% NA,
               control_type = m$control_type %||% NA,
               elapsed_time_since_excitation_s = m$elapsed_s %||% NA,
               s, check.names = FALSE)
  })
  do.call(rbind, rows)
}
