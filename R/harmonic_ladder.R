# ==============================================================================
# harmonic_ladder.R
# Harmonic-ladder statistic and descriptive ladder properties.
# ==============================================================================

# Standardized local power at each declared harmonic of f0 within the spectrum.
pd_harmonic_amplitudes <- function(spec, f0 = 60, n_harmonics = 16,
                                   tol_bins = 1) {
  df <- spec$frequency_Hz[2] - spec$frequency_Hz[1]
  rows <- lapply(seq_len(n_harmonics), function(n) {
    ft <- n * f0
    if (ft < min(spec$frequency_Hz) || ft > max(spec$frequency_Hz))
      return(NULL)
    idx <- which.min(abs(spec$frequency_Hz - ft))
    lo <- max(1, idx - 15); hi <- min(nrow(spec), idx + 15)
    local <- spec$amplitude[lo:hi]
    noise <- stats::median(local); nsd <- stats::mad(local) + 1e-12
    amp <- max(spec$amplitude[max(1, idx-tol_bins):min(nrow(spec), idx+tol_bins)])
    data.frame(harmonic = n, target_Hz = ft,
               bin_frequency_Hz = spec$frequency_Hz[idx],
               amplitude = amp,
               local_noise = noise,
               std_local_power = (amp - noise) / nsd,
               detected = amp > noise + 3 * nsd)
  })
  do.call(rbind, rows)
}

# Ladder statistic: sum of standardized local power over detected harmonics,
# penalized by frequency residual and missing harmonics.
pd_ladder_statistic <- function(spec, f0 = 60, n_harmonics = 16) {
  ha <- pd_harmonic_amplitudes(spec, f0, n_harmonics)
  if (is.null(ha) || !nrow(ha)) return(list(statistic = 0, detail = ha))
  det <- ha[ha$detected, , drop = FALSE]
  n_det <- nrow(det)
  # longest consecutive run of detected harmonics
  hs <- sort(ha$harmonic[ha$detected])
  run <- if (length(hs)) max(rle(c(1, diff(hs) == 1))$lengths[
                              rle(c(1, diff(hs) == 1))$values == 1], 0) + 0 else 0
  longest <- 0; cur <- 0; prev <- NA
  for (h in hs) { if (!is.na(prev) && h == prev + 1) cur <- cur + 1 else cur <- 1
                  longest <- max(longest, cur); prev <- h }
  total_power <- sum(pmax(0, det$std_local_power))
  broadband <- sum(pmax(0, ha$std_local_power)) + 1e-9
  resid_rms <- if (n_det) sqrt(mean((det$bin_frequency_Hz - det$harmonic * f0)^2)) else NA
  # odd/even balance
  odd <- sum(det$amplitude[det$harmonic %% 2 == 1])
  even <- sum(det$amplitude[det$harmonic %% 2 == 0])
  stat <- total_power * (longest / n_harmonics)
  list(
    statistic = stat,
    n_detected = n_det,
    longest_run = longest,
    fundamental = f0,
    freq_residual_rms = resid_rms,
    total_ladder_power = total_power,
    ladder_fraction = total_power / broadband,
    odd_sum = odd, even_sum = even,
    odd_even_balance = if (even > 0) odd / even else NA,
    detail = ha
  )
}
