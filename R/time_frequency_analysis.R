# ==============================================================================
# time_frequency_analysis.R
# STFT / spectrogram and time-resolved harmonic tracking.
# ==============================================================================

# Short-time Fourier transform. Returns a matrix [freq x time] of amplitude.
pd_stft <- function(y, sample_rate, seg_len = NULL, overlap = 0.75,
                    window = "hann") {
  N <- length(y)
  if (is.null(seg_len)) seg_len <- min(N, max(256, 2^floor(log2(N/8))))
  step <- max(1, floor(seg_len * (1 - overlap)))
  starts <- seq(1, N - seg_len + 1, by = step)
  w <- pd_window(seg_len, window)
  half <- seg_len %/% 2 + 1
  freq <- (0:(half - 1)) * sample_rate / seg_len
  tcent <- (starts + seg_len/2 - 1) / sample_rate
  S <- matrix(0, nrow = half, ncol = length(starts))
  for (j in seq_along(starts)) {
    seg <- y[starts[j]:(starts[j] + seg_len - 1)]
    seg <- (seg - mean(seg)) * w
    S[, j] <- Mod(stats::fft(seg))[1:half]
  }
  list(freq = freq, time = tcent, amplitude = S, seg_len = seg_len,
       freq_resolution_Hz = sample_rate / seg_len,
       time_resolution_s = seg_len / sample_rate)
}

# Track amplitude and phase of a target frequency across STFT frames.
pd_track_harmonic <- function(stft, target_hz) {
  idx <- which.min(abs(stft$freq - target_hz))
  data.frame(time_s = stft$time,
             target_Hz = target_hz,
             amplitude = stft$amplitude[idx, ])
}

# Time-local ladder stability: does the harmonic set persist across the capture?
pd_persistence_within_capture <- function(wf, fmin, fmax, harmonics = c(60,120,180,240),
                                          seg_len = NULL) {
  st <- pd_stft(wf$y, wf$sample_rate, seg_len = seg_len)
  tracks <- lapply(harmonics, function(h) {
    tr <- pd_track_harmonic(st, h); tr$harmonic <- h; tr
  })
  df <- do.call(rbind, tracks)
  # instantaneous RMS per frame
  list(stft = st, tracks = df,
       n_frames = length(st$time),
       capture_duration_s = wf$n / wf$sample_rate,
       note = paste("Persistence here means stability WITHIN the",
                    sprintf("%.3f s", wf$n / wf$sample_rate),
                    "capture only. No inference about seconds/minutes/days."))
}
