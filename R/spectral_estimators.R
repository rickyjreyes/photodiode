# ==============================================================================
# spectral_estimators.R
# Window functions, their coherent-gain / ENBW properties, and a family of
# spectral estimators (canonical unnormalized Hann FFT, normalized single-sided
# amplitude, PSD in V^2/Hz, Welch PSD, and a multitaper estimate when available).
# ==============================================================================

# ---- Window functions (length-N, symmetric, matching numpy conventions) ------
pd_window <- function(N, type = c("hann", "hamming", "blackmanharris", "flattop",
                                  "rectangular")) {
  type <- match.arg(type)
  n <- 0:(N - 1)
  switch(type,
    rectangular    = rep(1, N),
    # numpy.hanning: symmetric, divides by (N-1)
    hann           = 0.5 - 0.5 * cos(2 * pi * n / (N - 1)),
    hamming        = 0.54 - 0.46 * cos(2 * pi * n / (N - 1)),
    blackmanharris = {
      a <- c(0.35875, 0.48829, 0.14128, 0.01168)
      a[1] - a[2]*cos(2*pi*n/(N-1)) + a[3]*cos(4*pi*n/(N-1)) - a[4]*cos(6*pi*n/(N-1))
    },
    flattop        = {
      a <- c(0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368)
      a[1] - a[2]*cos(2*pi*n/(N-1)) + a[3]*cos(4*pi*n/(N-1)) -
        a[4]*cos(6*pi*n/(N-1)) + a[5]*cos(8*pi*n/(N-1))
    }
  )
}

# Window properties: coherent gain, ENBW (bins), amplitude correction.
pd_window_properties <- function(w) {
  N <- length(w)
  cg <- mean(w)                          # coherent gain
  enbw_bins <- N * sum(w^2) / (sum(w)^2) # equivalent noise bandwidth in bins
  list(
    N = N,
    coherent_gain = cg,
    amplitude_correction = 1 / cg,
    enbw_bins = enbw_bins,
    sum_w = sum(w),
    sum_w2 = sum(w^2)
  )
}

# ---- Canonical unnormalized Hann FFT amplitude (PARITY) ----------------------
# Reproduces compute_fft_spectrum() in scan_npy_logcos.py exactly.
pd_fft_canonical <- function(y, sample_rate) {
  N <- length(y)
  w <- pd_window(N, "hann")
  yw <- y * w
  ft <- stats::fft(yw)
  half <- N %/% 2 + 1
  spec <- Mod(ft)[1:half]
  freq <- (0:(half - 1)) * sample_rate / N
  data.frame(frequency_Hz = freq, amplitude = spec)
}

# ---- Correctly normalized single-sided amplitude spectrum --------------------
# Amplitude of a sinusoid of amplitude A is recovered as ~A.
pd_amplitude_spectrum <- function(y, sample_rate, window = "hann") {
  N <- length(y)
  w <- pd_window(N, window)
  cg <- mean(w)
  yw <- y * w
  ft <- stats::fft(yw)
  half <- N %/% 2 + 1
  mag <- Mod(ft)[1:half]
  amp <- mag / (N * cg) * 2
  amp[1] <- amp[1] / 2                       # DC not doubled
  if (N %% 2 == 0) amp[half] <- amp[half] / 2 # Nyquist not doubled
  freq <- (0:(half - 1)) * sample_rate / N
  data.frame(frequency_Hz = freq, amplitude = amp)
}

# ---- Power spectral density in V^2/Hz ----------------------------------------
pd_psd <- function(y, sample_rate, window = "hann") {
  N <- length(y)
  w <- pd_window(N, window)
  U <- sum(w^2)                              # window power normalization
  yw <- y * w
  ft <- stats::fft(yw)
  half <- N %/% 2 + 1
  psd <- (Mod(ft)[1:half])^2 / (sample_rate * U)
  psd[2:(half - (N %% 2 == 0))] <- psd[2:(half - (N %% 2 == 0))] * 2  # single-sided
  freq <- (0:(half - 1)) * sample_rate / N
  data.frame(frequency_Hz = freq, psd = psd)
}

# ---- Welch PSD ---------------------------------------------------------------
pd_welch <- function(y, sample_rate, nperseg = NULL, overlap = 0.5,
                     window = "hann") {
  N <- length(y)
  if (is.null(nperseg)) nperseg <- max(256, 2^floor(log2(N / 8)))
  nperseg <- min(nperseg, N)
  step <- max(1, floor(nperseg * (1 - overlap)))
  starts <- seq(1, N - nperseg + 1, by = step)
  w <- pd_window(nperseg, window)
  U <- sum(w^2)
  half <- nperseg %/% 2 + 1
  acc <- numeric(half)
  for (s in starts) {
    seg <- y[s:(s + nperseg - 1)]
    seg <- seg - mean(seg)
    ft <- stats::fft(seg * w)
    p <- (Mod(ft)[1:half])^2 / (sample_rate * U)
    p[2:(half - (nperseg %% 2 == 0))] <- p[2:(half - (nperseg %% 2 == 0))] * 2
    acc <- acc + p
  }
  psd <- acc / length(starts)
  freq <- (0:(half - 1)) * sample_rate / nperseg
  list(spectrum = data.frame(frequency_Hz = freq, psd = psd),
       nperseg = nperseg, n_segments = length(starts))
}

# ---- Multitaper PSD (DPSS via 'multitaper' pkg if present; else NA) ----------
pd_multitaper <- function(y, sample_rate, nw = 4, k = 7) {
  if (!requireNamespace("multitaper", quietly = TRUE)) {
    return(list(available = FALSE, spectrum = NULL))
  }
  ts <- stats::ts(y, frequency = sample_rate)
  sp <- multitaper::spec.mtm(ts, nw = nw, k = k, plot = FALSE)
  list(available = TRUE,
       spectrum = data.frame(frequency_Hz = sp$freq, psd = sp$spec))
}

# ---- Per-window summary table (coherent-gain-aware) --------------------------
pd_window_comparison <- function(y, sample_rate, fmin, fmax,
                                 windows = c("hann","hamming","blackmanharris",
                                             "flattop","rectangular")) {
  N <- length(y); df <- sample_rate / N
  rows <- lapply(windows, function(wt) {
    w <- pd_window(N, wt); props <- pd_window_properties(w)
    sp <- pd_amplitude_spectrum(y, sample_rate, wt)
    sel <- sp$frequency_Hz >= fmin & sp$frequency_Hz <= fmax & sp$frequency_Hz > 0
    spc <- sp[sel, ]
    ip <- which.max(spc$amplitude)
    pf <- spc$frequency_Hz[ip]; pa <- spc$amplitude[ip]
    noise <- stats::median(spc$amplitude)
    data.frame(
      window = wt,
      coherent_gain = props$coherent_gain,
      enbw_bins = props$enbw_bins,
      enbw_Hz = props$enbw_bins * df,
      freq_resolution_Hz = df,
      amplitude_correction = props$amplitude_correction,
      peak_frequency_Hz = pf,
      peak_amplitude = pa,
      local_noise_floor = noise,
      snr = pa / noise,
      resolution_limited = df > 1
    )
  })
  do.call(rbind, rows)
}
