test_that("pure-R npy reader round-trips float64", {
  tmp <- tempfile(fileext = ".npy")
  x <- rnorm(123)
  pd_test_write_npy(tmp, x)
  y <- pd_read_npy(tmp)
  expect_equal(length(y), 123)
  expect_equal(y, x, tolerance = 1e-12)
})

test_that("checksum is stable", {
  tmp <- tempfile(); writeLines("abc", tmp)
  expect_identical(pd_checksum(tmp), pd_checksum(tmp))
})

test_that("Hann coherent gain ~ 0.5 and ENBW ~ 1.5 bins", {
  w <- pd_window(1024, "hann"); p <- pd_window_properties(w)
  expect_equal(p$coherent_gain, 0.5, tolerance = 1e-3)
  expect_equal(p$enbw_bins, 1.5, tolerance = 2e-2)
})

test_that("single-sided amplitude recovers a known sine amplitude & frequency", {
  sr <- 1000; N <- 4000; tt <- (0:(N-1))/sr
  y <- 2.5 * cos(2*pi*77*tt)
  sp <- pd_amplitude_spectrum(y - mean(y), sr, "hann")
  i <- which.max(sp$amplitude)
  expect_equal(sp$frequency_Hz[i], 77, tolerance = 0.5)
  expect_equal(sp$amplitude[i], 2.5, tolerance = 0.05)
})

test_that("PSD normalization is positive and finite", {
  sr <- 1000; N <- 2048; y <- rnorm(N)
  ps <- pd_psd(y - mean(y), sr)
  expect_true(all(is.finite(ps$psd)))
  expect_true(all(ps$psd >= 0))
})

test_that("R reproduces canonical log-cos targets from the committed waveform", {
  skip_if_not(file.exists(file.path(PD_REPO, "captures/20260309_010802/waveform_raw.npy")))
  wf <- pd_load_waveform(file.path(PD_REPO, "captures/20260309_010802/waveform_raw.npy"), 20000)
  sp <- pd_fft_canonical(wf$y, 20000)
  sel <- sp$frequency_Hz >= 20 & sp$frequency_Hz <= 1000 & sp$frequency_Hz > 0
  ld <- pd_build_log_domain(sp[sel,], "quadratic")
  sc <- pd_scan_logcos(ld$ell_ln_f, ld$y_detrended, 0.5, 80, 4000)
  expect_equal(sc$best$k, 44.53413353338335, tolerance = 1e-9)
  expect_equal(sc$best$delta_chi2, 21.13908290898513, tolerance = 1e-6)
  expect_equal(sc$best$A, 0.2926851234475271, tolerance = 1e-6)
  expect_equal(sc$best$chi2_null, 491, tolerance = 1e-6)
})
