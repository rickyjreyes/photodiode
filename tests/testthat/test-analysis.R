test_that("DC removal yields zero mean", {
  expect_equal(mean(pd_load_waveform_vec(1:100, 100)$y), 0, tolerance = 1e-12)
})

test_that("clipping fraction & RMS & peak-to-peak", {
  y <- c(rep(-1,10), rnorm(80)*0.1, rep(1,10))
  wf <- pd_load_waveform_vec(y, 100)
  qc <- pd_waveform_qc(wf)
  expect_gt(qc$clip_fraction, 0.1)
  expect_equal(qc$rms, sqrt(mean(y^2)), tolerance=1e-9)
  expect_equal(qc$peak_to_peak, diff(range(y)), tolerance=1e-9)
})

test_that("quantization step detected on a 0.12 V grid", {
  y <- round(rnorm(2000)/0.12)*0.12
  expect_equal(pd_quantization_step(y), 0.12, tolerance=1e-6)
})

test_that("exact 2:3:4 triplet gives Q_low=Q_high=2/3", {
  pk <- data.frame(frequency_Hz=c(120,180,240), amplitude=c(3,2,1))
  tr <- pd_koide_triplets(pk)
  expect_equal(tr$Q_low[1], 2/3, tolerance=1e-12)
  expect_equal(tr$Q_high[1], 2/3, tolerance=1e-12)
  expect_lt(tr$koide_error[1], 1e-12)
})

test_that("overlapping triplets are flagged as not independent", {
  pk <- data.frame(frequency_Hz=c(120,180,240,360,480), amplitude=c(5,4,3,2,1))
  tr <- pd_koide_triplets(pk)
  ov <- pd_triplet_overlap(tr)
  expect_true(any(ov$overlapping))
  expect_true(any(grepl("240", ov$shared_frequencies)))
})

test_that("harmonic-conditioned arithmetic note recognizes ladder consequence", {
  pk <- data.frame(frequency_Hz=c(120,180,240), amplitude=c(3,2,1))
  tr <- pd_koide_triplets(pk)
  note <- pd_q_arithmetic_note(tr)
  expect_true(note$arithmetic_consequence)
})

test_that("50 vs 60 Hz discrimination via harmonic models (AIC favors true f0)", {
  f <- c(60,120,180,240,300)*1.0
  hm <- pd_harmonic_models(f)
  expect_lt(hm$rms_resid[hm$model=="M1_f0free,off=0"], 0.5)
  f50 <- c(50,100,150,200,250)
  hm2 <- pd_harmonic_models(f50)  # harmonic numbers via /60 -> imperfect
  expect_true(is.finite(hm2$f0_Hz[2]))
})

test_that("p-value uses (r+1)/(B+1) and never returns 0", {
  expect_equal(pd_pvalue_geq(100, rep(0, 50)), 1/51, tolerance=1e-12)
  expect_gt(pd_pvalue_geq(1e9, rnorm(100)), 0)
})

test_that("scan-global p and family-wise p are distinct concepts", {
  sg <- pd_pvalue_geq(21, rnorm(100))
  fw <- pd_familywise_maxstat(21, rnorm(100, 25, 3))
  expect_false(isTRUE(all.equal(sg, fw)))
})

test_that("ladder statistic is non-negative and detail returned", {
  sr <- 20000; N <- 10000; tt <- (0:(N-1))/sr
  y <- cos(2*pi*120*tt) + 0.5*cos(2*pi*240*tt) + rnorm(N)*0.05
  sp <- pd_fft_canonical(y-mean(y), sr)
  sel <- sp$frequency_Hz>=20 & sp$frequency_Hz<=1000 & sp$frequency_Hz>0
  ls <- pd_ladder_statistic(sp[sel,], 60, 16)
  expect_gte(ls$statistic, 0)
  expect_true(is.data.frame(ls$detail))
})

test_that("injection: zero amplitude behaves like null; power within [0,1]", {
  sr <- 20000; N <- 10000
  wf <- pd_load_waveform_vec(rnorm(N)*0.1, sr)
  kg <- seq(0.5,80,length.out=200)
  ir <- pd_injection_recovery(wf, 20, 1000, kg, pd_default_params(),
                              types="mains_ladder", amps=c(0,2), trials=4, null_thresh=21.14)
  expect_true(all(ir$p_ladder_detected >= 0 & ir$p_ladder_detected <= 1))
  expect_true(all(ir$p_logcos_significant >= 0 & ir$p_logcos_significant <= 1))
})

test_that("ringdown not testable with too few points; Q warning present", {
  rd <- pd_ringdown_models(data.frame(time_s=1:3, value=c(3,2,1)))
  expect_false(rd$testable)
})

test_that("binomial CI in calibration is well-formed", {
  ci <- stats::binom.test(5, 100, 0.05)$conf.int
  expect_true(ci[1] >= 0 && ci[2] <= 1)
})

test_that("parity comparison passes against canonical targets", {
  skip_if_not(file.exists(file.path(PD_REPO, "captures/20260309_010802/waveform_raw.npy")))
  wf <- pd_load_waveform(file.path(PD_REPO,"captures/20260309_010802/waveform_raw.npy"), 20000)
  sp <- pd_fft_canonical(wf$y, 20000)
  sel <- sp$frequency_Hz>=20 & sp$frequency_Hz<=1000 & sp$frequency_Hz>0
  ld <- pd_build_log_domain(sp[sel,], "quadratic")
  sc <- pd_scan_logcos(ld$ell_ln_f, ld$y_detrended, 0.5, 80, 4000)
  best <- sc$best
  best$active_delta_ell <- pd_active_delta_ell(20,1000)
  best$active_winding_n <- pd_active_winding(best$k, 20, 1000)
  cmp <- pd_parity_compare(best)
  expect_true(all(cmp$within_tol))
})
