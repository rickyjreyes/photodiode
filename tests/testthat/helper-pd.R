# Source all photodiode R modules for testing.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
PD_REPO <- Sys.getenv("PD_REPO_ROOT", unset = ".")
if (!dir.exists(file.path(PD_REPO, "R"))) {
  for (cand in c(".", "..", "../..")) if (dir.exists(file.path(cand, "R"))) { PD_REPO <- cand; break }
}
.pd_mods <- c("photodiode_io.R","spectral_estimators.R","logcos_scan.R",
              "peak_estimation.R","harmonic_ladder.R","ratio_geometry.R",
              "waveform_qc.R","logcos_nulls.R","time_frequency_analysis.R",
              "persistence_analysis.R","ringdown_models.R","control_comparison.R",
              "window_sensitivity.R","multiple_testing.R","null_calibration.R",
              "injection_recovery.R","python_r_parity.R","photodiode_registry.R",
              "build_claim_matrix.R","check_photodiode_dependencies.R")
for (m in .pd_mods) sys.source(file.path(PD_REPO, "R", m), envir = globalenv())

# small helper: build a waveform list from a numeric vector
pd_load_waveform_vec <- function(y, sr) {
  yd <- y - mean(y)
  list(t = seq(0, by = 1/sr, length.out = length(y)), y_raw = y, y = yd,
       n = length(y), sample_rate = sr, source = "vec")
}

# write a tiny float64 1-D .npy for IO tests
pd_test_write_npy <- function(path, x) {
  con <- file(path, "wb"); on.exit(close(con))
  writeBin(as.raw(c(0x93,0x4e,0x55,0x4d,0x50,0x59)), con)
  writeBin(as.raw(c(1,0)), con)
  hdr <- sprintf("{'descr': '<f8', 'fortran_order': False, 'shape': (%d,), }", length(x))
  total <- 10 + nchar(hdr) + 1
  pad <- (64 - (total %% 64)) %% 64
  hdr <- paste0(hdr, strrep(" ", pad), "\n")
  writeBin(as.integer(nchar(hdr)), con, size = 2L, endian = "little")
  writeChar(hdr, con, eos = NULL)
  writeBin(as.numeric(x), con, size = 8L, endian = "little")
}
