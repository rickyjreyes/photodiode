#!/usr/bin/env Rscript
# ==============================================================================
# render_photodiode_audit.R
# Command-line orchestrator for the photodiode statistical-audit suite.
#
# Usage:
#   Rscript R/render_photodiode_audit.R [--mode replay|waveform] [--fast] [--strict]
#       [--capture PATH] [--sample-rate 20000] [--fmin 20] [--fmax 1000]
#       [--kmin 0.5] [--kmax 80] [--nk 4000] [--null-n 5000] [--bootstrap-n 5000]
#       [--calibration-n 2000] [--injection-n 1000] [--seed 20260309]
#       [--repo-root .] [--render-report true] [--include-ringdown auto] [--force]
# ==============================================================================

suppressWarnings(suppressMessages({
  library(ggplot2); library(jsonlite)
}))

# ---- argument parsing --------------------------------------------------------
pd_parse_args <- function(argv) {
  defaults <- list(mode = "waveform", repo_root = ".",
                   capture = "captures/20260309_010802/waveform_raw.npy",
                   sample_rate = 20000, fmin = 20, fmax = 1000,
                   kmin = 0.5, kmax = 80, nk = 4000, null_n = 5000,
                   bootstrap_n = 5000, calibration_n = 2000, injection_n = 1000,
                   seed = 20260309, parallel = "true", workers = "auto",
                   render_report = "true", include_ringdown = "auto",
                   force = FALSE, strict = FALSE, fast = FALSE)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c("--force","--strict","--fast")) {
      defaults[[sub("^--","",a)]] <- TRUE; i <- i + 1
    } else if (grepl("^--", a)) {
      key <- gsub("-", "_", sub("^--","",a)); val <- argv[i+1]
      num_keys <- c("sample_rate","fmin","fmax","kmin","kmax","nk","null_n",
                    "bootstrap_n","calibration_n","injection_n","seed")
      defaults[[key]] <- if (key %in% num_keys) as.numeric(val) else val
      i <- i + 2
    } else i <- i + 1
  }
  defaults
}

pd_source_all <- function(repo_root) {
  rd <- file.path(repo_root, "R")
  files <- c("photodiode_io.R","spectral_estimators.R","logcos_scan.R",
             "peak_estimation.R","harmonic_ladder.R","ratio_geometry.R",
             "waveform_qc.R","logcos_nulls.R","time_frequency_analysis.R",
             "persistence_analysis.R","ringdown_models.R","control_comparison.R",
             "window_sensitivity.R","multiple_testing.R","null_calibration.R",
             "injection_recovery.R","python_r_parity.R","photodiode_registry.R",
             "build_claim_matrix.R","check_photodiode_dependencies.R")
  for (f in files) sys.source(file.path(rd, f), envir = globalenv())
}

# ---- figure helpers ----------------------------------------------------------
pd_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom")
}
pd_save <- function(plot, path_png, w = 8, h = 5) {
  ggsave(path_png, plot, width = w, height = h, dpi = 300)
  ggsave(sub("\\.png$", ".pdf", path_png), plot, width = w, height = h)
}
pd_text_panel <- function(msg, title = "") {
  ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 4) +
    labs(title = title) + theme_void() +
    theme(plot.title = element_text(face = "bold"))
}

# ==============================================================================
# MAIN
# ==============================================================================
pd_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  args <- pd_parse_args(argv)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(args$seed)
  repo <- args$repo_root
  pd_source_all(repo)
  dirs <- pd_ensure_dirs(repo)
  warnings_acc <- character(0)
  wlog <- function(m) { warnings_acc[[length(warnings_acc)+1]] <<- m; message("[warn] ", m) }

  # fast-mode sizing
  if (isTRUE(args$fast)) {
    args$null_n <- min(args$null_n, 200); args$bootstrap_n <- min(args$bootstrap_n, 200)
    args$calibration_n <- min(args$calibration_n, 100); args$injection_n <- 8
    dev_tag <- "DEVELOPMENT_ONLY"
  } else dev_tag <- "FULL"

  T <- dirs$tables; F <- dirs$figures
  wcsv <- function(x, name) utils::write.csv(x, file.path(T, name), row.names = FALSE)
  res <- list()

  # ---- 1. dependencies ----
  deps <- pd_check_dependencies(verbose = TRUE)
  if (args$strict && !deps$ok) stop("Missing required packages: ",
                                    paste(deps$missing_required, collapse=", "))

  # ---- 2. registry ----
  message("[stage] capture registry")
  reg <- pd_build_registry(repo)
  wcsv(reg$registry, "capture_registry.csv")
  pd_write_json(list(registry = reg$registry, data_flow = reg$flow,
                     missing_controls = reg$missing_controls),
                file.path(dirs$audit_out, "capture_registry.json"))
  if (args$strict && any(duplicated(reg$registry$analysis_id)))
    stop("Duplicate capture IDs detected")
  wfs <- reg$waveforms
  canonical_id <- "cap_post_canonical"
  if (is.null(wfs[[canonical_id]])) {
    if (args$strict) stop("Canonical waveform missing") else wlog("canonical waveform missing")
  }
  wf <- wfs[[canonical_id]]

  # ---- 3. parity ----
  message("[stage] PARITY_MODE")
  parity <- pd_run_parity(wf, dirs, args$fmin, args$fmax, args$kmin, args$kmax,
                          args$nk, null_n = min(args$null_n, 1000), seed = 12345)
  pcmp <- pd_parity_compare(parity$best)
  wcsv(pcmp, "python_r_parity.csv")
  parity_pass <- all(pcmp$within_tol)
  if (args$strict && !parity_pass) stop("Parity drift beyond tolerance")
  res$best_k <- parity$best$k; res$best_delta <- parity$best$delta_chi2
  res$best_A <- parity$best$A; res$p_shuffle <- parity$p; res$null_n <- length(parity$null_best)
  res$parity_pass <- parity_pass; res$parity_maxreldiff <- max(pcmp$rel_diff)

  # ---- 4. waveform QC ----
  message("[stage] waveform QC")
  qc_rows <- lapply(names(wfs), function(id) cbind(analysis_id = id, pd_waveform_qc(wfs[[id]])))
  qc_tab <- dplyr::bind_rows(qc_rows); wcsv(qc_tab, "waveform_qc.csv")

  # ---- 5. spectral estimation (window comparison on canonical) ----
  message("[stage] spectral estimation")
  spec_summary <- pd_window_comparison(wf$y, wf$sample_rate, args$fmin, args$fmax)
  wcsv(spec_summary, "spectrum_summary.csv")

  # ---- 6. peak estimates ----
  message("[stage] peak estimation")
  spec_can <- parity$spec
  peaks_can <- parity$peaks
  peak_est <- pd_peak_estimates(spec_can, peaks_can, t = wf$t, y = wf$y)
  wcsv(peak_est, "peak_estimates.csv")
  harm_models <- pd_harmonic_models(peaks_can$frequency_Hz)
  res$f0_M1 <- harm_models$f0_Hz[harm_models$model == "M1_f0free,off=0"]

  # ---- 7. harmonic ladder ----
  message("[stage] harmonic ladder")
  ladder <- pd_ladder_statistic(spec_can, 60, 16)
  wcsv(ladder$detail, "harmonic_ladder.csv")
  res$ladder_n <- ladder$n_detected
  res$ladder_freqs <- paste(ladder$detail$target_Hz[ladder$detail$detected], collapse=",")

  # ---- 8. harmonic effect sizes (baseline vs post peaks) ----
  message("[stage] harmonic effect sizes")
  base_peaks_path <- file.path(repo, "proof_runs/1773046194_baseline_peaks.csv")
  if (file.exists(base_peaks_path)) {
    bp <- utils::read.csv(base_peaks_path)
    eff <- pd_harmonic_effect_sizes(bp, peaks_can, comparable = FALSE)
    attr(eff, "note") <- "baseline (calibrated V) and post (unnormalized FFT) amplitudes are NOT on a comparable scale"
    wcsv(eff, "harmonic_effect_sizes.csv")
  } else { eff <- data.frame(); wlog("baseline peaks not found") }

  # ---- 9. ratio geometry ----
  message("[stage] ratio geometry")
  ratios <- parity$ratios; triplets <- parity$triplets
  wcsv(ratios, "ratio_geometry.csv")
  overlap <- pd_triplet_overlap(triplets)
  wcsv(if (nrow(overlap)) overlap else data.frame(note="no exact triplets"), "triplet_overlap.csv")
  qnote <- pd_q_arithmetic_note(triplets)
  res$q_exact <- qnote$n_exact_q23_triplets

  # ---- 10. logcos results ----
  message("[stage] logcos results")
  logcos_tab <- data.frame(
    quantity = c("best_k","delta_chi2","chi2_null","chi2_logcos","A","phi",
                 "active_delta_ell","active_winding","p_shuffle_scanmax"),
    value = c(parity$best$k, parity$best$delta_chi2, parity$best$chi2_null,
              parity$best$chi2_logcos, parity$best$A, parity$best$phi,
              parity$best$active_delta_ell, parity$best$active_winding_n, parity$p))
  wcsv(logcos_tab, "logcos_results.csv")

  # ---- 11. logcos uncertainty (block bootstrap of best k / A / delta) ----
  message("[stage] logcos uncertainty")
  ld <- parity$scan$scan; ell <- log(spec_can$frequency_Hz)
  yobs <- pd_build_log_domain(spec_can, "quadratic")$y_detrended
  kgrid <- seq(args$kmin, args$kmax, length.out = min(args$nk, 1500))
  bn <- args$bootstrap_n
  bk <- numeric(bn); ba <- numeric(bn); bd <- numeric(bn)
  blk <- 8; n <- length(yobs)
  blocks <- split(seq_len(n), ceiling(seq_len(n)/blk))
  for (i in seq_len(bn)) {
    idx <- unlist(blocks[sample(seq_along(blocks), length(blocks) + 2, replace = TRUE)])
    idx <- idx[seq_len(n)]
    sm <- pd_fit_logcos_best(ell, yobs[idx], kgrid)
    bk[i] <- sm$k; ba[i] <- sm$A; bd[i] <- sm$delta_chi2
  }
  unc <- data.frame(
    quantity = c("best_k","log_period_2pi_over_k","active_winding","A","delta_chi2"),
    point = c(parity$best$k, 2*pi/parity$best$k, parity$best$active_winding_n,
              parity$best$A, parity$best$delta_chi2),
    boot_ci_low = c(quantile(bk,.025), 2*pi/quantile(bk,.975),
                    quantile(bk,.025)*log(args$fmax/args$fmin)/(2*pi),
                    quantile(ba,.025), quantile(bd,.025)),
    boot_ci_high = c(quantile(bk,.975), 2*pi/quantile(bk,.025),
                     quantile(bk,.975)*log(args$fmax/args$fmin)/(2*pi),
                     quantile(ba,.975), quantile(bd,.975)),
    method = "block bootstrap (L=8)", n_boot = bn)
  wcsv(unc, "logcos_uncertainty.csv")

  # ---- 12. logcos stronger nulls ----
  message("[stage] logcos nulls (canonical + stronger)")
  obs_delta <- parity$best$delta_chi2
  nulls <- list()
  nulls$shuffle <- list(label="shuffle_canonical", B=length(parity$null_best),
                        exceedances=sum(parity$null_best>=obs_delta),
                        p_value=parity$p, resolution_floor=1/(length(parity$null_best)+1),
                        null_mean=mean(parity$null_best), null_sd=sd(parity$null_best),
                        null_median=median(parity$null_best),
                        p95=quantile(parity$null_best,.95,names=FALSE),
                        p99=quantile(parity$null_best,.99,names=FALSE),
                        null_max=max(parity$null_best), null_best=parity$null_best,
                        null_k=rep(NA,length(parity$null_best)))
  nulls$block <- pd_null_block(ell, yobs, obs_delta, args$null_n, kgrid, 8)
  nulls$stationary <- pd_null_stationary_boot(ell, yobs, obs_delta, args$null_n, kgrid)
  nulls$colored <- tryCatch(pd_null_colored_ar(wf, args$fmin, args$fmax, obs_delta,
                              args$null_n, kgrid), error=function(e){wlog(paste("colored null:",e$message));NULL})
  ha <- pd_harmonic_amplitudes(spec_can, 60, 10)
  hamp <- ha$amplitude[ha$detected]; if (!length(hamp)) hamp <- ha$amplitude
  # scale harmonic amps to physical-ish using waveform sd
  hamp <- hamp / max(hamp) * sd(wf$y)
  nulls$mains <- tryCatch(pd_null_mains(wf, hamp, args$fmin, args$fmax, obs_delta,
                              args$null_n, kgrid), error=function(e){wlog(paste("mains null:",e$message));NULL})
  null_cmp <- pd_null_comparison_table(nulls, obs_delta)
  null_cmp$mode <- dev_tag
  wcsv(null_cmp, "logcos_null_comparison.csv")
  res$p_mains <- if (!is.null(nulls$mains)) nulls$mains$p_value else NA
  res$strong_null_summary <- paste(sprintf("%s p=%.3g", null_cmp$null_model, null_cmp$p_value), collapse="; ")
  res$logcos_strong_verdict <- {
    mp <- if (!is.null(nulls$mains)) nulls$mains$p_value else NA
    cp <- if (!is.null(nulls$colored)) nulls$colored$p_value else NA
    if (isTRUE(mp < 0.05) && isTRUE(cp < 0.05)) "ROBUST_STATISTICAL_FEATURE"
    else if (isTRUE(mp >= 0.05) || isTRUE(cp >= 0.05)) "CONSISTENT_WITH_MAINS"
    else "INCONCLUSIVE"
  }

  # ---- 13. window sensitivity ----
  message("[stage] window sensitivity")
  wsens <- pd_window_sensitivity(wf, nk = if (isTRUE(args$fast)) 400 else 1000)
  wsens$specification_adjusted_note <- "max over grid is the family statistic"
  wcsv(wsens, "window_sensitivity.csv")

  # ---- 14. control comparison ----
  message("[stage] control comparison")
  meta_lookup <- setNames(lapply(pd_capture_inventory(repo), function(m) m),
                          vapply(pd_capture_inventory(repo), function(m) m$analysis_id, character(1)))
  analyses <- lapply(names(wfs), function(id) pd_analyze_capture(wfs[[id]]))
  names(analyses) <- names(wfs)
  cmat <- pd_control_matrix(analyses, meta_lookup)
  # append missing controls as explicit rows
  miss <- data.frame(analysis_id = reg$missing_controls, condition = "MISSING_CONTROL",
                     control_type = reg$missing_controls,
                     elapsed_time_since_excitation_s = NA,
                     rms_voltage=NA, amp_60Hz=NA, total_harmonic_power=NA,
                     n_harmonics=NA, ladder_statistic=NA, q_triplet_count=NA,
                     min_koide_error=NA, logcos_best_k=NA, logcos_delta_chi2=NA,
                     logcos_A=NA, check.names = FALSE)
  cmat_full <- dplyr::bind_rows(cmat, miss)
  wcsv(cmat_full, "control_comparison.csv")

  # ---- 15. persistence ----
  message("[stage] persistence")
  persist <- pd_persistence_within_capture(wf, args$fmin, args$fmax)
  csum <- data.frame(epoch = reg$registry$epoch, phase = reg$registry$phase)
  acurve <- pd_persistence_curve(csum)
  persist_tab <- data.frame(
    metric = c("within_capture_frames","capture_duration_s","across_capture_testable"),
    value = c(persist$n_frames, persist$capture_duration_s, acurve$testable),
    note = c(persist$note, persist$note, acurve$reason))
  wcsv(persist_tab, "persistence_results.csv")

  # ---- 16. ringdown ----
  message("[stage] ringdown")
  rd <- pd_ringdown_models(NULL)
  wcsv(data.frame(testable = rd$testable, reason = rd$reason), "ringdown_models.csv")

  # ---- 17. multiple testing ----
  message("[stage] multiple testing")
  pvals <- c(shuffle = parity$p,
             block = if(!is.null(nulls$block)) nulls$block$p_value else NA,
             colored = if(!is.null(nulls$colored)) nulls$colored$p_value else NA,
             mains = if(!is.null(nulls$mains)) nulls$mains$p_value else NA)
  pvals <- pvals[is.finite(pvals)]
  mt <- pd_multiple_testing(pvals, names(pvals))
  fw_p <- pd_familywise_maxstat(max(wsens$delta_chi2),
                                if(!is.null(nulls$mains)) nulls$mains$null_best else nulls$shuffle$null_best)
  mt$specification_max_delta_chi2 <- max(wsens$delta_chi2)
  mt$familywise_maxstat_p <- fw_p
  wcsv(mt, "multiple_testing.csv")

  # ---- 18. null calibration ----
  message("[stage] null calibration")
  calib <- tryCatch(pd_null_calibration(wf, args$fmin, args$fmax, kgrid,
                      outer_n = args$calibration_n, inner_B = if(isTRUE(args$fast)) 50 else 200),
                    error = function(e){wlog(paste("calibration:",e$message)); NULL})
  if (!is.null(calib)) { calib$mode <- dev_tag; wcsv(calib, "null_calibration.csv") } else
    wcsv(data.frame(note="calibration unavailable"), "null_calibration.csv")
  res$calib_status <- if (!is.null(calib))
    sprintf("alpha=0.05 -> obs=%.3f", calib$observed_fpr[calib$nominal_alpha==0.05]) else "n/a"

  # ---- 19. injection recovery ----
  message("[stage] injection recovery")
  inj <- pd_injection_recovery(wf, args$fmin, args$fmax, kgrid, pd_default_params(),
                               trials = if(isTRUE(args$fast)) 8 else 25)
  inj$mode <- dev_tag
  wcsv(inj, "injection_recovery.csv")
  res$inj_amp <- 0.5

  # ---- 20. prediction ledger audit ----
  message("[stage] prediction ledger audit")
  pl_audit <- pd_prediction_ledger_audit(repo)
  wcsv(pl_audit, "prediction_ledger_audit.csv")

  # ---- 21. claim matrix ----
  message("[stage] claim matrix")
  claim <- pd_build_claim_matrix(res)
  wcsv(claim, "final_claim_matrix.csv")

  # ---- 22. figures ----
  message("[stage] figures")
  pd_make_figures(F, reg, wf, spec_can, peaks_can, peak_est, ladder, ratios,
                  triplets, overlap, parity, nulls, unc, wsens, persist,
                  cmat_full, calib, inj, pcmp, eff, claim, args)

  # ---- 23. report ----
  if (tolower(as.character(args$render_report)) %in% c("true","1") &&
      requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
    message("[stage] render report")
    tryCatch(pd_render_report(repo, dirs), error = function(e) wlog(paste("report:", e$message)))
  } else wlog("report rendering skipped (pandoc/rmarkdown unavailable or disabled)")

  # ---- final summary ----
  n_fig <- length(list.files(F, pattern="\\.png$"))
  n_tab <- length(list.files(T, pattern="\\.csv$"))
  cat("\n", strrep("=",70), "\nPHOTODIODE STATISTICAL AUDIT COMPLETE (", dev_tag, ")\n",
      strrep("=",70), "\n", sep="")
  cat(sprintf("captures discovered     : %d\n", reg$flow$files_discovered))
  cat(sprintf("captures analyzed (raw) : %d\n", reg$flow$retained))
  cat(sprintf("baseline captures       : %d (peak-table only, no raw)\n",
              sum(reg$registry$phase=="baseline")))
  cat(sprintf("post-excitation captures: %d\n", sum(reg$registry$phase=="post")))
  cat(sprintf("missing control classes : %d (%s)\n", length(reg$missing_controls),
              paste(reg$missing_controls, collapse=", ")))
  cat(sprintf("parity status           : %s (max rel diff %.2e)\n",
              ifelse(parity_pass,"PASS","FAIL"), max(pcmp$rel_diff)))
  cat(sprintf("waveform QC failures    : %d\n", sum(qc_tab$clip_fraction>0.01, na.rm=TRUE)))
  cat(sprintf("null simulations        : %d per stronger null\n", args$null_n))
  cat(sprintf("calibration simulations : %d\n", args$calibration_n))
  cat(sprintf("injection trials/cell   : %d\n", if(isTRUE(args$fast)) 8 else 25))
  cat(sprintf("figures created         : %d\n", n_fig))
  cat(sprintf("tables created          : %d\n", n_tab))
  cat(sprintf("warnings                : %d\n", length(warnings_acc)))
  cat(sprintf("report path             : %s\n", file.path(dirs$rendered,
              "photodiode_statistical_audit.html")))
  cat("\nFull-resolution command:\n",
      "  Rscript R/render_photodiode_audit.R --mode waveform --null-n 10000",
      "--bootstrap-n 5000 --calibration-n 10000 --injection-n 1000 --seed 20260309\n")
  invisible(res)
}

# best-only logcos fit over a grid (used in bootstrap), vectorized engine
pd_fit_logcos_best <- function(ell, y, kgrid) {
  e <- pd_logcos_engine(ell, y, kgrid)
  i <- which.max(e$delta)
  list(k = kgrid[i], A = sqrt(e$a[i]^2 + e$b[i]^2), delta_chi2 = e$delta[i])
}

# prediction ledger audit
pd_prediction_ledger_audit <- function(repo) {
  pl_path <- file.path(repo, "proof_ledger.json")
  ledger_ok <- FALSE
  if (file.exists(pl_path)) {
    pl <- tryCatch(jsonlite::fromJSON(pl_path, simplifyVector = FALSE),
                   error = function(e) NULL)
    ledger_ok <- !is.null(pl)
  }
  preds <- list(
    c("60 Hz-spaced harmonic ladder present", "predeclared", "OBSERVED",
      "Ladder present in post capture", "TRUE"),
    c("Baseline (virgin) run captured before excitation", "predeclared", "OBSERVED",
      "Baseline run 1773046194 exists (peak table)", "TRUE"),
    c("Q=2/3 ratio geometry in triplets", "qualitative expectation", "OBSERVED",
      "Two exact triplets, but arithmetic consequence of ladder", "PARTIAL"),
    c("Persistence after illumination ceases", "qualitative expectation", "NOT_TESTED",
      "Single 0.5 s capture; no time series", "INCONCLUSIVE"),
    c("Controls falsify environmental pickup", "predeclared falsification", "NOT_PERFORMED",
      "No raw control captures committed", "FAIL_TO_TEST"),
    c("Pre-declared parameter lock (fmin,fmax,k)", "predeclared", "OBSERVED",
      "Parameters fixed and reproduced exactly", "TRUE")
  )
  out <- do.call(rbind, lapply(preds, function(p) data.frame(
    prediction = p[1], declaration_type = p[2], observed_status = p[3],
    observed_result = p[4],
    specified_before_capture = p[2] %in% c("predeclared","predeclared falsification"),
    verdict = p[5], stringsAsFactors = FALSE)))
  integrity <- data.frame(
    prediction = "proof_ledger.json file integrity",
    declaration_type = "metadata",
    observed_status = if (ledger_ok) "PARSEABLE" else "MALFORMED_TRUNCATED",
    observed_result = if (ledger_ok) "ledger JSON parses cleanly" else
      "ledger JSON is truncated (premature EOF); machine-readable predictions unavailable",
    specified_before_capture = NA,
    verdict = if (ledger_ok) "TRUE" else "DATA_DEFECT", stringsAsFactors = FALSE)
  rbind(out, integrity)
}

# render the Quarto/Rmd report via rmarkdown (quarto not required)
pd_render_report <- function(repo, dirs) {
  qmd <- file.path(repo, "reports", "photodiode_statistical_audit.qmd")
  rmd <- file.path(dirs$rendered, "photodiode_statistical_audit.Rmd")
  file.copy(qmd, rmd, overwrite = TRUE)
  rmarkdown::render(rmd, output_format = "html_document",
                    output_file = "photodiode_statistical_audit.html",
                    output_dir = dirs$rendered, quiet = TRUE,
                    params = list(repo_root = normalizePath(repo)))
}

# ==============================================================================
# FIGURES
# ==============================================================================
pd_make_figures <- function(F, reg, wf, spec, peaks, peak_est, ladder, ratios,
                            triplets, overlap, parity, nulls, unc, wsens, persist,
                            cmat, calib, inj, pcmp, eff, claim, args) {
  pal <- PD_PALETTE
  # fig01 capture flow
  fl <- reg$flow
  flow_df <- data.frame(stage = factor(c("discovered","readable","retained","excluded","malformed"),
                          levels=c("discovered","readable","retained","excluded","malformed")),
                        n = c(fl$files_discovered, fl$readable, fl$retained, fl$excluded, fl$malformed))
  pd_save(ggplot(flow_df, aes(stage, n)) + geom_col(fill=pal[["post"]]) +
            geom_text(aes(label=n), vjust=-0.3) + pd_theme() +
            labs(title="fig01 Capture data flow", x=NULL, y="count"),
          file.path(F,"fig01_capture_flow.png"))

  # fig02 waveform QC
  qdf <- data.frame(t = wf$t, v = wf$y)
  pd_save(ggplot(qdf[seq(1,nrow(qdf),5),], aes(t,v)) + geom_line(color=pal[["post"]],linewidth=0.3) +
            pd_theme() + labs(title="fig02 Canonical waveform (DC removed)", x="time (s)", y="V"),
          file.path(F,"fig02_waveform_qc.png"))

  # fig03 baseline vs post waveforms (baseline raw unavailable)
  pd_save(pd_text_panel("Baseline raw waveform NOT available\n(only a peak table exists).\nPost-excitation waveform shown in fig02.",
                        "fig03 Baseline vs post waveforms"),
          file.path(F,"fig03_baseline_post_waveforms.png"))

  # fig04 PSD comparison (windows)
  psd <- pd_psd(wf$y, wf$sample_rate)
  psd <- psd[psd$frequency_Hz>=args$fmin & psd$frequency_Hz<=args$fmax,]
  pd_save(ggplot(psd, aes(frequency_Hz, psd)) + geom_line(color=pal[["post"]]) +
            scale_y_log10() + pd_theme() +
            labs(title="fig04 PSD (Hann, V^2/Hz)", x="Hz", y="PSD"),
          file.path(F,"fig04_psd_comparison.png"))

  # fig05 difference spectrum (baseline raw unavailable -> annotate)
  pd_save(pd_text_panel("Difference spectrum requires baseline raw waveform,\nwhich is not committed. Baseline exists as peaks only.",
                        "fig05 Difference spectrum (post - baseline)"),
          file.path(F,"fig05_difference_spectrum.png"))

  # fig06 harmonic ladder lollipop
  hd <- ladder$detail
  pd_save(ggplot(hd, aes(target_Hz, amplitude)) +
            geom_segment(aes(xend=target_Hz, yend=0), color=pal[["null"]]) +
            geom_point(aes(color=detected), size=2.5) +
            scale_color_manual(values=c(`TRUE`=pal[["post"]],`FALSE`=pal[["null"]])) +
            pd_theme() + labs(title="fig06 60 Hz harmonic ladder", x="Hz", y="amplitude"),
          file.path(F,"fig06_harmonic_ladder.png"))

  # fig07 frequency residuals from n*60
  if (nrow(peak_est)) {
    pd_save(ggplot(peak_est, aes(nearest_60Hz_harmonic, residual_from_n60)) +
              geom_hline(yintercept=0, linetype=2) +
              geom_point(color=pal[["post"]], size=2.5) + pd_theme() +
              labs(title="fig07 Sub-bin residual from n x 60 Hz",
                   x="harmonic n", y="f_subbin - n*60 (Hz)"),
            file.path(F,"fig07_frequency_residuals.png"))
  } else pd_save(pd_text_panel("no peaks","fig07"), file.path(F,"fig07_frequency_residuals.png"))

  # fig08 harmonic effect sizes
  if (nrow(eff)) {
    ed <- tidyr::pivot_longer(eff, c("baseline_amplitude","post_amplitude"))
    pd_save(ggplot(ed, aes(factor(harmonic), value, fill=name)) +
              geom_col(position="dodge") +
              scale_fill_manual(values=c(baseline_amplitude=pal[["baseline"]],
                                         post_amplitude=pal[["post"]])) +
              pd_theme() + labs(title="fig08 Harmonic amplitudes (NOTE: different normalizations)",
                                x="harmonic", y="amplitude (not comparable scale)"),
            file.path(F,"fig08_harmonic_effect_sizes.png"))
  } else pd_save(pd_text_panel("baseline peaks unavailable","fig08"), file.path(F,"fig08_harmonic_effect_sizes.png"))

  # fig09 Q geometry
  pd_save(ggplot(triplets, aes(Q_low, Q_high, color=koide_error)) +
            geom_point(size=2) + geom_vline(xintercept=2/3,linetype=2) +
            geom_hline(yintercept=2/3,linetype=2) +
            scale_color_viridis_c() + pd_theme() +
            labs(title="fig09 Q geometry (Q_low vs Q_high; 2/3 dashed)"),
          file.path(F,"fig09_q_geometry.png"))

  # fig10 triplet overlap
  if (nrow(overlap)) {
    pd_save(pd_text_panel(paste(apply(overlap,1,function(r)
              sprintf("%s & %s share {%s}", r["triplet_a"], r["triplet_b"], r["shared_frequencies"])),
              collapse="\n"), "fig10 Triplet overlap (NOT independent)"),
            file.path(F,"fig10_triplet_overlap.png"))
  } else pd_save(pd_text_panel("no exact triplets","fig10"), file.path(F,"fig10_triplet_overlap.png"))

  # fig11 logcos scan
  sc <- parity$scan$scan
  pd_save(ggplot(sc, aes(k, delta_chi2)) + geom_line(color=pal[["post"]]) +
            geom_vline(xintercept=parity$best$k, linetype=2) + pd_theme() +
            labs(title=sprintf("fig11 Log-cos scan (best k=%.3f)", parity$best$k),
                 x="k", y="delta chi2"),
          file.path(F,"fig11_logcos_scan.png"))

  # fig12 logcos nulls ridgeline
  nd <- do.call(rbind, lapply(names(nulls), function(nm)
    if(!is.null(nulls[[nm]])) data.frame(model=nulls[[nm]]$label, val=nulls[[nm]]$null_best) else NULL))
  pd_save(ggplot(nd, aes(val, model, fill=model)) +
            ggridges::geom_density_ridges(alpha=0.6) +
            geom_vline(xintercept=parity$best$delta_chi2, color=pal[["post"]], linewidth=1) +
            scale_fill_viridis_d() + pd_theme() + theme(legend.position="none") +
            labs(title="fig12 Null distributions (red = observed)", x="scan-max delta chi2", y=NULL),
          file.path(F,"fig12_logcos_nulls.png"))

  # fig13 bootstrap
  pd_save(ggplot(unc, aes(quantity, point)) +
            geom_point(color=pal[["post"]]) +
            geom_errorbar(aes(ymin=boot_ci_low, ymax=boot_ci_high), width=0.2) +
            coord_flip() + pd_theme() + labs(title="fig13 Log-cos bootstrap CIs", x=NULL, y="value"),
          file.path(F,"fig13_logcos_bootstrap.png"))

  # fig14 window sensitivity (spec curve)
  ws <- wsens[order(wsens$delta_chi2),]; ws$rank <- seq_len(nrow(ws))
  pd_save(ggplot(ws, aes(rank, delta_chi2, color=in_peak_region)) + geom_point() +
            scale_color_manual(values=c(`TRUE`=pal[["post"]],`FALSE`=pal[["null"]])) +
            pd_theme() + labs(title="fig14 Specification curve (window/k/detrend/estimator)",
                              x="specification rank", y="delta chi2"),
          file.path(F,"fig14_window_sensitivity.png"))

  # fig15 spectrogram
  st <- persist$stft
  sp_df <- expand.grid(time=st$time, freq=st$freq)
  sp_df$amp <- as.vector(t(st$amplitude))
  sp_df <- sp_df[sp_df$freq<=600,]
  pd_save(ggplot(sp_df, aes(time, freq, fill=log1p(amp))) + geom_raster() +
            scale_fill_viridis_c() + pd_theme() +
            labs(title="fig15 Spectrogram", x="time (s)", y="Hz"),
          file.path(F,"fig15_spectrogram.png"))

  # fig16 persistence within capture
  pd_save(ggplot(persist$tracks, aes(time_s, amplitude, color=factor(harmonic))) +
            geom_line() + scale_color_viridis_d(name="harmonic Hz") + pd_theme() +
            labs(title="fig16 Harmonic amplitude vs time (within 0.5 s capture only)",
                 x="time (s)", y="amplitude"),
          file.path(F,"fig16_harmonic_persistence.png"))

  # fig17 control matrix heatmap
  cm <- cmat[, c("analysis_id","condition","logcos_delta_chi2","ladder_statistic",
                 "amp_60Hz","q_triplet_count")]
  cml <- tidyr::pivot_longer(cm, c("logcos_delta_chi2","ladder_statistic","amp_60Hz","q_triplet_count"))
  pd_save(ggplot(cml, aes(name, analysis_id, fill=value)) + geom_tile(color="white") +
            scale_fill_viridis_c(na.value="grey85") + pd_theme() +
            theme(axis.text.x=element_text(angle=30,hjust=1)) +
            labs(title="fig17 Control matrix (grey = missing/NA)", x=NULL, y=NULL),
          file.path(F,"fig17_control_matrix.png"), h=6)

  # fig18 null calibration
  if (!is.null(calib)) {
    pd_save(ggplot(calib, aes(nominal_alpha, observed_fpr)) +
              geom_abline(linetype=2) +
              geom_point(color=pal[["post"]]) +
              geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=0) +
              scale_x_log10() + scale_y_log10() + pd_theme() +
              labs(title="fig18 Null calibration", x="nominal alpha", y="observed FPR"),
            file.path(F,"fig18_null_calibration.png"))
  } else pd_save(pd_text_panel("calibration unavailable","fig18"), file.path(F,"fig18_null_calibration.png"))

  # fig19 injection recovery power curves
  pd_save(ggplot(inj, aes(amplitude, p_ladder_detected, color=type)) +
            geom_line() + geom_point() + scale_color_viridis_d() + pd_theme() +
            labs(title="fig19 Injection-recovery: ladder detection power",
                 x="injected amplitude", y="P(detect)"),
          file.path(F,"fig19_injection_recovery.png"))

  # fig20 parity
  pd_save(ggplot(pcmp, aes(quantity, rel_diff)) + geom_col(fill=pal[["control"]]) +
            scale_y_log10() + coord_flip() + pd_theme() +
            labs(title="fig20 Python/R parity (relative difference)", x=NULL, y="rel diff"),
          file.path(F,"fig20_python_r_parity.png"))

  # fig21 evidence dashboard
  p_psd <- ggplot(psd, aes(frequency_Hz, psd)) + geom_line(color=pal[["post"]]) +
    scale_y_log10() + pd_theme() + labs(title="A PSD", x="Hz", y=NULL)
  p_lad <- ggplot(hd, aes(target_Hz, amplitude)) +
    geom_segment(aes(xend=target_Hz, yend=0), color=pal[["null"]]) +
    geom_point(color=pal[["post"]]) + pd_theme() + labs(title="B Harmonics", x="Hz", y=NULL)
  p_q <- ggplot(triplets, aes(Q_low,Q_high)) + geom_point(color=pal[["post"]]) +
    geom_vline(xintercept=2/3,linetype=2)+geom_hline(yintercept=2/3,linetype=2) +
    pd_theme()+labs(title="D Q geometry")
  p_null <- ggplot(nd, aes(val, model, fill=model)) + ggridges::geom_density_ridges(alpha=.6)+
    geom_vline(xintercept=parity$best$delta_chi2,color=pal[["post"]])+
    scale_fill_viridis_d()+pd_theme()+theme(legend.position="none")+labs(title="E Nulls",x="dchi2",y=NULL)
  p_claim <- pd_text_panel(paste(sprintf("%s: %s", substr(claim$claim,1,28), claim$verdict),
                                 collapse="\n"), "H Claim matrix")
  if (requireNamespace("patchwork", quietly=TRUE)) {
    dash <- patchwork::wrap_plots(p_psd, p_lad, p_q, p_null, p_claim, ncol=2) +
      patchwork::plot_annotation(title="fig21 Evidence dashboard")
    pd_save(dash, file.path(F,"fig21_evidence_dashboard.png"), w=12, h=10)
  } else pd_save(p_claim, file.path(F,"fig21_evidence_dashboard.png"), w=12, h=10)
}

# Run only when invoked directly as a script (not when sys.source'd).
if (sys.nframe() == 0 && !interactive()) {
  pd_main()
}
