# ==============================================================================
# photodiode_registry.R
# Capture discovery, metadata accounting, and the data-flow ledger.
# Missing metadata is recorded explicitly and never inferred.
# ==============================================================================

# Hand-curated inventory of what the repository actually contains. Phase labels
# come from filenames / proof_ledger.json / README and are marked for whether
# they were declared before analysis. Raw availability is explicit.
pd_capture_inventory <- function(repo_root = ".") {
  rp <- function(...) file.path(repo_root, ...)
  inv <- list(
    list(analysis_id = "cap_post_canonical",
         path = rp("captures/20260309_010802/waveform_raw.npy"),
         run_id = "20260309_010802", phase = "post", control_type = NA,
         device_id = "SDSEV82D7R0087", scope_id = "Siglent SDS1104X HD",
         channel = "C3", vscale = 3.0, offset = 3.0, sample_rate = 20000,
         epoch = NA, declared = TRUE, raw_available = TRUE,
         notes = "canonical post-excitation artifact (README primary)"),
    list(analysis_id = "cap_post_1773046275",
         path = rp("csv/waveform_run_1773046275.csv"),
         run_id = "1773046275", phase = "post", control_type = NA,
         device_id = "SDSEV82D7R0087", scope_id = "Siglent SDS1104X HD",
         channel = "C3", vscale = NA, offset = NA, sample_rate = 20000,
         epoch = 1773046275, declared = NA, raw_available = TRUE,
         notes = "post-excitation harmonic run (csv waveform)"),
    list(analysis_id = "cap_run_1773047118",
         path = rp("csv/waveform_run_1773047118.csv"),
         run_id = "1773047118", phase = "unknown", control_type = NA,
         device_id = NA, scope_id = "Siglent SDS1104X HD",
         channel = NA, vscale = NA, offset = NA, sample_rate = 20000,
         epoch = 1773047118, declared = NA, raw_available = TRUE,
         notes = "low-amplitude run, phase label not documented"),
    list(analysis_id = "cap_usr_wf",
         path = rp("csv/usr_wf_data (1).csv"),
         run_id = "usr_wf", phase = "unknown", control_type = NA,
         device_id = NA, scope_id = NA, channel = NA, vscale = NA,
         offset = NA, sample_rate = 20000, epoch = NA, declared = NA,
         raw_available = TRUE,
         notes = "scope export, provenance/phase undocumented"),
    list(analysis_id = "cap_baseline_1773046194",
         path = rp("proof_runs/1773046194_baseline_peaks.csv"),
         run_id = "1773046194", phase = "baseline", control_type = NA,
         device_id = "SDSEV82D7R0087", scope_id = "Siglent SDS1104X HD",
         channel = "C3", vscale = 3.0, offset = 3.0, sample_rate = 20000,
         epoch = 1773046194, declared = TRUE, raw_available = FALSE,
         notes = "virgin diode baseline: PEAK TABLE ONLY, no raw waveform")
  )
  inv
}

# Control categories that the README declares necessary but for which NO raw
# data exists in the repository.
pd_missing_controls <- function() {
  c("sham_excited","disconnected_input","terminated_channel","dark_enclosure",
    "channel_swap","battery_isolated_supply","shielding_grounding_variation",
    "independent_device_or_day","mains_frequency_monitor")
}

# Build the full registry. Analyzable raw waveforms are loaded and profiled;
# metadata-only captures are recorded with NA analysis fields.
pd_build_registry <- function(repo_root = ".", params = pd_default_params()) {
  inv <- pd_capture_inventory(repo_root)
  flow <- list(files_discovered = length(inv), readable = 0L, malformed = 0L,
               excluded = 0L, retained = 0L, samples_loaded = 0,
               nonfinite = 0, clipped = 0, samples_used = 0)
  rows <- list(); wfs <- list()
  for (m in inv) {
    base <- data.frame(
      analysis_id = m$analysis_id, file_path = m$path,
      checksum = pd_checksum(m$path), run_id = m$run_id,
      epoch = m$epoch %||% NA, phase = m$phase,
      control_type = m$control_type %||% NA,
      device_id = m$device_id %||% NA, scope_id = m$scope_id %||% NA,
      channel = m$channel %||% NA, sample_rate = m$sample_rate %||% NA,
      declared_before_analysis = m$declared %||% NA,
      raw_available = m$raw_available, notes = m$notes,
      stringsAsFactors = FALSE)
    if (isTRUE(m$raw_available) && file.exists(m$path)) {
      wf <- try(pd_load_waveform(m$path, m$sample_rate), silent = TRUE)
      if (inherits(wf, "try-error")) {
        flow$malformed <- flow$malformed + 1L
        base$sample_count <- NA; base$status <- "MALFORMED"
        rows[[length(rows)+1]] <- base; next
      }
      flow$readable <- flow$readable + 1L; flow$retained <- flow$retained + 1L
      flow$samples_loaded <- flow$samples_loaded + length(wf$y_raw)
      flow$samples_used <- flow$samples_used + wf$n
      qc <- pd_waveform_qc(wf)
      nyq <- wf$sample_rate / 2
      base2 <- cbind(base, data.frame(
        sample_count = wf$n,
        duration_s = wf$n / wf$sample_rate,
        freq_resolution_Hz = wf$sample_rate / wf$n,
        nyquist_Hz = nyq,
        vscale = m$vscale %||% NA, offset = m$offset %||% NA,
        finite_samples = sum(is.finite(wf$y_raw)),
        missing_samples = sum(!is.finite(wf$y_raw)),
        clipped_samples = round(qc$clip_fraction * wf$n),
        quantization_step = qc$quantization_step,
        mean_voltage = qc$mean, rms_voltage = qc$rms, sd_voltage = qc$sd,
        min_voltage = qc$min, max_voltage = qc$max,
        peak_to_peak = qc$peak_to_peak,
        dc_removed = TRUE, status = "RETAINED",
        metadata_completeness = pd_meta_completeness(m)))
      wfs[[m$analysis_id]] <- wf
      rows[[length(rows)+1]] <- base2
    } else {
      flow$excluded <- flow$excluded + 1L
      base$sample_count <- NA
      base$status <- if (!file.exists(m$path)) "MISSING_FILE" else "METADATA_ONLY_NO_RAW"
      base$metadata_completeness <- pd_meta_completeness(m)
      rows[[length(rows)+1]] <- base
    }
  }
  reg <- dplyr::bind_rows(rows)
  list(registry = reg, waveforms = wfs, flow = flow,
       missing_controls = pd_missing_controls())
}

pd_meta_completeness <- function(m) {
  fields <- c("device_id","scope_id","channel","vscale","offset","sample_rate",
              "epoch","phase")
  present <- sum(vapply(fields, function(f) {
    v <- m[[f]]; !is.null(v) && !is.na(v)
  }, logical(1)))
  present / length(fields)
}
