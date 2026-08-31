# ==============================================================================
# build_claim_matrix.R
# Assemble the final claim-status matrix. Verdicts are conservative and tied to
# what the available data can actually support.
# ==============================================================================

PD_VERDICTS <- c("DIRECTLY_MEASURED","COMPUTATIONALLY_REPRODUCED",
                 "CONSISTENT_WITH_MAINS","POST_BASELINE_INCREASE",
                 "CONTROL_SEPARATED","NOT_CONTROL_SEPARATED",
                 "ARITHMETIC_CONSEQUENCE","ROBUST_STATISTICAL_FEATURE",
                 "INCONCLUSIVE","FAILED","NOT_TESTABLE_WITH_AVAILABLE_DATA")

pd_claim_row <- function(claim, required, observed, uncertainty, baseline,
                         control, global_p, calibration, injection, verdict,
                         limitation, output_ref) {
  data.frame(claim = claim, required_evidence = required, observed_result = observed,
             uncertainty = uncertainty, baseline_result = baseline,
             control_result = control, global_p = global_p,
             calibration_status = calibration, injection_sensitivity = injection,
             verdict = verdict, limitation = limitation,
             output_reference = output_ref, stringsAsFactors = FALSE)
}

pd_build_claim_matrix <- function(res) {
  p <- res
  rows <- list()
  rows[[1]] <- pd_claim_row(
    "Captured waveform contains a 60 Hz-spaced harmonic ladder",
    "Discrete peaks at integer multiples of 60 Hz in the post capture",
    sprintf("%d harmonics detected; peaks at %s Hz",
            p$ladder_n %||% NA, p$ladder_freqs %||% "60..600"),
    "FFT bin spacing = 2 Hz (0.5 s capture)",
    "Baseline peak table also shows 6 peaks coincident with 60 Hz harmonics",
    "No matched raw control waveforms available",
    NA, "n/a", sprintf("ladder detection power high for amp>=%s", p$inj_amp %||% "0.5"),
    "DIRECTLY_MEASURED",
    "Native 2 Hz resolution; sub-bin estimates required for exact frequency",
    "tables_r/.../harmonic_ladder.csv; figures fig06")
  rows[[2]] <- pd_claim_row(
    "Canonical Python log-cos result is reproducible in R",
    "best k, delta chi2, A, phi within tolerance of summary.json",
    sprintf("best k=%.6f, delta_chi2=%.6f, A=%.6f", p$best_k, p$best_delta, p$best_A),
    sprintf("max rel diff = %.2e", p$parity_maxreldiff %||% NA),
    "n/a", "n/a", NA, "n/a", "n/a",
    if (isTRUE(p$parity_pass)) "COMPUTATIONALLY_REPRODUCED" else "FAILED",
    "Deterministic part exact; shuffle null compared by distribution only",
    "outputs_r/parity/; tables_r/.../python_r_parity.csv")
  rows[[3]] <- pd_claim_row(
    "Exact 2:3:4 frequency geometry occurs in the detected peaks",
    "Q_low=Q_high=2/3 triplets among detected peaks",
    sprintf("%d exact Q=2/3 triplets", p$q_exact %||% NA),
    "Exact by construction on integer harmonics",
    "n/a", "n/a", NA, "n/a", "n/a",
    "ARITHMETIC_CONSEQUENCE",
    "Both triplets lie on the 60 Hz ladder and share the 240 Hz peak; they are not statistically independent",
    "tables_r/.../ratio_geometry.csv, triplet_overlap.csv")
  rows[[4]] <- pd_claim_row(
    "Log-cos modulation exceeds the canonical shuffle null",
    "Observed scan-max delta chi2 above shuffle scan-max distribution",
    sprintf("p_shuffle=%.4g (canonical)", p$p_shuffle %||% NA),
    sprintf("resolution floor 1/%d", (p$null_n %||% 1000)+1),
    "n/a", "n/a", p$p_shuffle %||% NA, p$calib_status %||% "see fig18", "n/a",
    "COMPUTATIONALLY_REPRODUCED",
    paste("Reproduces the canonical shuffle result; the shuffle null does not",
          "preserve all colored-noise or frequency-local correlations, so this",
          "exceedance is not by itself a source-attribution test."),
    "outputs_r/parity/logcos_null.csv; fig11")
  rows[[5]] <- pd_claim_row(
    "Log-cos modulation under structured-noise and harmonic-ladder nulls",
    "Compare the observed statistic with AR colored-noise and synthetic harmonic-ladder null distributions",
    p$strong_null_summary %||% "see logcos_null_comparison.csv",
    "Development-resolution nulls",
    "n/a", "n/a", p$p_mains %||% NA, "see fig18",
    "n/a",
    if (identical(p$logcos_strong_verdict, "ROBUST_STATISTICAL_FEATURE")) "ROBUST_STATISTICAL_FEATURE" else "INCONCLUSIVE",
    "These nulls test statistical compatibility with model classes; they do not identify the physical source of the measured waveform",
    "tables_r/.../logcos_null_comparison.csv; fig12")
  rows[[6]] <- pd_claim_row(
    "Harmonic ladder is an optically induced device state",
    "Matched control separation: baseline/sham/disconnected/dark differ from post under the same acquisition conditions",
    "Not fully evaluable: baseline is peak-table only and several matched raw controls are absent",
    "n/a", "Baseline peak table contains 60 Hz-coincident frequencies",
    "MISSING matched raw controls: sham, disconnected, terminated, dark, channel-swap, battery",
    NA, "n/a", "n/a",
    "NOT_TESTABLE_WITH_AVAILABLE_DATA",
    "The current archive does not establish a unique causal mechanism",
    "tables_r/.../control_comparison.csv; fig17")
  rows[[7]] <- pd_claim_row(
    "Structure persists beyond the recorded capture",
    "Multiple sequential timestamped post captures showing decay",
    "Single 0.5 s canonical post waveform; within-capture stability only",
    "n/a", "n/a", "n/a", NA, "n/a", "n/a",
    "NOT_TESTABLE_WITH_AVAILABLE_DATA",
    "No long-duration claim can be made from the canonical 0.5 s capture alone",
    "fig16; persistence_results.csv")
  rows[[8]] <- pd_claim_row(
    "A resonator quality factor (Q_resonator = pi f tau) is established",
    "Identified oscillation frequency + directly measured decay",
    "Not measured; no canonical ringdown time series",
    "n/a", "n/a", "n/a", NA, "n/a", "n/a",
    "NOT_TESTABLE_WITH_AVAILABLE_DATA",
    "Q_ratio (2/3 geometry) must not be confused with a resonator Q",
    "ringdown_models.csv")
  rows[[9]] <- pd_claim_row(
    "Physical origin of the 60 Hz-spaced ladder",
    "Matched controls that distinguish device-state dynamics from conducted, magnetic, environmental, and instrument coupling",
    "The frequency spacing is compatible with 60 Hz-related coupling, but frequency coincidence alone does not identify the source",
    "n/a", "Baseline peak table also contains 60 Hz-coincident frequencies",
    "Existing observations are condition-dependent, but decisive matched raw controls are incomplete",
    NA, "n/a", "n/a",
    "INCONCLUSIVE",
    "Mains-related coupling remains an unexcluded alternative; the available data neither establish nor uniquely favor it as the origin",
    "fig06, fig07; harmonic_ladder.csv; control_comparison.csv")
  do.call(rbind, rows)
}
