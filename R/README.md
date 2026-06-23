# Photodiode R statistical-audit suite

A self-contained R re-implementation, control audit, and signal-processing audit
of the photodiode post-excitation waveform analysis. It reproduces the canonical
Python `scan_npy_logcos.py` result and then stress-tests it against ordinary
explanations (mains pickup, colored noise, leakage, peak/window selection).

## Run

```bash
# fast development pass (reduced simulation counts, marked DEVELOPMENT_ONLY)
Rscript R/render_photodiode_audit.R --fast

# full-resolution audit
Rscript R/render_photodiode_audit.R --mode waveform --null-n 10000 \
  --bootstrap-n 5000 --calibration-n 10000 --injection-n 1000 --seed 20260309

# dependency check / tests
Rscript R/check_photodiode_dependencies.R
Rscript -e 'library(testthat); test_dir("tests/testthat")'
```

Outputs: `outputs_r/parity/`, `tables_r/statistical_audit/`,
`figures_r/statistical_audit/`, and the rendered report at
`reports/rendered/photodiode_statistical_audit.html`.

## Key findings (from the committed development run)

- **Parity:** R reproduces the canonical `best k`, `Δχ²`, `A`, `φ`, `χ²` to
  ~1e-15 (machine precision). R never calls the Python pipeline.
- **Log-cos nulls:** the modulation is significant only against the canonical
  *shuffle* null. It does **not** survive block-permutation, stationary-bootstrap,
  colored-noise (AR), or mains-harmonic nulls — the observed statistic sits
  *below* those null means → consistent with a finite mains ladder + colored noise.
- **Calibration:** the shuffle-null decision rule is anti-conservative on colored
  noise (observed FPR far above nominal α).
- **Q = 2/3 geometry:** an arithmetic consequence of the 60 Hz ladder; the two
  exact triplets share the 240 Hz peak and are not independent.
- **Controls:** no raw control waveforms exist; the causal "optically induced
  state" claim is `NOT_TESTABLE_WITH_AVAILABLE_DATA`.

## Modules

`photodiode_io` · `photodiode_registry` · `waveform_qc` · `spectral_estimators`
· `peak_estimation` · `harmonic_ladder` · `ratio_geometry` · `logcos_scan` ·
`logcos_nulls` · `time_frequency_analysis` · `control_comparison` ·
`persistence_analysis` · `ringdown_models` · `window_sensitivity` ·
`multiple_testing` · `null_calibration` · `injection_recovery` ·
`python_r_parity` · `build_claim_matrix` · `render_photodiode_audit` ·
`check_photodiode_dependencies`.
