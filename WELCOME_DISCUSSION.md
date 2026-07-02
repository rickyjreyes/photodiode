# Welcome to the Photodiode Experiment Discussions

This forum is for inspectable discussion of the photodiode waveform captures, FFT analysis, harmonic ladders, ratio geometry, persistence measurements, lock/relock events, controls, and possible physical interpretations.

The objective is not to defend a preferred explanation. It is to determine which observations reproduce, which controls survive, which mundane artifacts remain plausible, and what experiments would distinguish them.

## Good contributions include

- exact computational reproductions of committed captures;
- independent apparatus replications;
- failed or incomplete replications;
- dark, baseline, disconnected, dummy-load, sham, mains-monitor, and repeat-device controls;
- alternate FFT, spectral leakage, detrending, peak-detection, or null-model analyses;
- calibrated persistence, decay, and ringdown measurements;
- tests of lock/relock behavior under controlled perturbations;
- diagnoses of grounding, triggering, clipping, aliasing, saturation, or environmental pickup.

## Reporting standard

Identify the exact capture and repository commit. Preserve the raw waveform and complete acquisition metadata. Report the apparatus, wiring, grounding, shielding, optical protocol, oscilloscope settings, analysis command, parameters, dependency versions, random seed, machine-readable outputs, uncertainty, and all controls.

Separate four statements:

1. what the instrument recorded;
2. what the analysis measured in that recording;
3. what the controls showed;
4. what physical cause is being proposed.

The presence of a reproducible harmonic or log-cos feature in a stored waveform does not by itself establish its source. In particular, 60 Hz-spaced structure requires direct and repeated tests against mains and environmental pickup.

Null results, failed replications, and identified artifacts are valuable results and should remain public.
