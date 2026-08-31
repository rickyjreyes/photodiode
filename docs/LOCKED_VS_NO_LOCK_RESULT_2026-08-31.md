# Locked vs No-Lock Archived Run Contrast

**Status:** directly supported by committed run outputs, with provenance limitations stated below  
**Added:** 2026-08-31

## Question

Is the pronounced 60 Hz-spaced harmonic ladder continuously present in the archived measurement stream, or does the repository contain a later raw capture in which the lock is absent?

## Archived comparison

The repository contains a post-excitation locked run and a later retained raw run with sharply different spectral behavior.

| Run | Archived phase | RMS voltage | 60 Hz amplitude | Detected harmonics | Ladder statistic | Archived interpretation |
|---|---|---:|---:|---:|---:|---|
| `1773046275` | post | 0.0763359644 | 63.3030966 | 12 | 286.8951175 | locked harmonic run |
| `1773047118` | unknown | 0.0011999400 | 0.1081480 | 0 | 0 | later low-amplitude / no-lock run |

The epoch identifiers differ by **843 seconds (~14.1 minutes)**. The 60 Hz amplitude in the later run is approximately **585 times smaller** than in the locked `1773046275` run, while the detected harmonic count collapses from **12 to 0** and the ladder statistic from **286.9 to 0**.

Source files:

- `csv/waveform_run_1773046275.csv`
- `csv/waveform_run_1773047118.csv`
- `tables_r/statistical_audit/control_comparison.csv`
- `outputs_r/statistical_audit/capture_registry.json`

## What this establishes

The archived data support the following narrow statement:

> **The high-amplitude 60 Hz-spaced harmonic lock is not continuously present across the retained captures. A later raw capture has no detected harmonic ladder under the same audit pipeline.**

This is a stronger statement than merely observing that the amplitude changed. The state classification changes from a multi-harmonic ladder (12 detected harmonics; ladder statistic 286.9) to no detected ladder (0 harmonics; statistic 0).

This result is relevant to simple explanations in which the full locked ladder is assumed to be a continuously present, invariant feature of the oscilloscope/environment. Such an invariant-background model is not consistent with the archived run-to-run state contrast.

## What this does **not** establish

This comparison does **not** prove that all mains-related, grounding, conducted, magnetic, environmental, or instrument coupling is excluded. A coupling mechanism can itself be intermittent or condition-dependent.

The archive also does not presently establish all experimental details of run `1773047118`:

- its archived `phase` is `unknown`;
- its device identifier is not recorded in the current registry entry;
- the archive does not independently encode the exact shielding/restart chronology;
- the 0.5 s waveform does not by itself establish a delayed onset mechanism.

Accordingly, the result should be used as evidence of **state/condition dependence**, not as unique causal attribution.

## Experimental context annotation

The experimenter reports that the photodiode session used improvised electromagnetic shielding (a makeshift Faraday enclosure), that the acquisition system was restarted, and that a subsequent repeat showed a delayed response without entering the earlier locked state. This context is recorded here as a **retrospective experimental annotation added on 2026-08-31**. It is not treated as preregistered metadata and is not used to relabel the archived `phase: unknown` field without an independent contemporaneous record.

## Current interpretation

The strongest repository-supported interpretation is therefore:

1. A pronounced locked harmonic state was recorded and is reproducible computationally from the committed waveform.
2. A later retained raw capture does not contain that harmonic ladder.
3. The locked structure is therefore **condition-dependent rather than continuously present across archived runs**.
4. Its physical cause remains unresolved.
5. Mains-related coupling remains an unexcluded alternative, but the 60 Hz spacing alone does not establish mains as the source.

## Next decisive controls

Matched raw captures under identical acquisition settings remain the clean path to causal discrimination: disconnected/terminated input, sham excitation, battery isolation, channel swap, controlled shielding/grounding changes, mains phase/frequency monitoring, and independent device/day replication.
