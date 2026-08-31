# Diode-Dependent Harmonic State Observation

**Status:** retrospective experimental observation; dedicated matched raw control sequence not yet archived  
**Added:** 2026-08-31

## Observation

The experimenter reports that the pronounced harmonic state can be turned **on and off with the photodiode**, while the surrounding oscilloscope/acquisition environment is otherwise held fixed.

This is an important causal observation because it is different from merely observing 60 Hz in the laboratory environment. If the full harmonic ladder appears when the diode is in the signal path/state associated with the effect and disappears when the diode condition is removed or disabled, then the measured state is **diode-dependent** rather than an invariant feature of the oscilloscope or ambient background alone.

This statement is currently recorded as a **retrospective experimental annotation**. The present archive does not yet contain a dedicated, labeled ABAB or randomized ON/OFF raw-control sequence sufficient to quantify the switching probability directly from committed waveforms.

## Relation to the archived locked/no-lock result

The repository independently contains a directly analyzable run-to-run state contrast:

| Run | Archived phase | 60 Hz amplitude | Detected harmonics | Ladder statistic | State |
|---|---|---:|---:|---:|---|
| `1773046275` | post | 63.3030966 | 12 | 286.8951175 | **LOCKED** |
| `1773047118` | unknown | 0.1081480 | 0 | 0 | **NO LOCK** |

Those runs are separated by 843 seconds (~14.1 minutes). The later retained waveform demonstrates that the full locked ladder is not continuously present in the measurement stream.

The experimenter's diode-control observation supplies additional experimental context for that state dependence, but the repository does **not** retroactively relabel run `1773047118` as a formal diode-OFF control without contemporaneous metadata tying that run to the specific switch operation.

## What diode dependence would rule out

A reproducible diode ON/OFF association would strongly disfavor a simple model in which the full harmonic ladder is merely an invariant background generated independently of the diode.

It would **not**, by itself, establish an exotic mechanism or WCT. Conventional explanations would still include:

- nonlinear rectification by the photodiode junction;
- diode-dependent impedance or capacitance changing coupling into the scope input;
- carrier trapping, leakage, or surface-state changes after illumination;
- diode-mediated pickup or mixing of a 60 Hz environmental field;
- grounding or cable-path changes coupled to the diode condition.

The scientifically important distinction is therefore:

> **Ambient 60 Hz may provide a drive or reference frequency, while the diode may determine whether the high-amplitude multi-harmonic locked state exists.**

That is different from attributing the entire observation to continuously present mains pickup.

## Decisive matched control

The next high-value test is a predeclared diode-state switching sequence with no other acquisition changes.

Recommended minimum sequence:

| Trial | Diode condition | Analysis label |
|---|---|---|
| A1 | effect-enabled / ON condition | blinded |
| B1 | effect-disabled / OFF condition | blinded |
| A2 | ON | blinded |
| B2 | OFF | blinded |
| A3 | ON | blinded |
| B3 | OFF | blinded |

A stronger test uses 10-20 randomized ON/OFF trials with the waveform classifier blinded to condition.

For every trial, keep fixed as far as practical:

- scope channel and vertical scale;
- sample rate and capture duration;
- cable placement and grounding;
- shielding geometry;
- optical geometry except for the declared diode-state manipulation;
- analysis window and thresholds.

Commit the raw waveform, timestamp, condition label, scope metadata, and analysis output for every trial, including failures and ambiguous trials.

## Primary statistic

For a blinded ON/OFF control, define a lock classifier before inspecting the condition labels. A simple version can use:

- detected harmonic count;
- ladder statistic;
- 60 Hz amplitude;
- optionally the predeclared log-cos statistic as a secondary diagnostic.

The key causal quantity is the separation between conditions, for example:

`P(lock | diode ON)` versus `P(lock | diode OFF)`.

A near-deterministic ON/OFF separation reproduced across repeated trials would be substantially stronger evidence of diode dependence than a single high-significance spectral statistic from one locked capture.

## Current claim level

**Directly supported by committed data:** the full harmonic lock is not continuously present across retained captures; a later raw run contains no detected ladder under the same audit pipeline.

**Experimenter-reported:** the harmonic state can be switched on and off with the diode under otherwise unchanged measurement conditions.

**Not yet established from a dedicated archived control series:** the quantitative causal probability of lock given diode ON versus diode OFF, or the physical mechanism producing the diode-dependent state.

The physical origin therefore remains unresolved. Mains-related coupling remains an unexcluded possible drive/coupling mechanism, but the diode-control observation motivates treating **diode-dependent state formation or nonlinear transduction** as a distinct hypothesis rather than collapsing the result into a generic mains-pickup attribution.