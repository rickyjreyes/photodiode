# Harmonic/statistical semantics (#5)

The new \`HarmonicRatios.lean\` module extends the existing
\`PhotodiodeProtocol.Spectral\` exact-real model. It represents:

* \`binFreq\`: \`k*fs/N\`, with multiple/step identities;
* \`lowQ\`: \`f1/f2\` and \`highQ\`: \`f3/(2*f2)\` for ordered triples;
* the squared Koide-style two-thirds mismatch;
* \`logCosModel\` and \`logCosRSS\` intercept-shift invariance.

Important: the source's independent 120/240/360 reference triplet gives
\`Q_low=1/2\`, \`Q_high=3/4\`; an exact two-thirds match instead has
\`f1:f2:f3=2:3:4\`. This is **not** a claim of experimental Koide
mass physics or a photonic origin. The actual recorded spectrum and its
frequency-window/peak detection remain subject to thresholding and
floating-point effects. In particular division by zero returns zero
in Lean but is not interpreted as a valid physical ratio.

\`tests/test_harmonic_ratio_contracts.py\` runs deterministic fixtures
against the actual source implementation (not a rewritten formula).

Build: \`cd formal && lake exe cache get && lake build\`.
Python regression:
\`python -m unittest discover -s tests -p 'test_harmonic_ratio_contracts.py' -v\`.

The code and theorem correspondence are **subject to a real Lean build**;
a failed or step-less Actions run cannot certify these theorems.
