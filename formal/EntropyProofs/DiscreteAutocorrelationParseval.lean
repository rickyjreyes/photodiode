import Mathlib

/-!
# Exact 2×2 periodic real-grid DFT, Parseval, and Wiener–Khinchin

Source: src/wct_entropy/diagnostics.py:centered_autocorrelation.

For a real 2×2 array with row-major entries a, b, c, d, numpy.fft.fft2
has four *real* Fourier values (the phases are ±1). The default inverse
FFT divides by 4. The following algebra uses exactly these conventions.

This is a nontrivial fully expanded finite-grid witness, *not* a theorem
for arbitrary rectangular dimensions. The full n×m theorem is still open.
-/

namespace EntropyProofs.DiscreteAutocorrelationParseval

def fft00 (a b c d : ℝ) : ℝ := a + b + c + d
def fft01 (a b c d : ℝ) : ℝ := a - b + c - d
def fft10 (a b c d : ℝ) : ℝ := a + b - c - d
def fft11 (a b c d : ℝ) : ℝ := a - b - c + d

/-- The unnormalized forward DFT power at each frequency. -/
def spectral00 (a b c d : ℝ) : ℝ := (fft00 a b c d) ^ 2
def spectral01 (a b c d : ℝ) : ℝ := (fft01 a b c d) ^ 2
def spectral10 (a b c d : ℝ) : ℝ := (fft10 a b c d) ^ 2
def spectral11 (a b c d : ℝ) : ℝ := (fft11 a b c d) ^ 2

/-- Discrete field L² energy; the array is real, so |x|²=x². -/
def energy (a b c d : ℝ) : ℝ := a ^ 2 + b ^ 2 + c ^ 2 + d ^ 2

/-- Inverse FFT of the four power-spectrum coefficients at lag (0,0). -/
def invPower00 (a b c d : ℝ) : ℝ :=
  (spectral00 a b c d + spectral01 a b c d +
    spectral10 a b c d + spectral11 a b c d) / 4

def invPower01 (a b c d : ℝ) : ℝ :=
  (spectral00 a b c d - spectral01 a b c d +
    spectral10 a b c d - spectral11 a b c d) / 4

def invPower10 (a b c d : ℝ) : ℝ :=
  (spectral00 a b c d + spectral01 a b c d -
    spectral10 a b c d - spectral11 a b c d) / 4

def invPower11 (a b c d : ℝ) : ℝ :=
  (spectral00 a b c d - spectral01 a b c d -
    spectral10 a b c d + spectral11 a b c d) / 4

/-- Four-point Parseval identity, unnormalized forward transform. -/
theorem parseval_2x2 (a b c d : ℝ) :
    spectral00 a b c d + spectral01 a b c d +
        spectral10 a b c d + spectral11 a b c d =
      4 * energy a b c d := by
  dsimp [spectral00, spectral01, spectral10, spectral11,
    fft00, fft01, fft10, fft11, energy]
  ring

/-- Discrete Wiener–Khinchin identity at zero lag. -/
theorem wiener_khinchin_zero_lag (a b c d : ℝ) :
    invPower00 a b c d = energy a b c d := by
  dsimp [invPower00, spectral00, spectral01, spectral10,
    spectral11, fft00, fft01, fft10, fft11, energy]
  ring

/-- Lag (0,1): one spatial step in the column direction on the periodic grid. -/
theorem wiener_khinchin_column_lag (a b c d : ℝ) :
    invPower01 a b c d = 2 * (a * b + c * d) := by
  dsimp [invPower01, spectral00, spectral01, spectral10,
    spectral11, fft00, fft01, fft10, fft11]
  ring

/-- Lag (1,0): one spatial step in the row direction. -/
theorem wiener_khinchin_row_lag (a b c d : ℝ) :
    invPower10 a b c d = 2 * (a * c + b * d) := by
  dsimp [invPower10, spectral00, spectral01, spectral10,
    spectral11, fft00, fft01, fft10, fft11]
  ring

/-- Lag (1,1): diagonal wrap-around on the periodic grid. -/
theorem wiener_khinchin_diagonal_lag (a b c d : ℝ) :
    invPower11 a b c d = 2 * (a * d + b * c) := by
  dsimp [invPower11, spectral00, spectral01, spectral10,
    spectral11, fft00, fft01, fft10, fft11]
  ring

/-- Field intensity is nonnegative for every real four-point configuration. -/
theorem energy_nonneg (a b c d : ℝ) :
    0 ≤ energy a b c d := by
  dsimp [energy]
  positivity

/-- Zero energy means every real field entry is zero. -/
theorem energy_eq_zero_iff (a b c d : ℝ) :
    energy a b c d = 0 ↔
      a = 0 ∧ b = 0 ∧ c = 0 ∧ d = 0 := by
  constructor
  · intro h
    dsimp [energy] at h
    have ha : a ^ 2 = 0 := by
      nlinarith [sq_nonneg b, sq_nonneg c, sq_nonneg d]
    have hb : b ^ 2 = 0 := by
      nlinarith [sq_nonneg a, sq_nonneg c, sq_nonneg d]
    have hc : c ^ 2 = 0 := by
      nlinarith [sq_nonneg a, sq_nonneg b, sq_nonneg d]
    have hd : d ^ 2 = 0 := by
      nlinarith [sq_nonneg a, sq_nonneg b, sq_nonneg c]
    constructor
    · nlinarith [sq_nonneg a]
    constructor
    · nlinarith [sq_nonneg b]
    constructor
    · nlinarith [sq_nonneg c]
    · nlinarith [sq_nonneg d]
  · rintro ⟨rfl, rfl, rfl, rfl⟩
    simp [energy]

end EntropyProofs.DiscreteAutocorrelationParseval
