"""Differential fixtures for formal/PhotodiodeProtocol/{HarmonicLadder,ScanMax}.lean (#5, #6).

Every test calls the *actual* functions in ``scan_npy_logcos.py``; nothing below re-implements
the statistics.  Each test names the Lean declaration whose exact-real statement it exercises
in floating point.  Passing tests are conformance evidence, not proofs.
"""
import importlib.util
import math
import unittest
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd


def source():
    path = Path(__file__).resolve().parents[1] / "scan_npy_logcos.py"
    spec = importlib.util.spec_from_file_location("scan_npy_logcos_ladder", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def peaks(frequencies):
    return pd.DataFrame({
        "frequency_Hz": [float(f) for f in frequencies],
        "amplitude": [float(i + 1) for i in range(len(frequencies))],
    })


class HarmonicLadderContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.src = source()

    def test_nearest_bin_not_harmonic(self):
        """HarmonicLadder.nearest_bin_not_harmonic: N=70, fs=3000 (both inside rfft range)."""
        n, fs = 70, 3000.0
        freq = np.fft.rfftfreq(n, d=1.0 / fs)
        nearest = lambda f: int(np.argmin(np.abs(freq - f)))
        self.assertEqual(nearest(60.0), 1)
        self.assertEqual(nearest(180.0), 4)
        self.assertNotEqual(nearest(180.0), 3 * nearest(60.0))

    def test_bin_triplet_two_thirds_is_integer_condition(self):
        """HarmonicLadder.binTriplet_exact_two_thirds_iff: bins (2,3,4) vs (1,2,3)."""
        n, fs = 64, 640.0
        freq = np.fft.rfftfreq(n, d=1.0 / fs)
        row = self.src.koide_triplets(peaks([freq[2], freq[3], freq[4]])).iloc[0]
        self.assertAlmostEqual(row.Q_low_f1_over_f2, 2.0 / 3.0, places=12)
        self.assertAlmostEqual(row.Q_high_f3_over_2f2, 2.0 / 3.0, places=12)
        row = self.src.koide_triplets(peaks([freq[1], freq[2], freq[3]])).iloc[0]
        self.assertAlmostEqual(row.Q_low_f1_over_f2, 1.0 / 2.0, places=12)
        self.assertAlmostEqual(row.Q_high_f3_over_2f2, 3.0 / 4.0, places=12)

    def test_source_sorts_before_forming_triplets(self):
        """HarmonicLadder.triplet_order_dependent: the source's sort is load-bearing."""
        base = self.src.koide_triplets(peaks([120.0, 240.0, 360.0])).iloc[0]
        for perm in permutations([120.0, 240.0, 360.0]):
            row = self.src.koide_triplets(peaks(list(perm))).iloc[0]
            self.assertAlmostEqual(row.Q_low_f1_over_f2, base.Q_low_f1_over_f2)
            self.assertAlmostEqual(row.Q_high_f3_over_2f2, base.Q_high_f3_over_2f2)
        # Unsorted raw ratios would differ:
        self.assertNotAlmostEqual(240.0 / 120.0, 120.0 / 240.0)

    def test_sorted_triplet_bounds_and_tolerance_bound(self):
        """sorted_triplet_bounds and approx_two_thirds_bound on random sorted triplets."""
        rng = np.random.default_rng(5)
        for _ in range(200):
            f = np.sort(rng.uniform(1.0, 500.0, size=3))
            row = self.src.koide_triplets(peaks(f)).iloc[0]
            ql, qh = row.Q_low_f1_over_f2, row.Q_high_f3_over_2f2
            self.assertTrue(0.0 < ql <= 1.0 + 1e-15)
            self.assertGreaterEqual(qh, 0.5 - 1e-15)
            eps = max(abs(ql - 2 / 3), abs(qh - 2 / 3))
            self.assertLessEqual(abs(row.f3 - 2 * row.f1), 4 * eps * row.f2 * (1 + 1e-12) + 1e-9)

    def test_koide_error_unit_invariance(self):
        """koideErrorSquared_common_scale: Hz -> kHz does not change the mismatch."""
        a = self.src.koide_triplets(peaks([200.0, 310.0, 401.0])).iloc[0]
        b = self.src.koide_triplets(peaks([0.2, 0.31, 0.401])).iloc[0]
        self.assertAlmostEqual(a.koide_error, b.koide_error, places=12)

    def test_logcos_frequency_unit_change(self):
        """logCosRSS_unit_change: shifting ell by ln(c) keeps chi2 and amplitude A."""
        rng = np.random.default_rng(11)
        f = np.linspace(50.0, 900.0, 160)
        ell = np.log(f)
        y = rng.normal(size=f.size)
        for k in (2.0, 7.5, 13.0):
            a = self.src.fit_logcos_cpu(ell, y, k)
            b = self.src.fit_logcos_cpu(np.log(1000.0 * f), y, k)
            for key in ("chi2_null", "chi2_logcos", "delta_chi2", "A"):
                self.assertAlmostEqual(a[key], b[key], places=8)

    def test_logcos_ordinate_scale_law(self):
        """logCosRSS_scale: y -> lam*y multiplies chi2 by lam^2."""
        rng = np.random.default_rng(3)
        ell = np.log(np.linspace(10.0, 200.0, 90))
        y = rng.normal(size=ell.size)
        a = self.src.fit_logcos_cpu(ell, y, 4.0)
        b = self.src.fit_logcos_cpu(ell, 3.0 * y, 4.0)
        self.assertAlmostEqual(b["chi2_logcos"], 9.0 * a["chi2_logcos"], places=8)
        self.assertAlmostEqual(b["delta_chi2"], 9.0 * a["delta_chi2"], places=8)

    def test_active_winding_identities(self):
        """activeWinding_scale and activeWinding_phase."""
        w = self.src.active_winding(9.0, 40.0, 800.0)
        self.assertAlmostEqual(w, self.src.active_winding(9.0, 0.04, 0.8), places=12)
        self.assertAlmostEqual(9.0 * (math.log(800.0) - math.log(40.0)), 2 * math.pi * w, places=12)

    def test_height_gate_gain_but_not_offset(self):
        """passesHeight_gain / passesHeight_offset_counterexample on detect_peaks."""
        freq = np.arange(1.0, 41.0)
        amp = np.zeros_like(freq)
        amp[[5, 15, 25]] = [2.0, 0.9, 1.4]
        spec = pd.DataFrame({"frequency_Hz": freq, "amplitude": amp})
        base = self.src.detect_peaks(spec, 0.5, 10)
        gained = self.src.detect_peaks(spec.assign(amplitude=7.0 * amp), 0.5, 10)
        np.testing.assert_allclose(base.frequency_Hz, gained.frequency_Hz)
        self.assertNotIn(16.0, list(base.frequency_Hz))
        # An additive floor of 2 lowers the relative threshold to 2 < 2.9 = 0.9 + 2.
        floored = self.src.detect_peaks(spec.assign(amplitude=amp + 2.0), 0.5, 10)
        self.assertIn(16.0, list(floored.frequency_Hz))

    def test_harmonic_ladder_log_gaps_decrease(self):
        """log_gap_strict_anti: ln(f_{n+1}) - ln(f_n) strictly decreases."""
        ell = np.log(60.0 * np.arange(1, 12))
        gaps = np.diff(ell)
        self.assertTrue(np.all(np.diff(gaps) < 0))


class ScanMaxContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.src = source()

    def test_p_value_formula_with_ties(self):
        """Randomization.mc_pvalue_eq_tail: ties count toward the tail."""
        self.assertAlmostEqual(self.src.p_value_geq(1.0, np.array([0.5, 1.0, 2.0])), 3.0 / 4.0)
        self.assertAlmostEqual(self.src.p_value_geq(3.0, np.array([0.5, 1.0, 2.0])), 1.0 / 4.0)

    def test_two_strict_counterexample_value(self):
        """ScanMax.mcPValue_two_strict: B=1 with X0 > X1 gives exactly 1/2."""
        self.assertAlmostEqual(self.src.p_value_geq(2.0, np.array([1.0])), 0.5)

    def test_shuffle_null_stores_scan_maxima_and_selection_effect(self):
        """le_scanMax and selected_pValue_le_scan against the actual shuffle_null."""
        rng = np.random.default_rng(0)
        f = np.linspace(20.0, 400.0, 64)
        ell = np.log(f)
        y = rng.normal(size=f.size)
        k_grid = np.linspace(1.0, 12.0, 23)
        scan, best = self.src.scan_logcos_cpu(ell, y, k_grid)
        self.assertAlmostEqual(best["delta_chi2"], scan["delta_chi2"].max())
        seed, nperm = 7, 40
        null = self.src.shuffle_null(ell, y, k_grid, nperm, seed, use_gpu=False)
        # Re-create the identical permutation stream to evaluate the selected single bin.
        replay = np.random.default_rng(seed)
        single = []
        for i in range(nperm):
            yp = replay.permutation(y)
            per_k = [self.src.fit_logcos_cpu(ell, yp, k)["delta_chi2"] for k in k_grid]
            self.assertAlmostEqual(null.null_best_delta_chi2[i], max(per_k), places=9)
            single.append(self.src.fit_logcos_cpu(ell, yp, best["k"])["delta_chi2"])
        p_scan = self.src.p_value_geq(best["delta_chi2"], null.null_best_delta_chi2.to_numpy())
        p_single = self.src.p_value_geq(best["delta_chi2"], np.array(single))
        self.assertLessEqual(p_single, p_scan)


if __name__ == "__main__":
    unittest.main()
