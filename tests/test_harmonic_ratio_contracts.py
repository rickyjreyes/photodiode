"""Source-level fixtures for the harmonic ratio and log-cos Lean contracts."""
import importlib.util
import math
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


def source():
    path = Path(__file__).resolve().parents[1] / "scan_npy_logcos.py"
    spec = importlib.util.spec_from_file_location("scan_npy_logcos_contracts", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class HarmonicRatioRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.src = source()

    @staticmethod
    def peaks(frequencies):
        return pd.DataFrame({
            "frequency_Hz": frequencies,
            "amplitude": [float(i + 1) for i in range(len(frequencies))],
        })

    def test_exact_two_thirds_triplet_is_2_to_3_to_4(self):
        row = self.src.koide_triplets(self.peaks([200., 300., 400.])).iloc[0]
        self.assertAlmostEqual(row.Q_low_f1_over_f2, 2. / 3.)
        self.assertAlmostEqual(row.Q_high_f3_over_2f2, 2. / 3.)
        self.assertAlmostEqual(row.koide_error, 0.)
        self.assertEqual(row.f3, 2 * row.f1)

    def test_120_240_360_is_distinct_reference(self):
        row = self.src.koide_triplets(self.peaks([120., 240., 360.])).iloc[0]
        self.assertAlmostEqual(row.Q_low_f1_over_f2, 0.5)
        self.assertAlmostEqual(row.Q_high_f3_over_2f2, 0.75)
        self.assertGreater(row.koide_error, 0.0)
        self.assertAlmostEqual(row.integer_120_240_360_error, 0.0)

    def test_fft_frequency_bins(self):
        fs, n = 80.0, 8
        y = np.arange(n, dtype=float)
        out = self.src.compute_fft_spectrum(y, sample_rate=fs, use_cupy=False)
        np.testing.assert_allclose(out.frequency_Hz, [0., 10., 20., 30., 40.])

    def test_logcos_intercept_offset_invariance(self):
        ell = np.linspace(0.0, 2.0, 100)
        y = 0.7 + 0.23 * np.cos(3.1 * ell) - 0.11 * np.sin(3.1 * ell)
        before = self.src.fit_logcos_cpu(ell, y, 3.1)
        after = self.src.fit_logcos_cpu(ell, y + 9., 3.1)
        for k in ("chi2_null", "chi2_logcos", "delta_chi2", "a", "b"):
            self.assertAlmostEqual(before[k], after[k], places=9)
        self.assertAlmostEqual(after["C"] - before["C"], 9., places=9)


if __name__ == "__main__":
    unittest.main()
