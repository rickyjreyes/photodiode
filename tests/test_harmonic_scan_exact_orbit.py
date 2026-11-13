"""Additional source-level audit fixtures for photodiode #5/#6.

Finite-orbit permutation enumeration tests scan-max selection using the real
released scanner. This validates a *synthetic exchangeable orbit only* and
does not establish exchangeability of correlated experimental FFT bins.
"""
import importlib.util
import itertools
import math
from pathlib import Path
import unittest

import numpy as np


def actual_scanner():
    source = Path(__file__).resolve().parents[1] / "scan_npy_logcos.py"
    spec = importlib.util.spec_from_file_location("photodiode_scan_orbit", source)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class FiniteOrbitSourceChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.src = actual_scanner()

    def test_scanmax_exact_orbit_superuniform(self):
        """Every permutation's scan maximum is calibrated over the *whole* orbit."""
        f = np.array([21., 46., 95., 193., 382.])
        ell = np.log(f)
        y = np.array([-2., -1., 0.5, 1.8, 3.7])
        kgrid = np.array([1.5, 3.1, 7.0])
        stats = []
        for perm in itertools.permutations(y.tolist()):
            scan, best = self.src.scan_logcos_cpu(ell, np.array(perm), kgrid)
            self.assertAlmostEqual(float(best["delta_chi2"]),
                                   float(scan["delta_chi2"].max()), places=9)
            stats.append(float(best["delta_chi2"]))
        stats = np.asarray(stats)
        pvalues = np.array([np.mean(stats >= s - 1e-10) for s in stats])
        for alpha in (0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0):
            self.assertLessEqual(np.mean(pvalues <= alpha), alpha + 1e-12)

    def test_preselected_bin_is_not_familywise_calibration(self):
        """Single-bin calibration must not be promoted to scan significance."""
        scores = np.array([0., 1., 2., 3.])
        selected = scores[2]
        scan = float(scores.max())
        self.assertLess(selected, scan)
        # An exactly finite null: two marginally valid coin p-values.
        cases = [(0.5, 0.5), (0.5, 1.0), (1.0, 0.5), (1.0, 1.0)]
        self.assertEqual(sum(min(a, b) <= 0.5 for a, b in cases), 3)
        self.assertEqual(len(cases), 4)

    def test_nonexchangeable_two_draw_null_is_invalid(self):
        """If observed is *always* larger, a one-draw randomization p is always 1/2."""
        for observed, null in ((2.0, 1.0), (4.0, -3.0), (100.0, 0.0)):
            self.assertEqual(self.src.p_value_geq(observed, np.array([null])), 0.5)

    def test_nyquist_bin_index_and_harmonic_ladder(self):
        """All represented harmonics must stay within the one-sided FFT range."""
        for n in (33, 34, 69, 70, 128):
            fs = 3000.0
            bins = np.fft.rfftfreq(n, 1/fs)
            self.assertEqual(len(bins), n//2+1)
            self.assertTrue(np.all(np.diff(bins) > 0))
            self.assertLessEqual(bins[-1], fs/2 + 1e-10)
            for h in (1, 2, 3):
                target = h*60.0
                if target <= bins[-1]:
                    ix = int(np.argmin(abs(bins-target)))
                    self.assertLessEqual(abs(bins[ix]-target),
                                         fs/n/2 + 1e-9)


if __name__ == "__main__":
    unittest.main()
