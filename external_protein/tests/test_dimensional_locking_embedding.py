"""Source-level falsification fixture for dimension-specific unsigned turning.

Tests actual wct_dimensional_locking.py with NumPy substituted for CuPy.
The square is closed and nondegenerate. Its 3D and 4D scores must agree.
"""
import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


def load_source():
    sys.modules["cupy"] = np
    path = Path(__file__).resolve().parents[1] / "wct_dimensional_locking.py"
    spec = importlib.util.spec_from_file_location("dimensional_locking_under_test", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class DimensionEmbeddingRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = load_source()

    def test_closed_square_has_identical_3d_and_4d_turning(self):
        square = np.array([[0., 0., 0.], [1., 0., 0.],
                           [1., 1., 0.], [0., 1., 0.]])
        lifted = np.pad(square, ((0, 0), (0, 1)))
        a3 = self.source.transport_angle(square)
        a4 = self.source.transport_angle(lifted)
        np.testing.assert_allclose(a3, 2 * np.pi, atol=1e-12)
        np.testing.assert_allclose(a4, a3, atol=1e-12)
        self.assertEqual(bool(self.source.transport_closes(square)),
                         bool(self.source.transport_closes(lifted)))

    def test_random_closed_polygon_embedding_preserves_score(self):
        rng = np.random.default_rng(340)
        points3 = rng.normal(size=(32, 3))
        points4 = np.pad(points3, ((0, 0), (0, 1)))
        np.testing.assert_allclose(self.source.transport_angle(points3),
                                   self.source.transport_angle(points4),
                                   atol=1e-12, rtol=1e-12)

    def test_one_sided_modulo_is_asymmetric_near_two_pi(self):
        tau = 2 * np.pi
        delta = self.source.TOL * 0.5
        near_below = tau - delta
        near_above = tau + delta
        self.assertFalse(abs(np.mod(near_below, tau)) < self.source.TOL)
        self.assertTrue(abs(np.mod(near_above, tau)) < self.source.TOL)
        # A symmetric circular distance treats both as equally close to 0.
        dist = lambda x: min(np.mod(x, tau), tau - np.mod(x, tau))
        np.testing.assert_allclose(dist(near_below), dist(near_above), atol=1e-12)


if __name__ == "__main__":
    unittest.main()
