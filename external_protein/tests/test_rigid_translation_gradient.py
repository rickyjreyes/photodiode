"""Deterministic CPU regression for the CuPy topology-flow finite-difference bug.

The module is loaded with NumPy as a CuPy-compatible stand-in, avoiding a GPU
requirement. This tests the actual source functions, not a rewritten energy.
"""
from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


def load_source():
    sys.modules["cupy"] = np
    path = Path(__file__).resolve().parents[1] / "wct_topology_flow.py"
    spec = importlib.util.spec_from_file_location("topology_flow_under_test", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class RigidTranslationRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = load_source()
        cls.points = np.random.default_rng(2405).normal(size=(9, 3))
        cls.points[:, 0] *= 1.7
        cls.points[:, 2] *= 0.6

    def test_curvature_energy_is_translation_invariant(self):
        base = self.source.curvature_energy(self.points)
        for t in ([2.5, -0.2, 1.7], [0.0, 0.0, 0.0], [-1.2, 4, 2]):
            shifted = self.source.curvature_energy(self.points + t)
            np.testing.assert_allclose(shifted, base, rtol=1e-11, atol=1e-11)

    def test_pairwise_penalty_is_translation_invariant(self):
        base = self.source.self_avoidance_energy(self.points, eps=2.0)
        shifted = self.source.self_avoidance_energy(
            self.points + np.array([3.7, 1.2, -1.3]), eps=2.0
        )
        np.testing.assert_allclose(shifted, base, rtol=1e-11, atol=1e-11)

    def test_current_bulk_gradient_is_nearly_zero(self):
        grad = self.source.energy_gradient(self.points.copy(), sa_weight=4.0, h=1e-4)
        self.assertEqual(grad.shape, self.points.shape)
        # Floating point cancellation is not identically zero.
        np.testing.assert_allclose(grad, 0.0, rtol=0.0, atol=1e-6)
        np.testing.assert_allclose(grad, np.tile(grad[0], (len(self.points), 1)), atol=0)

    def test_single_vertex_difference_is_not_bulk_translation(self):
        h = 1e-4
        base = self.points.copy()
        plus, minus = base.copy(), base.copy()
        plus[2, 0] += h
        minus[2, 0] -= h
        single = (self.source.total_energy(plus, sa_weight=4.0) -
                  self.source.total_energy(minus, sa_weight=4.0)) / (2 * h)
        self.assertGreater(abs(single), 1e-5)


if __name__ == "__main__":
    unittest.main()
