"""Domain and modular-test regressions against the source PDB proxy.

The source script imports cupy. NumPy provides the exact array API used by
the inspected dihedral and closure functions; the GPU experiment is not run.
"""
import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


def load_source():
    sys.modules["cupy"] = np
    path = Path(__file__).resolve().parents[1] / "wct_pdb_transport.py"
    spec = importlib.util.spec_from_file_location("pdb_transport_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


class DihedralContractRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.src = load_source()

    def test_translation_preserves_pseudodihedrals(self):
        rng = np.random.default_rng(625)
        p = rng.normal(size=(11, 3))
        a = self.src.dihedral_angles(p.copy())
        b = self.src.dihedral_angles(p + [2.5, -1.3, 0.8])
        np.testing.assert_allclose(a, b, atol=1e-9, rtol=1e-9)

    def test_collinear_chain_is_a_proxy_false_positive(self):
        # Every discrete cross-product normal is zero, but the regularized
        # atan2(0,0) still produces a numerical angle of zero.
        p = np.stack([np.arange(8.0), np.zeros(8), np.zeros(8)], axis=1)
        a = self.src.dihedral_angles(p)
        np.testing.assert_allclose(a, 0.0, atol=1e-12)
        self.assertTrue(bool(self.src.dihedral_transport_closes(p)))
        # A Frenet normal or geometric dihedral is undefined here.

    def test_one_sided_modulo_discards_nearby_negative_phase(self):
        T = 2 * np.pi
        d = self.src.TOL / 2
        below, above = T - d, T + d
        self.assertFalse(abs(np.mod(below, T)) < self.src.TOL)
        self.assertTrue(abs(np.mod(above, T)) < self.src.TOL)
        symmetric = lambda x: min(np.mod(x, T), T - np.mod(x, T))
        np.testing.assert_allclose(symmetric(below), symmetric(above), atol=1e-12)


if __name__ == "__main__":
    unittest.main()
