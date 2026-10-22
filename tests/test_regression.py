"""Regression test: the vectorised solver reproduces the original (2023) loop-based solver.

The reference was produced by the original code for the 3-blade rotor of Examples.ipynb
(22 panels per blade, 10 wake steps), on macOS. Results agree to ~1e-15 on the machine that
made the reference; across platforms (other BLAS and libm) they differ by up to ~3e-7 relative,
far below the solver's own convergence tolerance (0.5 % change in circulation), so the test
allows 1e-6. A change to the algorithm shows up well above that.
Run with:  python -m unittest discover tests
"""
import os
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


class RegressionTest(unittest.TestCase):
    def test_three_blade_rotor_matches_original(self):
        cwd = os.getcwd()
        os.chdir(ROOT)
        try:
            from Libs.LL_ConstructObjs import construstGeometry, createRadius
            from Libs.LL_MAIN import liftingLineSolver
            r, dr = createRadius(10, 50, 22, "cosine")
            Omega = 8 * 10 / 50
            panels, rings, sol = construstGeometry([3, 1, 1], [-14, 1, 2], r, dr, 0.1, 10, Omega, 3, 22, 10)
            sol = liftingLineSolver(sol, panels, rings, 10, Omega, 50, "./Resources/DU95W180.csv")
        finally:
            os.chdir(cwd)
        ref = np.load(ROOT / "tests" / "data" / "reference_3blades_ndt10.npz")
        for key in ref.files:
            np.testing.assert_allclose(getattr(sol, key), ref[key], rtol=1e-6, atol=1e-9, err_msg=key)


if __name__ == "__main__":
    unittest.main()
