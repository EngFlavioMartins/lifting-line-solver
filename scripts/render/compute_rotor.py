"""Solve the three-bladed rotor of examples/three_blade_rotor.py and write its geometry for ParaView.

Writes outputs/render/rotor.npz with
  blades     (nb * ncp, 4, 3)  panel corners, chord along z, span along y, axis along x
  gamma      (nb * ncp,)       bound circulation of each panel
  filaments  (nb * (ncp + 1), n, 3)  trailing filaments from the panel edges, cut at X_MAX
  trailed    (nb * (ncp + 1),) circulation trailed by each filament (the jump in bound
             circulation across the edge: strong at the tips and roots, weak in between)
  gamma_ref  pi U^2 / (nb Omega), the usual circulation scale
"""
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
from Libs.LL_ConstructObjs import construstGeometry, createRadius  # noqa: E402
from Libs.LL_MAIN import liftingLineSolver  # noqa: E402

nb, ncp, R, Rroot, TSR, dt, ndt, Uwind, a = 3, 22, 50.0, 10.0, 8.0, 0.1, 400, 10.0, 1 / 3
X_MAX = 2.0 * R                                    # length of wake to draw
Omega = TSR * Uwind / R
r, dr = createRadius(Rroot, R, ncp, "cosine")
panels, rings, sol = construstGeometry([3, 1, 1], [-14, 1, 2], r, dr, dt, Uwind * (1 - a), -Omega, nb, ncp, ndt)
sol = liftingLineSolver(sol, panels, rings, Uwind, Omega, R, "./Resources/DU95W180.csv")

gamma = np.asarray(sol.Gamma, float)
blades = np.array([[p.p1, p.p2, p.p3, p.p4] for p in panels], float)
filaments, trailed = [], []
for b in range(nb):
    g = gamma[b * ncp:(b + 1) * ncp]
    jumps = np.r_[g[0], np.diff(g), -g[-1]]          # root edge, inner edges, tip edge
    for k in range(ncp + 1):
        ring = rings[b * ncp + min(k, ncp - 1)]
        pts = np.asarray(ring.p1 if k < ncp else ring.p2, float).T
        pts = np.delete(pts, 2, axis=0)              # the wake's first point repeats the trailing edge
        filaments.append(pts[pts[:, 0] <= X_MAX])
        trailed.append(jumps[k])
n = min(len(f) for f in filaments)
out = ROOT / "outputs" / "render"
out.mkdir(parents=True, exist_ok=True)
np.savez(out / "rotor.npz", blades=blades, gamma=gamma, filaments=np.array([f[:n] for f in filaments]),
         trailed=np.array(trailed), gamma_ref=np.pi * Uwind**2 / (nb * Omega), R=R)
print(f"wrote {out / 'rotor.npz'}: {len(filaments)} filaments of {n} points")
