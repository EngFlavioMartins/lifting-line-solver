"""Three-bladed wind-turbine rotor with the lifting-line method (the main case of Examples.ipynb).

Rotor radius 50 m, root 10 m, tip-speed ratio 8, DU95W180 airfoil, 22 panels per blade and
a 400-step helical wake convected at U(1 - a) with a = 1/3, as in Examples.ipynb. Prints the thrust and power coefficients and saves the spanwise
circulation and induction to outputs/three_blade_rotor.png.
"""
import os
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
from Libs.LL_ConstructObjs import construstGeometry, createRadius  # noqa: E402
from Libs.LL_MAIN import liftingLineSolver  # noqa: E402

nb, ncp, R, Rroot, TSR, dt, ndt, Uwind, a = 3, 22, 50.0, 10.0, 8.0, 0.1, 400, 10.0, 1 / 3
Omega = TSR * Uwind / R
r, dr = createRadius(Rroot, R, ncp, "cosine")

t0 = time.time()
# the wake is convected at U(1 - a) and turns opposite to the rotor (-Omega)
panels, rings, sol = construstGeometry([3, 1, 1], [-14, 1, 2], r, dr, dt, Uwind * (1 - a), -Omega, nb, ncp, ndt)
sol = liftingLineSolver(sol, panels, rings, Uwind, Omega, R, "./Resources/DU95W180.csv")
elapsed = time.time() - t0

# rotor coefficients from the blade loads (one blade's spanwise loads x number of blades)
dr_all = np.tile(dr, nb)
r_all = np.tile(r, nb)
thrust = np.sum(sol.Fnorm * dr_all)
power = np.sum(sol.Ftan * dr_all * r_all) * Omega
q = 0.5 * Uwind**2 * np.pi * R**2
print(f"Solved in {elapsed:.2f} s: CT = {thrust / q:.3f}, CP = {power / (q * Uwind):.3f}")

out = ROOT / "outputs"
out.mkdir(exist_ok=True)
fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
s = slice(0, ncp)                                  # first blade
axs[0].plot(sol.r_R[s], sol.Gamma[s] / (np.pi * Uwind**2 / (nb * Omega)), "o-", ms=3)
axs[0].set(xlabel="$r/R$", ylabel=r"$\Gamma \,/\, (\pi U_\infty^2 / N_b\Omega)$", title="Circulation")
axs[1].plot(sol.r_R[s], sol.a[s], "o-", ms=3, label="axial, $a$")
axs[1].plot(sol.r_R[s], sol.aline[s], "s-", ms=3, label="azimuthal, $a'$")
axs[1].set(xlabel="$r/R$", ylabel="induction", title="Induction factors")
axs[1].legend(frameon=False)
fig.savefig(out / "three_blade_rotor.png", dpi=150)
print(f"Figure: {out / 'three_blade_rotor.png'}")
