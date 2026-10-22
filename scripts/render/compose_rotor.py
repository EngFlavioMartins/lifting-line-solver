"""Add the LaTeX colour bar to the rotor render.

Usage: python compose_rotor.py <raw png> <output stem>
The bar is the circulation trailed by each filament in units of pi U^2 / (N_b Omega), on the
range of render_animation.py, with the tip vortices positive.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from overlay import compose  # noqa: E402

raw, stem = Path(sys.argv[1]), Path(sys.argv[2])
compose([(raw, None)], stem, bar=dict(name="signed_dark", vmin=-0.3, vmax=0.3, title=r"$\Gamma/\Gamma_0$",
                                     ticks=[(0.0, "$-$0.3"), (0.5, "0"), (1.0, "0.3")]))
