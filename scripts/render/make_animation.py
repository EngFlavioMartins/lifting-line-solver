"""Animation for the README: the three-bladed rotor of examples/three_blade_rotor.py turning in its wake.

Solves the rotor with this solver (compute_rotor.py), renders one blade passage with ParaView
(render_animation.py) and writes docs/assets/animation.webp.
Run:  python scripts/render/make_animation.py   (needs ParaView; set PVPYTHON if it is not found)
"""
import subprocess
import sys
import tempfile
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from gif import crop_union, pvpython, save_gif  # noqa: E402

FRAMES = 60

subprocess.run([sys.executable, str(HERE / "compute_rotor.py")], check=True, cwd=ROOT)
with tempfile.TemporaryDirectory() as tmp:
    subprocess.run([pvpython(), str(HERE / "render_animation.py"), str(ROOT / "outputs" / "render" / "rotor.npz"),
                    tmp, str(FRAMES)], check=True)
    frames = crop_union([Image.open(p) for p in sorted(Path(tmp).glob("frame_*.png"))])
save_gif(frames, ROOT / "docs" / "assets" / "animation.webp", fps=20, width=840)
