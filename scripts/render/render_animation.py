"""ParaView (pvpython) frames of the three-bladed rotor and its helical wake, turning.

Usage: pvpython render_animation.py <rotor.npz> <frame directory> [frames]
Reads compute_rotor.py's output. Black background. Blades and hub are brushed metal; the
trailing filaments that carry at least 4 % of the reference circulation are tubes coloured by
that circulation with signed_dark (sign chosen so the tip vortices are red and the root
vortices teal; weaker filaments of the inner sheet are left out). The
rotor and wake turn together through one blade passage (120 degrees), so the loop is seamless.
With frames = 1 a single high-resolution still is written instead.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import onda_style as style  # noqa: E402
from paraview.simple import Show, Sphere, Transform, TrivialProducer, Tube  # noqa: E402
from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray  # noqa: E402
from vtkmodules.vtkCommonCore import vtkPoints  # noqa: E402
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData  # noqa: E402

src, out = Path(sys.argv[1]), Path(sys.argv[2])
frames = int(sys.argv[3]) if len(sys.argv) > 3 else 60
width, height = (3600, 1800) if frames == 1 else (1600, 800)
d = np.load(src)
R = float(d["R"])
T_RANGE = 0.3                                     # colour range of trailed circulation / gamma_ref


def polydata(points, cells, kind, arrays=()):
    pd = vtkPolyData()
    p = vtkPoints()
    p.SetData(numpy_to_vtk(np.ascontiguousarray(points.reshape(-1, 3)), deep=True))
    pd.SetPoints(p)
    ca = vtkCellArray()
    conn = np.hstack([np.full((len(cells), 1), cells.shape[1]), cells]).ravel().astype(np.int64)
    ca.SetCells(len(cells), numpy_to_vtkIdTypeArray(conn, deep=True))
    (pd.SetPolys if kind == "polys" else pd.SetLines)(ca)
    for name, values in arrays:
        a = numpy_to_vtk(np.ascontiguousarray(values, float), deep=True)
        a.SetName(name)
        pd.GetPointData().AddArray(a)
    return pd


def producer(pd):
    t = TrivialProducer()
    t.GetClientSideObject().SetOutput(pd)
    return t


strength = -d["trailed"] / d["gamma_ref"]          # tip vortices positive
keep = np.abs(strength) >= 0.04
fil = d["filaments"][keep]                        # (n_filaments, n_points, 3)
nf, npts = fil.shape[:2]
trailed = np.repeat(strength[keep], npts)
wake = producer(polydata(fil, np.arange(nf * npts).reshape(nf, npts), "lines", [("trailed", trailed)]))
blades = d["blades"]
quads = producer(polydata(blades, np.arange(len(blades) * 4).reshape(-1, 4), "polys"))
hub = Sphere(Radius=0.08 * R, ThetaResolution=64, PhiResolution=32)

view = style.setup_view(width, height)
lut = style.lut("trailed", "signed_dark", -T_RANGE, T_RANGE)
turn = [Transform(Input=s) for s in (wake, quads, hub)]
for t in turn:
    t.Transform = "Transform"
tubes = Tube(Input=turn[0])
tubes.Radius = 0.006 * R
tubes.NumberofSides = 12
style.mapped(Show(tubes, view), "trailed", lut, roughness=0.35)
lut.RescaleTransferFunction(-T_RANGE, T_RANGE)
style.metal(turn[1], view, mapping="plane", repeat=(1.0, 6.0))
style.metal(turn[2], view, mapping="sphere", repeat=(2.0, 1.0), rgb=(0.62, 0.64, 0.67))

view.CameraPosition = [-2.3 * R, 1.3 * R, 2.6 * R]
view.CameraFocalPoint = [0.85 * R, 0.0, 0.0]
view.CameraViewUp = [0, 1, 0]
view.CameraViewAngle = 33
out.mkdir(parents=True, exist_ok=True) if frames > 1 else out.parent.mkdir(parents=True, exist_ok=True)
for k in range(frames):
    angle = 120.0 * k / frames                     # the wake lags the blades, so it travels downstream
    for t in turn:
        t.Transform.Rotate = [angle, 0.0, 0.0]
    style.save(view, out if frames == 1 else out / f"frame_{k:03d}.png", width, height)
    print(k, end=" ", flush=True)
print()
