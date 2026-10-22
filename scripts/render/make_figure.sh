#!/bin/sh
# Make the README figure in outputs/render/ (PNG, WebP and PDF).
# Needs Python with the packages in the README, ParaView's pvpython and LaTeX.
# Set PVPYTHON if pvpython is not on your PATH, e.g.
#   PVPYTHON=/Applications/ParaView-6.2.0.app/Contents/bin/pvpython sh scripts/render/make_figure.sh
set -e
cd "$(dirname "$0")/../.."
PV=${PVPYTHON:-pvpython}
R=scripts/render
O=outputs/render
mkdir -p $O
python $R/compute_rotor.py
$PV $R/render_animation.py $O/rotor.npz $O/rotor-raw.png 1
python $R/compose_rotor.py $O/rotor-raw.png $O/lifting-line-rotor
