#!/usr/bin/env bash
set -euo pipefail

cd 03_GROMACS/05_PROC-GROMACS

gmx trjconv -s md.tpr -f md_center.xtc -n index.ndx -o md_fit.xtc -center -fit rot+trans

# 11
# 11
# 0

cd ../..
