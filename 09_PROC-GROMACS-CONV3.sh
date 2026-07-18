#!/usr/bin/env bash
set -euo pipefail

cd 03_GROMACS/05_PROC-GROMACS

gmx trjconv -s md.tpr -f md_cluster.xtc -n index.ndx -o md_center.xtc -center -pbc mol -ur compact

# 11
# 0

cd ../..
