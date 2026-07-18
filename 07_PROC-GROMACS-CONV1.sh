#!/usr/bin/env bash
set -euo pipefail

cd 03_GROMACS/05_PROC-GROMACS

gmx trjconv -s md.tpr -f md.xtc -n index.ndx -o md_whole.xtc -pbc whole

# 0

cd ../..
