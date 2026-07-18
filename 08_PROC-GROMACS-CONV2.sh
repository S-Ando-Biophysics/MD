#!/usr/bin/env bash
set -euo pipefail

cd 03_GROMACS/05_PROC-GROMACS

gmx trjconv -s md.tpr -f md_whole.xtc -n index.ndx -o md_cluster.xtc -pbc cluster

# 11
# 0

cd ../..
