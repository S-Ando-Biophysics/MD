#!/usr/bin/env bash
set -euo pipefail

cd 03_GROMACS/04_MD

gmx mdrun -v -deffnm md -nb auto -cpi md.cpt

cd ../..

echo "=== Finished ==="