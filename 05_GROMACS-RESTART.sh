#!/usr/bin/env bash
set -euo pipefail

cd 03_GROMACS/04_MD

gmx mdrun -v -deffnm md -nb auto -nbfe auto -pme auto -pmefft auto -bonded auto -update auto -cpi md.cpt

cd ../..

echo "=== Finished ==="
