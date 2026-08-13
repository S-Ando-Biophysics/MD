#!/usr/bin/env bash
set -euo pipefail

mkdir -p 03_GROMACS/05_PROC-GROMACS

cd 03_GROMACS/05_PROC-GROMACS

cp ../04_MD/md.gro .

gmx make_ndx -f ../04_MD/md.tpr -o index.ndx

# l
# ri 1-○○
# name 11 Nucleic
# q

cd ../..
