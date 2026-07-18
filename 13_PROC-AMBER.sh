#!/usr/bin/env bash
set -euo pipefail

FIT_RES="1-○○"

mkdir -p 03_GROMACS/05_PROC-AMBER

cd 03_GROMACS/05_PROC-AMBER

cp ../04_MD/md.xtc .
cp ../../01_AMBER/leap.prmtop .

cat > trajfit.in << EOF
parm leap.prmtop
trajin md.xtc

unwrap :${FIT_RES}
center :${FIT_RES} mass origin
autoimage

rms first :${FIT_RES}&!@H= out rmsd_fit_nosol.dat

strip :WAT,SOL,Na+,NA+,NA,Cl-,CL-,CL

trajout md_fit_nosol.nc netcdf
trajout md_fit_nosol.xtc xtc
trajout md_nosol.pdb pdb onlyframes 1

run
quit
EOF

cpptraj -i trajfit.in

cd ../..