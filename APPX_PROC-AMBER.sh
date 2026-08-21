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

rms first :${FIT_RES}&!@H= out rmsd_target.dat

strip :WAT,SOL,Na+,NA+,NA,K+,K,Cl-,CL-,CL

trajout md_fit_target.nc netcdf
trajout md_fit_target.xtc xtc
trajout md_target.pdb pdb onlyframes 1

run
quit
EOF

cpptraj -i trajfit.in

cd ../..
