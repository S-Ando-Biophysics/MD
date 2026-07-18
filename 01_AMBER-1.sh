#!/usr/bin/env bash
set -euo pipefail

mkdir -p 01_AMBER
cd 01_AMBER

cp ../model.pdb .

if [[ -f ../00_LIG/LIG.mol2 ]]; then
  cp ../00_LIG/LIG.mol2 .
fi
if [[ -f ../00_LIG/LIG.frcmod ]]; then
  cp ../00_LIG/LIG.frcmod .
fi

pdb4amber -i model.pdb -o clean.pdb --nohyd

cp ../pre-leap.in .
tleap -f pre-leap.in

cd ..