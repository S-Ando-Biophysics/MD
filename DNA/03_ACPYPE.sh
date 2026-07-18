#!/usr/bin/env bash
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh

mkdir -p 02_ACPYPE
cd 02_ACPYPE

cp ../01_AMBER/leap.prmtop .
cp ../01_AMBER/leap.inpcrd .

cd ..

cd 02_ACPYPE
conda activate acpype
acpype -p leap.prmtop -x leap.inpcrd
conda deactivate

cd ..