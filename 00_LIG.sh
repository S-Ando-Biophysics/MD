#!/usr/bin/env bash
set -euo pipefail

NC=0

mkdir -p 00_LIG
cd 00_LIG

cp ../LIG.pdb .

antechamber -i LIG.pdb -fi pdb -o LIG.mol2 -fo mol2 -c bcc -s 2 -nc "${NC}" -rn LIG -at gaff2
parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod -s gaff2

cd ..
