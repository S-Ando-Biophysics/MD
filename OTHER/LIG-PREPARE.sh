#!/usr/bin/env bash
set -euo pipefail

NC=0

# This script belongs in MD/OTHER/.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MD_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
INPUT_DIR="${MD_DIR}/INPUT"
WORK_DIR="${SCRIPT_DIR}/LIG"

if [[ ! -f "${INPUT_DIR}/LIG.pdb" ]]; then
  printf 'Missing input file: %s\n' "${INPUT_DIR}/LIG.pdb" >&2
  exit 1
fi

mkdir -p "${INPUT_DIR}" "${WORK_DIR}"
cd "${WORK_DIR}"

cp "${INPUT_DIR}/LIG.pdb" .

antechamber -i LIG.pdb -fi pdb -o "${INPUT_DIR}/LIG.mol2" -fo mol2 -c bcc -s 2 -nc "${NC}" -rn LIG -at gaff2
parmchk2 -i "${INPUT_DIR}/LIG.mol2" -f mol2 -o "${INPUT_DIR}/LIG.frcmod" -s gaff2
