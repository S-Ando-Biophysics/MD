#!/usr/bin/env bash
set -euo pipefail

# Resolve paths relative to this script, regardless of the launch directory.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${SCRIPT_DIR}/INPUT"
WORK_DIR="${SCRIPT_DIR}/01_AMBER"

for file in model.pdb leap.in; do
  if [[ ! -f "${INPUT_DIR}/${file}" ]]; then
    printf 'Missing input file: %s\n' "${INPUT_DIR}/${file}" >&2
    exit 1
  fi
done

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

cp "${INPUT_DIR}/model.pdb" .

# Optional ligand parameters for leap.in files that use a ligand.
for file in LIG.mol2 LIG.frcmod; do
  if [[ -f "${INPUT_DIR}/${file}" ]]; then
    cp "${INPUT_DIR}/${file}" .
  fi
done

pdb4amber -i model.pdb -o clean.pdb --nohyd

cp "${INPUT_DIR}/leap.in" .
tleap -f leap.in
