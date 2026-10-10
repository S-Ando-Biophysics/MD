#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

PROC_DIR="04_PROCESSING-AMBER"

PDB="${PROC_DIR}/md_target.pdb"
INDEX="${PROC_DIR}/index.ndx"

# ------------------------------------------------------------
# Check input file
# ------------------------------------------------------------

if [ ! -f "$PDB" ]; then
    echo "ERROR: PDB file not found:"
    echo "  $PDB"
    exit 1
fi

# ------------------------------------------------------------
# Create index
# ------------------------------------------------------------

echo
echo "========================================"
echo " Create index "
echo "========================================"
echo

echo "Input:"
echo "  PDB   : $PDB"
echo

echo "Output:"
echo "  INDEX : $INDEX"
echo

gmx make_ndx \
    -f "$PDB" \
    -o "$INDEX"

echo
echo "========================================"
echo " Index creation completed "
echo "========================================"
echo

echo "Created:"
echo "  $INDEX"
echo