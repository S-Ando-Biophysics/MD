#!/usr/bin/env bash
set -euo pipefail

cd 03_GROMACS/05_PROC-GROMACS

gmx trjconv -s md.tpr -f md_fit.xtc -n index.ndx -o md_fit_nosol.xtc

awk '
function trim(s) {
    gsub(/^ +| +$/, "", s)
    return s
}

function norm_resname(r) {
    r = trim(r)

    if (r ~ /[53]$/) {
        return substr(r, 1, length(r)-1)
    }

    return r
}

function is_nucleic(r) {
    return (r == "A"  || r == "C"  || r == "G"  || r == "U"  || r == "T"  || \
            r == "DA" || r == "DC" || r == "DG" || r == "DU" || r == "DT")
}

function is_ligand(r) {
    return (r == "LIG")
}

(/^ATOM  / || /^HETATM/) {
    res_raw = substr($0, 18, 3)
    res = trim(res_raw)
    newres = norm_resname(res)

    if (is_nucleic(newres) || is_ligand(newres)) {
        print
    }

    next
}

/^TER/ {
    next
}

/^END/ {
    print
    next
}
' md.pdb > md_nosol.pdb

cd ../..