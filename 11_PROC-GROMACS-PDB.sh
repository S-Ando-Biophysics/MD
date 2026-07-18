#!/usr/bin/env bash
set -euo pipefail

cd 03_GROMACS/05_PROC-GROMACS

gmx trjconv -s md.tpr -f md_fit.xtc -n index.ndx -o md.pdb -dump 0

awk '
BEGIN {
    chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
}

function trim(s) {
    gsub(/^ +| +$/, "", s)
    return s
}

function pad80(s) {
    while (length(s) < 80) s = s " "
    return s
}

function chain_char(i) {
    return substr(chains, i, 1)
}

function norm_resname(r) {
    r = trim(r)

    if (r == "WAT" || r == "SOL") return "HOH"
    if (r == "NA+" || r == "Na+" || r == "NA") return "NA"
    if (r == "CL-" || r == "Cl-" || r == "CL") return "CL"

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

function is_solvent_or_ion(r) {
    return (r == "HOH" || r == "NA" || r == "CL")
}

function guess_element(atom, res) {
    atom = trim(atom)
    gsub(/[0-9]/, "", atom)
    gsub(/\047/, "", atom)
    gsub(/\*/, "", atom)

    up = toupper(atom)

    if (res == "NA") return "NA"
    if (res == "CL") return "CL"

    if (substr(up,1,2) == "NA") return "NA"
    if (substr(up,1,2) == "CL") return "CL"

    return substr(up,1,1)
}

FNR == NR {
    if ($0 ~ /^ATOM  / || $0 ~ /^HETATM/) {
        line = pad80($0)

        res_raw = substr(line, 18, 3)
        resseq  = substr(line, 23, 4)
        icode   = substr(line, 27, 1)

        res = trim(res_raw)
        newres = norm_resname(res)
        res_uid = resseq icode

        if (res_uid != first_prev_res_uid) {
            if (res ~ /5$/ && is_nucleic(newres)) {
                nuc_chain_count++
            }

            if (is_ligand(newres)) {
                lig_count++
            }

            first_prev_res_uid = res_uid
        }
    }

    next
}

BEGIN {
    chain_index = 0
    current_chain = " "
    in_chain = 0
    prev_res_uid = ""
    prev_res_was_3 = 0
    lig_seen = 0
}

(/^ATOM  / || /^HETATM/) {
    line = pad80($0)

    atom    = substr(line, 13, 4)
    res_raw = substr(line, 18, 3)
    resseq  = substr(line, 23, 4)
    icode   = substr(line, 27, 1)

    res = trim(res_raw)
    newres = norm_resname(res)
    res_uid = resseq icode

    if (res_uid != prev_res_uid) {
        if (prev_res_was_3) {
            in_chain = 0
            current_chain = " "
            prev_res_was_3 = 0
        }

        if (res ~ /5$/ && is_nucleic(newres)) {
            chain_index++
            current_chain = chain_char(chain_index)
            in_chain = 1
        }

        if (res ~ /3$/ && is_nucleic(newres)) {
            prev_res_was_3 = 1
        }

        if (is_ligand(newres)) {
            lig_seen++
            ligand_chain = chain_char(nuc_chain_count + lig_seen)
        }

        prev_res_uid = res_uid
    }

    chain = " "

    if (in_chain && is_nucleic(newres)) {
        chain = current_chain
    } else if (is_ligand(newres)) {
        chain = ligand_chain
    } else if (is_solvent_or_ion(newres)) {
        chain = chain_char(nuc_chain_count + lig_count + 1)
    }

    elem = guess_element(atom, newres)

    printf "%s%3s%s%s%s%2s%s\n", \
        substr(line,1,17), \
        newres, \
        substr(line,21,1), \
        chain, \
        substr(line,23,54), \
        elem, \
        substr(line,79,2)

    next
}

{
    print
}
' md.pdb md.pdb > md.pdb.tmp

mv md.pdb.tmp md.pdb

cd ../..