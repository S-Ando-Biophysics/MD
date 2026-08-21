#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

INPUT_DIR="${ROOT_DIR}/04_ANALYSIS/00_ALL-FRAMES"

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found:"
    echo "  $INPUT_DIR"
    exit 1
fi

NUMBER_OF_PDB_FILES=$(
    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.pdb" \
    | wc -l
)

if [ "$NUMBER_OF_PDB_FILES" -eq 0 ]; then
    echo "ERROR: No PDB files found in:"
    echo "  $INPUT_DIR"
    exit 1
fi

while true; do

    read -rp "Number of chains: " NUMBER_OF_CHAINS

    if [[ "$NUMBER_OF_CHAINS" =~ ^[1-9][0-9]*$ ]]; then
        break
    fi

    echo "ERROR: Please enter a positive integer."

done

declare -A USED_CHAIN_IDS=()
declare -A RESIDUE_TO_CHAIN=()

CHAIN_IDS=()
CHAIN_RESIDUE_SPECS=()

for ((CHAIN_INDEX=1; CHAIN_INDEX<=NUMBER_OF_CHAINS; CHAIN_INDEX++)); do

    echo
    echo "----------------------------------------"
    echo " Chain ${CHAIN_INDEX}/${NUMBER_OF_CHAINS}"
    echo "----------------------------------------"

    while true; do

        read -rp "Chain ID (one character): " CHAIN_ID

        if [[ ! "$CHAIN_ID" =~ ^[A-Za-z0-9]$ ]]; then
            echo "ERROR: Chain ID must be one alphanumeric character."
            continue
        fi

        if [[ -n "${USED_CHAIN_IDS[$CHAIN_ID]+x}" ]]; then
            echo "ERROR: Chain ID '${CHAIN_ID}' is already used."
            continue
        fi

        break

    done

    while true; do

        read -rp "Residues (e.g. 1-12,15,18-20): " RESIDUE_SPEC

        RESIDUE_SPEC="${RESIDUE_SPEC// /}"

        if [ -z "$RESIDUE_SPEC" ]; then
            echo "ERROR: Residue selection cannot be empty."
            continue
        fi

        VALID_SELECTION=true
        CANDIDATE_RESIDS=()

        unset CANDIDATE_SEEN
        declare -A CANDIDATE_SEEN=()

        IFS=',' read -ra TOKENS <<< "$RESIDUE_SPEC"

        for TOKEN in "${TOKENS[@]}"; do

            if [[ "$TOKEN" =~ ^([0-9]+)-([0-9]+)$ ]]; then

                START_RESID=$((10#${BASH_REMATCH[1]}))
                END_RESID=$((10#${BASH_REMATCH[2]}))

                if [ "$START_RESID" -gt "$END_RESID" ]; then
                    echo "ERROR: Invalid residue range:"
                    echo "  $TOKEN"
                    VALID_SELECTION=false
                    break
                fi

                for ((RESID=START_RESID; RESID<=END_RESID; RESID++)); do

                    if [[ -n "${RESIDUE_TO_CHAIN[$RESID]+x}" ]]; then
                        echo "ERROR: Residue ${RESID} is already assigned to chain ${RESIDUE_TO_CHAIN[$RESID]}."
                        VALID_SELECTION=false
                        break 2
                    fi

                    if [[ -z "${CANDIDATE_SEEN[$RESID]+x}" ]]; then
                        CANDIDATE_RESIDS+=("$RESID")
                        CANDIDATE_SEEN[$RESID]=1
                    fi

                done

            elif [[ "$TOKEN" =~ ^[0-9]+$ ]]; then

                RESID=$((10#$TOKEN))

                if [[ -n "${RESIDUE_TO_CHAIN[$RESID]+x}" ]]; then
                    echo "ERROR: Residue ${RESID} is already assigned to chain ${RESIDUE_TO_CHAIN[$RESID]}."
                    VALID_SELECTION=false
                    break
                fi

                if [[ -z "${CANDIDATE_SEEN[$RESID]+x}" ]]; then
                    CANDIDATE_RESIDS+=("$RESID")
                    CANDIDATE_SEEN[$RESID]=1
                fi

            else

                echo "ERROR: Invalid residue selection:"
                echo "  $TOKEN"
                VALID_SELECTION=false
                break

            fi

        done

        if [ "$VALID_SELECTION" = true ]; then
            break
        fi

    done

    USED_CHAIN_IDS[$CHAIN_ID]=1

    CHAIN_IDS+=("$CHAIN_ID")
    CHAIN_RESIDUE_SPECS+=("$RESIDUE_SPEC")

    for RESID in "${CANDIDATE_RESIDS[@]}"; do
        RESIDUE_TO_CHAIN[$RESID]="$CHAIN_ID"
    done

done

MAP_FILE="$(mktemp)"

trap 'rm -f "$MAP_FILE"' EXIT

{
    for RESID in "${!RESIDUE_TO_CHAIN[@]}"; do
        printf "%s %s\n" \
            "$RESID" \
            "${RESIDUE_TO_CHAIN[$RESID]}"
    done
} | sort -n > "$MAP_FILE"

echo
echo "========================================"
echo " PDB format cleanup "
echo "========================================"
echo
echo "Chain assignments:"

for ((CHAIN_INDEX=0; CHAIN_INDEX<NUMBER_OF_CHAINS; CHAIN_INDEX++)); do
    echo "  Chain ${CHAIN_IDS[$CHAIN_INDEX]} : ${CHAIN_RESIDUE_SPECS[$CHAIN_INDEX]}"
done

echo
echo "Number of PDB files: $NUMBER_OF_PDB_FILES"
echo

CURRENT_FILE=0

while IFS= read -r -d '' PDB_FILE; do

    CURRENT_FILE=$((CURRENT_FILE + 1))

    PDB_NAME="$(basename "$PDB_FILE")"
    TEMP_FILE="${PDB_FILE}.tmp"

    echo "${CURRENT_FILE}/${NUMBER_OF_PDB_FILES} : ${PDB_NAME}"

    if ! awk -v map_file="$MAP_FILE" '

        function trim(text) {
            gsub(/^[[:space:]]+/, "", text)
            gsub(/[[:space:]]+$/, "", text)
            return text
        }

        function pad80(line) {
            while (length(line) < 80) {
                line = line " "
            }
            return line
        }

        function infer_element(atom_name, residue_name, atom, residue, first, i, character) {

            atom = toupper(trim(atom_name))
            residue = toupper(trim(residue_name))

            if (atom == "") {
                return ""
            }

            if (atom == residue && atom ~ /^(NA|MG|CL|CA|ZN|FE|MN|CU|CO|NI|BR|K|LI|AL|SI|SE|CD|HG|PB)$/) {
                return atom
            }

            first = ""

            for (i = 1; i <= length(atom); i++) {

                character = substr(atom, i, 1)

                if (character ~ /[A-Z]/) {
                    first = character
                    break
                }

            }

            if (first ~ /^(H|C|N|O|P|S|F|I|B)$/) {
                return first
            }

            return ""
        }

        BEGIN {

            while ((getline < map_file) > 0) {
                residue_chain[$1] = $2
            }

            close(map_file)

        }

        {

            line = $0
            record = substr(line, 1, 6)

            if (record == "ATOM  " || record == "HETATM") {

                line = pad80(line)

                resid_text = trim(substr(line, 23, 4))

                if (resid_text ~ /^-?[0-9]+$/) {

                    resid = resid_text + 0

                    if (resid in residue_chain) {
                        line = substr(line, 1, 21) residue_chain[resid] substr(line, 23)
                    }

                }

                element = trim(substr(line, 77, 2))

                if (element == "") {

                    element = infer_element(substr(line, 13, 4), substr(line, 18, 3))

                    if (element != "") {
                        line = substr(line, 1, 76) sprintf("%2s", element) substr(line, 79)
                    }

                }

            }

            else if (substr(line, 1, 3) == "TER") {

                line = pad80(line)

                resid_text = trim(substr(line, 23, 4))

                if (resid_text ~ /^-?[0-9]+$/) {

                    resid = resid_text + 0

                    if (resid in residue_chain) {
                        line = substr(line, 1, 21) residue_chain[resid] substr(line, 23)
                    }

                }

            }

            print line

        }

    ' "$PDB_FILE" > "$TEMP_FILE"; then

        rm -f "$TEMP_FILE"

        echo "ERROR: Failed to process:"
        echo "  $PDB_FILE"

        exit 1

    fi

    mv "$TEMP_FILE" "$PDB_FILE"

done < <(
    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.pdb" \
        -print0 \
    | sort -z
)

echo
echo "========================================"
echo " Finished "
echo "========================================"
echo
echo "Processed PDB files: $NUMBER_OF_PDB_FILES"
echo
echo "Updated directory:"
echo "  $INPUT_DIR"
echo