#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ACPYPE_DIR="02_ACPYPE/leap.amb2gmx"
TOP_IN="${ACPYPE_DIR}/leap_GMX.top"
GRO_IN="${ACPYPE_DIR}/leap_GMX.gro"
OUT_DIR="03_GROMACS/00_TOPOLOGY"

if [[ ! -f "$TOP_IN" ]]; then
  echo "ERROR: $TOP_IN が見つかりません。" >&2
  echo "先に 03_ACPYPE.sh を実行してください。" >&2
  exit 1
fi

if [[ ! -f "$GRO_IN" ]]; then
  echo "ERROR: $GRO_IN が見つかりません。" >&2
  echo "先に 03_ACPYPE.sh を実行してください。" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

cp "$TOP_IN" "$OUT_DIR/leap_GMX.top.original"
cp "$GRO_IN" "$OUT_DIR/leap_GMX.gro"

TOP_ORIGINAL="$OUT_DIR/leap_GMX.top.original"
TOP_POSRES="$OUT_DIR/leap_GMX_posres.top"

FIRST_MOLNAME="$(awk '
BEGIN {
  mol = 0
  read_name = 0
}

/^[[:space:]]*\[[[:space:]]*moleculetype[[:space:]]*\]/ {
  mol++

  if (mol == 1) {
    read_name = 1
  }

  next
}

read_name && $0 !~ /^[[:space:]]*;/ && $0 !~ /^[[:space:]]*$/ {
  print $1
  exit
}
' "$TOP_ORIGINAL")"

MOLECULETYPE_COUNT="$(awk '
/^[[:space:]]*\[[[:space:]]*moleculetype[[:space:]]*\]/ {
  count++
}

END {
  print count + 0
}
' "$TOP_ORIGINAL")"

SOLUTE_ATOMS="$(awk '
BEGIN {
  mol = 0
  in_atoms = 0
}

/^[[:space:]]*\[[[:space:]]*moleculetype[[:space:]]*\]/ {
  mol++
}

/^[[:space:]]*\[[[:space:]]*atoms[[:space:]]*\]/ && mol == 1 {
  in_atoms = 1
  next
}

in_atoms && /^[[:space:]]*\[/ {
  exit
}

in_atoms && $1 ~ /^[0-9]+$/ {
  last_atom = $1
}

END {
  print last_atom
}
' "$TOP_ORIGINAL")"

FIRST_MOL_COUNT="$(awk -v target="$FIRST_MOLNAME" '
BEGIN {
  in_molecules = 0
}

/^[[:space:]]*\[[[:space:]]*molecules[[:space:]]*\]/ {
  in_molecules = 1
  next
}

in_molecules && /^[[:space:]]*\[/ {
  exit
}

in_molecules && $0 !~ /^[[:space:]]*;/ && NF >= 2 && $1 == target {
  print $2
  exit
}
' "$TOP_ORIGINAL")"

if [[ -z "$FIRST_MOLNAME" ]]; then
  echo "ERROR: 最初のmoleculetype名を取得できません。" >&2
  exit 1
fi

if [[ "$MOLECULETYPE_COUNT" -lt 2 ]]; then
  echo "ERROR: moleculetypeが2つ未満です。" >&2
  echo "想定したACPYPEトポロジーではありません。" >&2
  exit 1
fi

if [[ -z "$SOLUTE_ATOMS" || ! "$SOLUTE_ATOMS" =~ ^[0-9]+$ ]]; then
  echo "ERROR: 最初のmoleculetypeの原子数を取得できません。" >&2
  exit 1
fi

if [[ "$FIRST_MOL_COUNT" != "1" ]]; then
  echo "ERROR: [ molecules ] における $FIRST_MOLNAME の分子数が1ではありません。" >&2
  echo "取得値: ${FIRST_MOL_COUNT:-未取得}" >&2
  exit 1
fi

make_posres() {
  local force_constant="$1"
  local output_file="$2"

  awk -v fc="$force_constant" '
  BEGIN {
    mol = 0
    in_atoms = 0

    print "[ position_restraints ]"
    print "; atom  funct       fc_x       fc_y       fc_z"
  }

  /^[[:space:]]*\[[[:space:]]*moleculetype[[:space:]]*\]/ {
    mol++
  }

  /^[[:space:]]*\[[[:space:]]*atoms[[:space:]]*\]/ && mol == 1 {
    in_atoms = 1
    next
  }

  in_atoms && /^[[:space:]]*\[/ {
    exit
  }

  in_atoms && $1 ~ /^[0-9]+$/ && ($8 + 0.0) > 2.0 {
    printf "%8d %6d %12.3f %12.3f %12.3f\n", $1, 1, fc, fc, fc
  }
  ' "$TOP_ORIGINAL" > "$output_file"
}

make_posres 1000 "$OUT_DIR/posres_1000.itp"
make_posres 500 "$OUT_DIR/posres_500.itp"
make_posres 100 "$OUT_DIR/posres_100.itp"
make_posres 50 "$OUT_DIR/posres_50.itp"

awk '
BEGIN {
  mol = 0
}

/^[[:space:]]*\[[[:space:]]*moleculetype[[:space:]]*\]/ {
  mol++

  if (mol == 2) {
    print ""
    print "; Position restraints for the first moleculetype"

    print "#ifdef POSRES_1000"
    print "#include \"posres_1000.itp\""
    print "#endif"
    print ""

    print "#ifdef POSRES_500"
    print "#include \"posres_500.itp\""
    print "#endif"
    print ""

    print "#ifdef POSRES_100"
    print "#include \"posres_100.itp\""
    print "#endif"
    print ""

    print "#ifdef POSRES_50"
    print "#include \"posres_50.itp\""
    print "#endif"
    print ""
  }
}

{
  print
}
' "$TOP_ORIGINAL" > "$TOP_POSRES"

RESTRAINED_ATOMS="$(grep -cE '^[[:space:]]*[0-9]+' "$OUT_DIR/posres_1000.itp")"

LAST_RESTRAINED_ATOM="$(awk '
$1 ~ /^[0-9]+$/ {
  last = $1
}

END {
  print last
}
' "$OUT_DIR/posres_1000.itp")"

for force_constant in 1000 500 100 50; do
  current_count="$(grep -cE '^[[:space:]]*[0-9]+' "$OUT_DIR/posres_${force_constant}.itp")"

  if [[ "$current_count" != "$RESTRAINED_ATOMS" ]]; then
    echo "ERROR: posres_${force_constant}.itp の拘束原子数が一致しません。" >&2
    exit 1
  fi
done

if [[ -z "$LAST_RESTRAINED_ATOM" || "$LAST_RESTRAINED_ATOM" -gt "$SOLUTE_ATOMS" ]]; then
  echo "ERROR: 拘束原子番号が溶質原子数を超えています。" >&2
  exit 1
fi

for force_constant in 1000 500 100 50; do
  if ! grep -q "POSRES_${force_constant}" "$TOP_POSRES"; then
    echo "ERROR: POSRES_${force_constant} の挿入に失敗しました。" >&2
    exit 1
  fi
done

echo
echo "========================================"
echo "Position-restraint topology prepared"
echo "========================================"
echo "First moleculetype : $FIRST_MOLNAME"
echo "Number of copies   : $FIRST_MOL_COUNT"
echo "Total solute atoms : $SOLUTE_ATOMS"
echo "Restrained atoms   : $RESTRAINED_ATOMS"
echo "Last restrained ID : $LAST_RESTRAINED_ATOM"
echo "Force constants    : 1000, 500, 100, 50"
echo "Output topology    : $TOP_POSRES"
echo

echo "Insertion check:"
grep -n -B2 -A22 'POSRES_1000' "$TOP_POSRES"