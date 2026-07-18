#!/usr/bin/env bash
set -euo pipefail

# model.pdb の残基名から DNA/RNA と LIG の有無を判定する。
[[ -f model.pdb ]] || {
  echo "Error: model.pdb が見つかりません。" >&2
  exit 1
}

read -r has_dna has_rna has_lig < <(
  awk '
    /^(ATOM  |HETATM)/ {
      resname = substr($0, 18, 3)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", resname)

      if (resname == "DA" || resname == "DC" ||
          resname == "DG" || resname == "DT") {
        has_dna = 1
      }
      if (resname == "A" || resname == "C" ||
          resname == "G" || resname == "T") {
        has_rna = 1
      }
      if (resname == "LIG") {
        has_lig = 1
      }
    }
    END {
      print has_dna + 0, has_rna + 0, has_lig + 0
    }
  ' model.pdb
)

if (( has_dna && has_rna )); then
  echo "Error: model.pdb に DNA 残基と RNA 残基の両方が含まれています。" >&2
  exit 1
elif (( has_dna )); then
  molecule_type="DNA"
elif (( has_rna )); then
  molecule_type="RNA"
else
  echo "Error: model.pdb から DNA または RNA を判定できませんでした。" >&2
  exit 1
fi

if (( has_lig )); then
  molecule_type+="-LIG"
fi

pre_leap_source="pre-leap_${molecule_type}.in"
leap_source="leap_${molecule_type}.in"

for source_file in "$pre_leap_source" "$leap_source"; do
  [[ -f "$source_file" ]] || {
    echo "Error: $source_file が見つかりません。" >&2
    exit 1
  }
done

cp -- "$pre_leap_source" pre-leap.in
cp -- "$leap_source" leap.in
printf 'Detected type: %s\n' "$molecule_type"
printf 'Copied: %s -> pre-leap.in\n' "$pre_leap_source"
printf 'Copied: %s -> leap.in\n' "$leap_source"

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
