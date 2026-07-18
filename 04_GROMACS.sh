#!/usr/bin/env bash
set -euo pipefail

TOP_DIR="03_GROMACS/00_TOPOLOGY"
TOP_NAME="leap_GMX_posres.top"

EM_POSRES_DIR="03_GROMACS/01_EM/01_POSRES1000"
EM_FREE_DIR="03_GROMACS/01_EM/02_FREE"

NVT_POSRES_DIR="03_GROMACS/02_NVT/01_POSRES1000"

NPT_1000_DIR="03_GROMACS/03_NPT/01_POSRES1000"
NPT_500_DIR="03_GROMACS/03_NPT/02_POSRES500"
NPT_100_DIR="03_GROMACS/03_NPT/03_POSRES100"
NPT_50_DIR="03_GROMACS/03_NPT/04_POSRES50"
NPT_FREE_DIR="03_GROMACS/03_NPT/05_FREE"

MD_DIR="03_GROMACS/04_MD"

mkdir -p "$EM_POSRES_DIR"
mkdir -p "$EM_FREE_DIR"
mkdir -p "$NVT_POSRES_DIR"
mkdir -p "$NPT_1000_DIR"
mkdir -p "$NPT_500_DIR"
mkdir -p "$NPT_100_DIR"
mkdir -p "$NPT_50_DIR"
mkdir -p "$NPT_FREE_DIR"
mkdir -p "$MD_DIR"

bash prepare_posres.sh

copy_topology() {
  local destination="$1"

  cp "$TOP_DIR/$TOP_NAME" "$destination/"
  cp "$TOP_DIR/posres_1000.itp" "$destination/"
  cp "$TOP_DIR/posres_500.itp" "$destination/"
  cp "$TOP_DIR/posres_100.itp" "$destination/"
  cp "$TOP_DIR/posres_50.itp" "$destination/"
}

echo "=== EM: POSRES 1000 ==="

copy_topology "$EM_POSRES_DIR"
cp "$TOP_DIR/leap_GMX.gro" "$EM_POSRES_DIR/"
cp em_posres1000.mdp "$EM_POSRES_DIR/"

(
  cd "$EM_POSRES_DIR"
  gmx grompp -f em_posres1000.mdp -c leap_GMX.gro -r leap_GMX.gro -p "$TOP_NAME" -o em_posres1000.tpr
  gmx mdrun -v -deffnm em_posres1000
)

echo "=== EM: FREE ==="

copy_topology "$EM_FREE_DIR"
cp "$EM_POSRES_DIR/em_posres1000.gro" "$EM_FREE_DIR/"
cp em_free.mdp "$EM_FREE_DIR/"

(
  cd "$EM_FREE_DIR"
  gmx grompp -f em_free.mdp -c em_posres1000.gro -p "$TOP_NAME" -o em_free.tpr
  gmx mdrun -v -deffnm em_free
)

REFERENCE_STRUCTURE="$EM_FREE_DIR/em_free.gro"

echo "=== NVT: POSRES 1000 ==="

copy_topology "$NVT_POSRES_DIR"
cp "$REFERENCE_STRUCTURE" "$NVT_POSRES_DIR/start.gro"
cp "$REFERENCE_STRUCTURE" "$NVT_POSRES_DIR/reference.gro"
cp nvt_posres1000.mdp "$NVT_POSRES_DIR/"

(
  cd "$NVT_POSRES_DIR"
  gmx grompp -f nvt_posres1000.mdp -c start.gro -r reference.gro -p "$TOP_NAME" -o nvt_posres1000.tpr
  gmx mdrun -v -deffnm nvt_posres1000
)

run_restrained_npt() {
  local stage_directory="$1"
  local mdp_file="$2"
  local input_gro="$3"
  local input_cpt="$4"
  local output_name="$5"

  copy_topology "$stage_directory"

  cp "$input_gro" "$stage_directory/start.gro"
  cp "$input_cpt" "$stage_directory/start.cpt"
  cp "$REFERENCE_STRUCTURE" "$stage_directory/reference.gro"
  cp "$mdp_file" "$stage_directory/"

  (
    cd "$stage_directory"
    gmx grompp -f "$(basename "$mdp_file")" -c start.gro -r reference.gro -t start.cpt -p "$TOP_NAME" -o "${output_name}.tpr"
    gmx mdrun -v -deffnm "$output_name"
  )
}

echo "=== NPT: POSRES 1000 ==="

run_restrained_npt "$NPT_1000_DIR" "npt_posres1000.mdp" "$NVT_POSRES_DIR/nvt_posres1000.gro" "$NVT_POSRES_DIR/nvt_posres1000.cpt" "npt_posres1000"

echo "=== NPT: POSRES 500 ==="

run_restrained_npt "$NPT_500_DIR" "npt_posres500.mdp" "$NPT_1000_DIR/npt_posres1000.gro" "$NPT_1000_DIR/npt_posres1000.cpt" "npt_posres500"

echo "=== NPT: POSRES 100 ==="

run_restrained_npt "$NPT_100_DIR" "npt_posres100.mdp" "$NPT_500_DIR/npt_posres500.gro" "$NPT_500_DIR/npt_posres500.cpt" "npt_posres100"

echo "=== NPT: POSRES 50 ==="

run_restrained_npt "$NPT_50_DIR" "npt_posres50.mdp" "$NPT_100_DIR/npt_posres100.gro" "$NPT_100_DIR/npt_posres100.cpt" "npt_posres50"

echo "=== NPT: FREE ==="

copy_topology "$NPT_FREE_DIR"
cp "$NPT_50_DIR/npt_posres50.gro" "$NPT_FREE_DIR/start.gro"
cp "$NPT_50_DIR/npt_posres50.cpt" "$NPT_FREE_DIR/start.cpt"
cp npt_free.mdp "$NPT_FREE_DIR/"

(
  cd "$NPT_FREE_DIR"
  gmx grompp -f npt_free.mdp -c start.gro -t start.cpt -p "$TOP_NAME" -o npt_free.tpr
  gmx mdrun -v -deffnm npt_free
)

echo "=== Production MD ==="

copy_topology "$MD_DIR"
cp "$NPT_FREE_DIR/npt_free.gro" "$MD_DIR/start.gro"
cp "$NPT_FREE_DIR/npt_free.cpt" "$MD_DIR/start.cpt"
cp md.mdp "$MD_DIR/"

(
  cd "$MD_DIR"
  gmx grompp -f md.mdp -c start.gro -t start.cpt -p "$TOP_NAME" -o md.tpr
  gmx mdrun -v -deffnm md -nb auto -nbfe auto -pme auto -pmefft auto -bonded auto -update auto -cpi md.cpt
)

echo "=== Finished ==="