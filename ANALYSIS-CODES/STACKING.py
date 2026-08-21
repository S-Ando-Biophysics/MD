#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import re

import MDAnalysis as mda
import numpy as np
import pandas as pd



# ============================================================
# User setting 1/2: Target
# ============================================================

# "RESIDUE_SEQUENCES" defines ordered residue sequences.
#
#   Stacking is evaluated between consecutive residues in each sequence.
#   For example, if a sequence is [1, 2, 3, 4], the following base steps are analyzed: 1-2, 2-3, 3-4.
#
#   Residues do not need to be consecutive residue numbers. 
#   For example, if a sequence is [1, 3, 5], the following base steps are analyzed: 1-3, 3-5.


# "BASE_STEPS" defines individual base steps explicitly.
#
#    Each entry contains exactly two residue numbers.
#    BASE_STEPS can be analyzed in addition to RESIDUE_SEQUENCES.
#

# To use only "RESIDUE_SEQUENCES", set "BASE_STEPS" to an empty dictionary.
# To use only "BASE_STEPS", set "RESIDUE_SEQUENCES" to an empty dictionary.

# Note: Python "range(start, stop)" does not include stop.

RESIDUE_SEQUENCES = {
    "strand_1": list(range(1, 13)),
    "strand_2": list(range(13, 25)),
}

BASE_STEPS = {

}
 
# Examples:

# 1) Consecutive residue sequences
#
# RESIDUE_SEQUENCES = {
#     "strand_1": list(range(1, 13)),
#     "strand_2": list(range(13, 25)),
# }
#
# BASE_STEPS = {
# }
#
#  This analyzes:
#     strand_1 : 1-2, 2-3, 3-4, ..., 11-12
#     strand_2 : 13-14, 14-15, 15-16, ..., 23-24
#

# 2) Manually specified residue sequence
#
# RESIDUE_SEQUENCES = {
#     "selected_region": [3, 4, 5, 6, 7],
# }
#
# BASE_STEPS = {
# }
#
# This analyzes:
#     3-4, 4-5, 5-6, 6-7
#

# 3) Explicit base steps
#
# RESIDUE_SEQUENCES = {
# }
#
# BASE_STEPS = {
#     "step_1": (1, 5),
#     "step_2": (3, 8),
#     "step_3": (12, 13),
# }
#
# This analyzes:
#     1-5, 3-8, 12-13
#

# 4) Combination
#
# RESIDUE_SEQUENCES = {
#     "strand_1": list(range(1, 13)),
# }
#
# BASE_STEPS = {
#     "additional_step": (3, 15),
# }
#
# This analyzes:
#     1-2, 2-3, 3-4, ..., 11-12
#     and 3-15
#

# ============================================================



# ============================================================
# User setting 2/2: Trajectory
# ============================================================

# Python uses zero-based indexing. (0 = Frame 1, 1 = Frame 2, ...)

# Frame 1 is excluded here because it is the initial structure.
# If you do not want to exclude it, please set to 0.
START_FRAME_INDEX = 1

# None = analyze through the final frame
STOP_FRAME_INDEX = None

# 1 = all frames
# 10 = every 10th frame
STRIDE = 1

# ============================================================





SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent

INPUT_DIR = ROOT_DIR / "03_GROMACS" / "06_OUTPUT"
ANALYSIS_ROOT = ROOT_DIR / "04_ANALYSIS"

def latest_file(directory: Path, pattern: str) -> Path:

    files = [path for path in directory.glob(pattern) if path.is_file()]

    if not files:
        raise FileNotFoundError(
            f"Input file not found.\n"
            f"Directory: {directory}\n"
            f"Pattern: {pattern}"
        )

    return max(files, key=lambda path: path.stat().st_mtime)

if not INPUT_DIR.is_dir():
    raise FileNotFoundError(
        f"Input directory not found.\n"
        f"{INPUT_DIR}"
    )

TOPOLOGY_FILE = latest_file(INPUT_DIR, "*.pdb")
TRAJECTORY_FILE = latest_file(INPUT_DIR, "*.xtc")

ANALYSIS_ROOT.mkdir(parents=True, exist_ok=True)

analysis_numbers = []

for path in ANALYSIS_ROOT.iterdir():
    if not path.is_dir():
        continue

    match = re.match(r"^(\d+)_", path.name)

    if match:
        analysis_numbers.append(int(match.group(1)))

next_analysis_number = (
    max(analysis_numbers) + 1
    if analysis_numbers
    else 1
)

OUTPUT_DIR = (
    ANALYSIS_ROOT
    / f"{next_analysis_number:02d}_STACKING"
)

SUMMARY_CSV = OUTPUT_DIR / "stacking_summary.csv"
TIMESERIES_CSV = OUTPUT_DIR / "stacking_timeseries.csv"
DIAGNOSTICS_CSV = OUTPUT_DIR / "stacking_diagnostics.csv"

# ============================================================
# Stacking criteria
# ============================================================
#
# A frame is classified as stacked only when all three criteria are satisfied.
#
# 1. Minimum heavy-atom distance between the two bases < 4 A
# 2. Center-of-mass distance between the two bases < 5 A
# 3. Angle between the two base-plane normal vectors:
#       0-45 degrees or 135-180 degrees
#
# (Brown RF et al. J. Chem. Theory Comput. 2015, 11(5), 2315-2328. https://doi.org/10.1021/ct501170h)

MIN_HEAVY_ATOM_DISTANCE_CUTOFF_A = 4.0
COM_DISTANCE_CUTOFF_A = 5.0
NORMAL_ANGLE_CUTOFF_DEG = 45.0

# ============================================================
# Base atom definitions
# ============================================================

PURINE_RING_ATOMS = (
    "N9", "C8", "N7", "C5", "C6",
    "N1", "C2", "N3", "C4",
)

PYRIMIDINE_RING_ATOMS = (
    "N1", "C2", "N3", "C4", "C5", "C6",
)

BASE_HEAVY_ATOM_NAMES = {
    # Common ring atoms
    "N1", "C2", "N3", "C4", "C5", "C6",

    # Purine-specific ring atoms
    "N7", "C8", "N9",

    # Exocyclic atoms
    "N2", "N4", "N6",
    "O2", "O4", "O6",

    # Thymine 5-methyl carbon: C7 or C5M depending on the force field
    "C7", "C5M",
}

# ============================================================
# Atomic masses
# ============================================================

ELEMENT_MASS = {
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
}

# ============================================================
# Analysis
# ============================================================

def check_input_file(file_path: Path) -> None:
    if not file_path.exists():
        raise FileNotFoundError(
            "Input file not found.\n"
            f"{file_path}"
        )

check_input_file(TOPOLOGY_FILE)
check_input_file(TRAJECTORY_FILE)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("Loading trajectory")
print(f"Topology   : {TOPOLOGY_FILE}")
print(f"Trajectory : {TRAJECTORY_FILE}")

universe = mda.Universe(
    str(TOPOLOGY_FILE),
    str(TRAJECTORY_FILE),
)

number_of_total_frames = len(universe.trajectory)

stop_frame = (
    number_of_total_frames
    if STOP_FRAME_INDEX is None
    else min(STOP_FRAME_INDEX, number_of_total_frames)
)

analysis_frame_indices = np.arange(
    START_FRAME_INDEX,
    stop_frame,
    STRIDE,
    dtype=int,
)

number_of_analysis_frames = len(analysis_frame_indices)

if number_of_analysis_frames == 0:
    raise ValueError("No frames are available for analysis.")

print(f"Total frames : {number_of_total_frames}")
print(f"Used frames  : {number_of_analysis_frames}")

def atom_element_from_name(atom_name: str) -> str:
    """
    Infer the element symbol of a nucleobase heavy atom from its atom name.
    Supported elements are C, N, and O.
    """
    name = atom_name.strip().upper()

    for char in name:
        if char in ELEMENT_MASS:
            return char

    raise ValueError(
        f"Could not determine element from atom name: {atom_name}"
    )

def base_heavy_atoms(residue: mda.core.groups.Residue) -> mda.core.groups.AtomGroup:
    """
    Return nucleobase heavy atoms from a residue.
    """
    selected = residue.atoms[
        np.array(
            [
                atom.name.strip().upper() in BASE_HEAVY_ATOM_NAMES
                for atom in residue.atoms
            ],
            dtype=bool,
        )
    ]

    if len(selected) < 5:
        available_names = ", ".join(
            atom.name for atom in residue.atoms
        )
        raise ValueError(
            "Could not select enough base heavy atoms.\n"
            f"resid={residue.resid}, resname={residue.resname}\n"
            f"selected atoms={len(selected)}\n"
            f"available atoms={available_names}"
        )

    return selected

def detect_ring_atom_names(
    base_atoms: mda.core.groups.AtomGroup,
) -> tuple[str, ...]:
    """
    Determine whether the base is a purine or pyrimidine from atom names,
    and return the ring atom names used for plane fitting.

    If N9, N7, and C8 are present, classify the base as a purine.
    """
    names = {
        atom.name.strip().upper()
        for atom in base_atoms
    }

    purine_signature = {"N9", "C8", "N7"}

    if purine_signature.issubset(names):
        ring_names = PURINE_RING_ATOMS
    else:
        ring_names = PYRIMIDINE_RING_ATOMS

    missing = [
        name
        for name in ring_names
        if name not in names
    ]

    if missing:
        raise ValueError(
            "Could not identify all required base-ring atoms.\n"
            f"available={sorted(names)}\n"
            f"missing={missing}"
        )

    return ring_names

def ring_atoms(
    base_atoms: mda.core.groups.AtomGroup,
    ring_names: tuple[str, ...],
) -> mda.core.groups.AtomGroup:
    """
    Return only ring atoms from the base heavy-atom group.
    """
    atom_by_name = {
        atom.name.strip().upper(): atom
        for atom in base_atoms
    }

    atoms = [
        atom_by_name[name]
        for name in ring_names
    ]

    return mda.core.groups.AtomGroup(atoms)

def center_of_mass_manual(
    atoms: mda.core.groups.AtomGroup,
) -> np.ndarray:
    """
    Calculate the center of mass using standard atomic masses for C/N/O.
    """
    positions = atoms.positions.astype(np.float64)

    masses = np.array(
        [
            ELEMENT_MASS[
                atom_element_from_name(atom.name)
            ]
            for atom in atoms
        ],
        dtype=np.float64,
    )

    return np.average(
        positions,
        axis=0,
        weights=masses,
    )

def minimum_heavy_atom_distance(
    atoms_1: mda.core.groups.AtomGroup,
    atoms_2: mda.core.groups.AtomGroup,
) -> float:
    """
    Minimum heavy-atom distance between two bases (A).
    """
    xyz_1 = atoms_1.positions.astype(np.float64)
    xyz_2 = atoms_2.positions.astype(np.float64)

    difference = (
        xyz_1[:, np.newaxis, :]
        - xyz_2[np.newaxis, :, :]
    )

    squared_distance = np.sum(
        difference * difference,
        axis=2,
    )

    return float(
        np.sqrt(
            np.min(squared_distance)
        )
    )

def center_of_mass_distance(
    atoms_1: mda.core.groups.AtomGroup,
    atoms_2: mda.core.groups.AtomGroup,
) -> float:
    """
    Center-of-mass distance between two bases (A).
    """
    com_1 = center_of_mass_manual(atoms_1)
    com_2 = center_of_mass_manual(atoms_2)

    return float(
        np.linalg.norm(com_1 - com_2)
    )

def base_plane_normal(
    atoms: mda.core.groups.AtomGroup,
) -> np.ndarray:
    """
    Fit a least-squares plane to the ring atoms
    and return the unit normal vector.

    The minimum-variance direction from SVD is used as the normal.
    """
    xyz = atoms.positions.astype(np.float64)

    center = np.mean(
        xyz,
        axis=0,
    )

    centered = xyz - center

    _, _, vh = np.linalg.svd(
        centered,
        full_matrices=False,
    )

    normal = vh[-1].astype(np.float64)

    norm = np.linalg.norm(normal)

    if norm == 0:
        raise ValueError(
            "Could not calculate the base-plane normal."
        )

    return normal / norm

def normal_angle_deg(
    ring_1: mda.core.groups.AtomGroup,
    ring_2: mda.core.groups.AtomGroup,
) -> float:
    """
    Return the angle between two base-plane normal vectors in the range 0-180 degrees.
    """
    normal_1 = base_plane_normal(ring_1)
    normal_2 = base_plane_normal(ring_2)

    cosine = float(
        np.dot(normal_1, normal_2)
    )

    cosine = float(
        np.clip(cosine, -1.0, 1.0)
    )

    return float(
        np.degrees(
            np.arccos(cosine)
        )
    )

def is_brown_stacked(
    minimum_distance_a: float,
    com_distance_a: float,
    angle_deg: float,
) -> bool:
    """
    Apply the three Brown et al. stacking criteria.
    """
    criterion_1 = (
        minimum_distance_a
        < MIN_HEAVY_ATOM_DISTANCE_CUTOFF_A
    )

    criterion_2 = (
        com_distance_a
        < COM_DISTANCE_CUTOFF_A
    )

    criterion_3 = (
        angle_deg <= NORMAL_ANGLE_CUTOFF_DEG
        or angle_deg >= (
            180.0 - NORMAL_ANGLE_CUTOFF_DEG
        )
    )

    return bool(
        criterion_1
        and criterion_2
        and criterion_3
    )

step_records = []

def select_unique_residue(resid: int) -> mda.core.groups.Residue:
    selected_residues = universe.select_atoms(
        f"resid {resid}"
    ).residues

    if len(selected_residues) != 1:
        raise ValueError(
            f"Could not uniquely select resid {resid}.\n"
            f"Selected residues = {len(selected_residues)}"
        )

    return selected_residues[0]

def add_step_record(
    group_name: str,
    step_index: int,
    resid_1: int,
    resid_2: int,
) -> None:
    residue_1 = select_unique_residue(
        resid_1
    )
    residue_2 = select_unique_residue(
        resid_2
    )

    base_1 = base_heavy_atoms(
        residue_1
    )
    base_2 = base_heavy_atoms(
        residue_2
    )

    ring_names_1 = detect_ring_atom_names(
        base_1
    )
    ring_names_2 = detect_ring_atom_names(
        base_2
    )

    ring_1 = ring_atoms(
        base_1,
        ring_names_1,
    )

    ring_2 = ring_atoms(
        base_2,
        ring_names_2,
    )

    step_records.append(
        {
            "sequence": group_name,
            "step_index_in_sequence": step_index,
            "residue_1": residue_1,
            "residue_2": residue_2,
            "base_1": base_1,
            "base_2": base_2,
            "ring_1": ring_1,
            "ring_2": ring_2,
        }
    )

for sequence_name, residue_ids in RESIDUE_SEQUENCES.items():

    for index in range(
        len(residue_ids) - 1
    ):
        add_step_record(
            sequence_name,
            index + 1,
            residue_ids[index],
            residue_ids[index + 1],
        )

for step_index, (step_name, residue_pair) in enumerate(
    BASE_STEPS.items(),
    start=1,
):
    if len(residue_pair) != 2:
        raise ValueError(
            f"BASE_STEPS['{step_name}'] must contain exactly two residue IDs."
        )

    add_step_record(
        step_name,
        step_index,
        residue_pair[0],
        residue_pair[1],
    )

number_of_steps = len(step_records)

print()
print("Base steps")

for step_index, record in enumerate(
    step_records,
    start=1,
):
    residue_1 = record["residue_1"]
    residue_2 = record["residue_2"]

    print(
        f"{step_index:2d}  "
        f"{record['sequence']:10s}  "
        f"{residue_1.resid}:{residue_1.resname}"
        f" - "
        f"{residue_2.resid}:{residue_2.resname}"
    )

stacked_matrix = np.zeros(
    (
        number_of_analysis_frames,
        number_of_steps,
    ),
    dtype=np.uint8,
)

minimum_distance_matrix = np.empty(
    (
        number_of_analysis_frames,
        number_of_steps,
    ),
    dtype=np.float32,
)

com_distance_matrix = np.empty(
    (
        number_of_analysis_frames,
        number_of_steps,
    ),
    dtype=np.float32,
)

normal_angle_matrix = np.empty(
    (
        number_of_analysis_frames,
        number_of_steps,
    ),
    dtype=np.float32,
)

times_ps = np.empty(
    number_of_analysis_frames,
    dtype=np.float64,
)

print()
print("Running Brown stacking analysis")

for analysis_row, trajectory_index in enumerate(
    analysis_frame_indices
):
    timestep = universe.trajectory[
        trajectory_index
    ]

    times_ps[analysis_row] = float(
        timestep.time
    )

    for step_column, record in enumerate(
        step_records
    ):
        minimum_distance_a = (
            minimum_heavy_atom_distance(
                record["base_1"],
                record["base_2"],
            )
        )

        com_distance_a = (
            center_of_mass_distance(
                record["base_1"],
                record["base_2"],
            )
        )

        angle_deg = normal_angle_deg(
            record["ring_1"],
            record["ring_2"],
        )

        stacked = is_brown_stacked(
            minimum_distance_a,
            com_distance_a,
            angle_deg,
        )

        minimum_distance_matrix[
            analysis_row,
            step_column,
        ] = minimum_distance_a

        com_distance_matrix[
            analysis_row,
            step_column,
        ] = com_distance_a

        normal_angle_matrix[
            analysis_row,
            step_column,
        ] = angle_deg

        stacked_matrix[
            analysis_row,
            step_column,
        ] = int(stacked)

    if (
        analysis_row == 0
        or (analysis_row + 1) % 1000 == 0
        or analysis_row + 1
        == number_of_analysis_frames
    ):
        print(
            f"\r"
            f"{analysis_row + 1}"
            f"/{number_of_analysis_frames}",
            end="",
            flush=True,
        )

print()

summary_rows = []

for step_column, record in enumerate(
    step_records
):
    residue_1 = record["residue_1"]
    residue_2 = record["residue_2"]

    stacked_frames = int(
        np.sum(
            stacked_matrix[:, step_column]
        )
    )

    stacking_percentage = (
        stacked_frames
        / number_of_analysis_frames
        * 100.0
    )

    summary_rows.append(
        {
            "Sequence": record["sequence"],
            "Step_in_sequence": (
                record["step_index_in_sequence"]
            ),
            "Step_label": (
                f"{residue_1.resid}-{residue_2.resid}"
            ),
            "Resid_1": int(residue_1.resid),
            "Resname_1": str(residue_1.resname),
            "Resid_2": int(residue_2.resid),
            "Resname_2": str(residue_2.resname),
            "Stacked_frames": stacked_frames,
            "Total_frames": (
                number_of_analysis_frames
            ),
            "Stacking_percentage": (
                stacking_percentage
            ),
            "Mean_minimum_distance_A": float(
                np.mean(
                    minimum_distance_matrix[
                        :,
                        step_column,
                    ]
                )
            ),
            "Mean_COM_distance_A": float(
                np.mean(
                    com_distance_matrix[
                        :,
                        step_column,
                    ]
                )
            ),
            "Mean_normal_angle_deg": float(
                np.mean(
                    normal_angle_matrix[
                        :,
                        step_column,
                    ]
                )
            ),
        }
    )

summary_df = pd.DataFrame(
    summary_rows
)

summary_df.to_csv(
    SUMMARY_CSV,
    index=False,
    encoding="utf-8-sig",
)

timeseries_data = {
    "Trajectory_index_0based": (
        analysis_frame_indices
    ),
    "Original_frame_1based": (
        analysis_frame_indices + 1
    ),
    "Time_ps": times_ps,
    "Time_ns": times_ps / 1000.0,
}

for step_column, record in enumerate(
    step_records
):
    residue_1 = record["residue_1"]
    residue_2 = record["residue_2"]

    column_name = (
        f"{record['sequence']}"
        f"_{residue_1.resid}-{residue_2.resid}"
    )

    timeseries_data[
        column_name
    ] = stacked_matrix[
        :,
        step_column,
    ]

timeseries_df = pd.DataFrame(
    timeseries_data
)

timeseries_df.to_csv(
    TIMESERIES_CSV,
    index=False,
    encoding="utf-8-sig",
)

WRITE_DIAGNOSTICS = True

if WRITE_DIAGNOSTICS:
    diagnostics_rows = []

    for step_column, record in enumerate(
        step_records
    ):
        residue_1 = record["residue_1"]
        residue_2 = record["residue_2"]

        for analysis_row in range(
            number_of_analysis_frames
        ):
            diagnostics_rows.append(
                {
                    "Trajectory_index_0based": int(
                        analysis_frame_indices[
                            analysis_row
                        ]
                    ),
                    "Original_frame_1based": int(
                        analysis_frame_indices[
                            analysis_row
                        ]
                        + 1
                    ),
                    "Time_ps": float(
                        times_ps[
                            analysis_row
                        ]
                    ),
                    "Time_ns": float(
                        times_ps[
                            analysis_row
                        ]
                        / 1000.0
                    ),
                    "Sequence": record["sequence"],
                    "Step_label": (
                        f"{residue_1.resid}"
                        f"-{residue_2.resid}"
                    ),
                    "Resid_1": int(
                        residue_1.resid
                    ),
                    "Resname_1": str(
                        residue_1.resname
                    ),
                    "Resid_2": int(
                        residue_2.resid
                    ),
                    "Resname_2": str(
                        residue_2.resname
                    ),
                    "Minimum_distance_A": float(
                        minimum_distance_matrix[
                            analysis_row,
                            step_column,
                        ]
                    ),
                    "COM_distance_A": float(
                        com_distance_matrix[
                            analysis_row,
                            step_column,
                        ]
                    ),
                    "Normal_angle_deg": float(
                        normal_angle_matrix[
                            analysis_row,
                            step_column,
                        ]
                    ),
                    "Stacked": int(
                        stacked_matrix[
                            analysis_row,
                            step_column,
                        ]
                    ),
                }
            )

    diagnostics_df = pd.DataFrame(
        diagnostics_rows
    )

    diagnostics_df.to_csv(
        DIAGNOSTICS_CSV,
        index=False,
        encoding="utf-8-sig",
    )

print()
print("Analysis completed")

print()
print(
    summary_df[
        [
            "Sequence",
            "Step_label",
            "Resname_1",
            "Resname_2",
            "Stacked_frames",
            "Total_frames",
            "Stacking_percentage",
        ]
    ].to_string(
        index=False
    )
)

print()
print(f"Output directory : {OUTPUT_DIR}")
print(f"Summary          : {SUMMARY_CSV}")
print(f"Time series      : {TIMESERIES_CSV}")

if WRITE_DIAGNOSTICS:
    print(f"Diagnostics      : {DIAGNOSTICS_CSV}")

