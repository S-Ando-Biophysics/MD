#!/usr/bin/env python3

from __future__ import annotations

import gc
import warnings
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis import align
from matplotlib.ticker import MaxNLocator
from sklearn.cluster import HDBSCAN
from sklearn.decomposition import IncrementalPCA



# ============================================================
# User settings
# ============================================================


# ------------------------------------------------------------
# Residue selections
# ------------------------------------------------------------

FIT_SELECTION = (
    "resid 1:○○ and not name H*"
)

PCA_SELECTION = (
    "resid 1:○○ and not name H*"
)

OUTPUT_SELECTION = "resid 1:○○"


# ------------------------------------------------------------
# Trajectory settings
# ------------------------------------------------------------

# Python uses zero-based indexing. (0 = Frame 1, 1 = Frame 2, ...)

# Frame 1 is excluded here because it is the initial structure.
# If you do not want to exclude it, please set to 0.
START_FRAME_INDEX = 1

# None = analyze through the final frame
STOP_FRAME_INDEX = None

# 1 = all frames
# 10 = every 10th frame
STRIDE = 1

# Reference frame used for structural alignment
# 0 = Frame 1 of the original trajectory
REFERENCE_FRAME_INDEX = 0


# ------------------------------------------------------------
# PCA settings
# ------------------------------------------------------------

# Number of principal components
N_PCS = 10

# Batch size for incremental PCA
PCA_BATCH_SIZE = 1000


# ------------------------------------------------------------
# HDBSCAN settings
# ------------------------------------------------------------

# Minimum number of frames required to form a cluster
MIN_CLUSTER_SIZE = 500

# Controls the strictness of density estimation
MIN_SAMPLES = 5

# "eom"  = extract major clusters
# "leaf" = produce finer clusters
CLUSTER_SELECTION_METHOD = "eom"

# Whether to allow a single cluster
# True = allow
# False = do not allow
ALLOW_SINGLE_CLUSTER = True


# ============================================================





SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent

INPUT_DIR = ROOT_DIR / "03_GROMACS" / "06_OUTPUT"
ANALYSIS_ROOT = ROOT_DIR / "04_ANALYSIS"

ANALYSIS_ROOT.mkdir(parents=True, exist_ok=True)

analysis_numbers = []

for path in ANALYSIS_ROOT.iterdir():
    if not path.is_dir():
        continue

    match = re.match(r"^(\d+)_", path.name)

    if match:
        analysis_numbers.append(
            int(match.group(1))
        )

next_analysis_number = (
    max(analysis_numbers) + 1
    if analysis_numbers
    else 1
)

OUTPUT_DIR = (
    ANALYSIS_ROOT
    / f"{next_analysis_number:02d}_CLUSTER_PCA-HDBSCAN"
)

REPRESENTATIVE_DIR = OUTPUT_DIR / "representative-structures"

COORDINATE_CACHE = OUTPUT_DIR / "aligned_coordinates.float32.dat"
ASSIGNMENT_CSV = OUTPUT_DIR / "frame_assignments.csv"
SUMMARY_CSV = OUTPUT_DIR / "cluster_summary.csv"
PCA_VARIANCE_CSV = OUTPUT_DIR / "PCA_explained_variance.csv"
PCA_SCORES_FILE = OUTPUT_DIR / "PCA_scores.npy"

def find_latest_file(directory: Path, pattern: str) -> Path:
    """Return the most recently modified file matching the pattern."""
    if not directory.is_dir():
        raise FileNotFoundError(
            f"Input directory not found:\n{directory}"
        )

    candidates = [
        path
        for path in directory.glob(pattern)
        if path.is_file()
    ]

    if not candidates:
        raise FileNotFoundError(
            f"No file matching {pattern!r} was found in:\n{directory}"
        )

    return max(
        candidates,
        key=lambda path: path.stat().st_mtime,
    )

TOPOLOGY_FILE = find_latest_file(INPUT_DIR, "*.pdb")
TRAJECTORY_FILE = find_latest_file(INPUT_DIR, "*.xtc")

# ============================================================
# Figure configuration
# ============================================================

FIGURE_WIDTH_MM = 182
FIGURE_HEIGHT_MM = 130
OUTPUT_DPI = 600


def mm_to_inch(value_mm: float) -> float:
    return value_mm / 25.4


FIGURE_SIZE = (
    mm_to_inch(FIGURE_WIDTH_MM),
    mm_to_inch(FIGURE_HEIGHT_MM),
)

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 10,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

# ============================================================
# Preparing
# ============================================================

def check_input_file(file_path: Path) -> None:
    if not file_path.exists():
        raise FileNotFoundError(
            f"Input file not found.\n"
            f"{file_path}"
        )


check_input_file(TOPOLOGY_FILE)
check_input_file(TRAJECTORY_FILE)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
REPRESENTATIVE_DIR.mkdir(parents=True, exist_ok=True)

print("Loading trajectory")
print(f"Topology   : {TOPOLOGY_FILE}")
print(f"Trajectory : {TRAJECTORY_FILE}")

universe = mda.Universe(
    str(TOPOLOGY_FILE),
    str(TRAJECTORY_FILE),
)

number_of_total_frames = len(universe.trajectory)

if number_of_total_frames < 2:
    raise ValueError(
        "The trajectory must contain at least two frames."
    )

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

number_of_analysis_frames = len(
    analysis_frame_indices
)

if number_of_analysis_frames == 0:
    raise ValueError(
        "No frames are available for analysis."
    )

fit_atoms = universe.select_atoms(
    FIT_SELECTION
)

pca_atoms = universe.select_atoms(
    PCA_SELECTION
)

output_atoms = universe.select_atoms(
    OUTPUT_SELECTION
)

if len(fit_atoms) == 0:
    raise ValueError(
        "No atoms were selected by FIT_SELECTION.\n"
        f"Selection: {FIT_SELECTION}\n"
        "Check residue IDs in the PDB file."
    )

if len(pca_atoms) == 0:
    raise ValueError(
        "No atoms were selected by PCA_SELECTION.\n"
        f"Selection: {PCA_SELECTION}"
    )

if len(output_atoms) == 0:
    raise ValueError(
        "No atoms were selected by OUTPUT_SELECTION."
    )

print()
print("Atom selections")
print(f"Fitting atoms : {len(fit_atoms)}")
print(f"PCA atoms     : {len(pca_atoms)}")
print(f"Output atoms  : {len(output_atoms)}")
print(f"Total frames  : {number_of_total_frames}")
print(f"Used frames   : {number_of_analysis_frames}")

universe.trajectory[
    REFERENCE_FRAME_INDEX
]

reference_fit_positions = (
    fit_atoms.positions
    .astype(np.float64)
    .copy()
)

reference_center = np.mean(
    reference_fit_positions,
    axis=0,
)

reference_fit_centered = (
    reference_fit_positions
    - reference_center
)

# ============================================================
# Align
# ============================================================

def align_current_frame() -> float:

    mobile_positions = (
        fit_atoms.positions
        .astype(np.float64)
        .copy()
    )

    mobile_center = np.mean(
        mobile_positions,
        axis=0,
    )

    mobile_centered = (
        mobile_positions
        - mobile_center
    )

    rotation_matrix, _ = align.rotation_matrix(
        mobile_centered,
        reference_fit_centered,
    )

    # Apply the same transformation to all atoms
    universe.atoms.translate(
        -mobile_center
    )

    universe.atoms.rotate(
        rotation_matrix
    )

    universe.atoms.translate(
        reference_center
    )

    difference = (
        fit_atoms.positions.astype(np.float64)
        - reference_fit_positions
    )

    fitted_rmsd = np.sqrt(
        np.mean(
            np.sum(
                difference**2,
                axis=1,
            )
        )
    )

    return float(fitted_rmsd)

# ============================================================
# Extract coordinates for PCA
# ============================================================

number_of_features = (
    len(pca_atoms) * 3
)

if N_PCS > min(
    number_of_analysis_frames,
    number_of_features,
):
    raise ValueError(
        f"N_PCS={N_PCS} is too large.\n"
        f"Number of analysis frames: {number_of_analysis_frames}\n"
        f"Number of coordinate features: {number_of_features}"
    )

coordinates = np.memmap(
    COORDINATE_CACHE,
    dtype=np.float32,
    mode="w+",
    shape=(
        number_of_analysis_frames,
        number_of_features,
    ),
)

times_ps = np.empty(
    number_of_analysis_frames,
    dtype=np.float64,
)

fitting_rmsd = np.empty(
    number_of_analysis_frames,
    dtype=np.float64,
)

print()
print("Extracting aligned coordinates")

for analysis_row, trajectory_index in enumerate(
    analysis_frame_indices
):
    timestep = universe.trajectory[
        trajectory_index
    ]

    fitting_rmsd[analysis_row] = (
        align_current_frame()
    )

    coordinates[
        analysis_row,
        :,
    ] = (
        pca_atoms.positions
        .astype(np.float32)
        .reshape(-1)
    )

    times_ps[analysis_row] = float(
        timestep.time
    )

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

coordinates.flush()

# ============================================================
# Incremental PCA
# ============================================================

print()
print("Running PCA")
print(f"Number of PCs : {N_PCS}")
print(f"Batch size    : {PCA_BATCH_SIZE}")

pca_model = IncrementalPCA(
    n_components=N_PCS,
    batch_size=PCA_BATCH_SIZE,
)

pca_model.fit(
    coordinates
)

pca_scores = np.empty(
    (
        number_of_analysis_frames,
        N_PCS,
    ),
    dtype=np.float32,
)

for batch_start in range(
    0,
    number_of_analysis_frames,
    PCA_BATCH_SIZE,
):
    batch_end = min(
        batch_start + PCA_BATCH_SIZE,
        number_of_analysis_frames,
    )

    pca_scores[
        batch_start:batch_end
    ] = pca_model.transform(
        coordinates[
            batch_start:batch_end
        ]
    ).astype(np.float32)

np.save(
    PCA_SCORES_FILE,
    pca_scores,
)

explained_variance_percent = (
    pca_model.explained_variance_ratio_
    * 100.0
)

cumulative_variance_percent = np.cumsum(
    explained_variance_percent
)

variance_df = pd.DataFrame(
    {
        "PC": np.arange(1, N_PCS + 1),
        "Explained_variance_percent": (
            explained_variance_percent
        ),
        "Cumulative_variance_percent": (
            cumulative_variance_percent
        ),
    }
)

variance_df.to_csv(
    PCA_VARIANCE_CSV,
    index=False,
    encoding="utf-8-sig",
)

print(
    f"PC1-PC{N_PCS} cumulative variance: "
    f"{cumulative_variance_percent[-1]:.2f}%"
)

# ============================================================
# Temporary coordinate cache setting
# ============================================================

# False = delete the aligned-coordinate cache after PCA.
# True  = keep the cache file.
KEEP_COORDINATE_CACHE = False

# ============================================================
# Close and remove temporary coordinate cache
# ============================================================

del coordinates
gc.collect()

if (
    not KEEP_COORDINATE_CACHE
    and COORDINATE_CACHE.exists()
):
    try:
        COORDINATE_CACHE.unlink()
    except PermissionError:
        warnings.warn(
            "Could not delete the temporary coordinate cache.\n"
            f"{COORDINATE_CACHE}",
            stacklevel=2,
        )

# ============================================================
# HDBSCAN execution setting
# ============================================================

# Number of CPU cores. -1 = use all available cores.
N_JOBS = -1

# ============================================================
# HDBSCAN
# ============================================================

if MIN_CLUSTER_SIZE >= number_of_analysis_frames:
    raise ValueError(
        "MIN_CLUSTER_SIZE must be smaller than the number of analysis frames."
    )

print()
print("Running HDBSCAN")
print(
    f"min_cluster_size       : "
    f"{MIN_CLUSTER_SIZE}"
)
print(
    f"min_samples            : "
    f"{MIN_SAMPLES}"
)
print(
    f"cluster_selection      : "
    f"{CLUSTER_SELECTION_METHOD}"
)

cluster_model = HDBSCAN(
    min_cluster_size=MIN_CLUSTER_SIZE,
    min_samples=MIN_SAMPLES,
    metric="euclidean",
    algorithm="auto",
    cluster_selection_method=(
        CLUSTER_SELECTION_METHOD
    ),
    allow_single_cluster=ALLOW_SINGLE_CLUSTER,
    store_centers="medoid",
    n_jobs=N_JOBS,
)

raw_labels = cluster_model.fit_predict(
    pca_scores
)

membership_probabilities = (
    cluster_model.probabilities_
)

raw_cluster_labels = np.sort(
    np.unique(
        raw_labels[
            raw_labels >= 0
        ]
    )
)

number_of_clusters = len(
    raw_cluster_labels
)

noise_mask = raw_labels == -1
number_of_noise_frames = int(
    np.sum(noise_mask)
)

if number_of_clusters == 0:
    raise RuntimeError(
        "No clusters were detected.\n"
        "Reduce MIN_CLUSTER_SIZE or MIN_SAMPLES."
    )

# ============================================================
# Renumber clusters as states in descending population order
# ============================================================

cluster_sizes = {
    int(raw_label): int(
        np.sum(raw_labels == raw_label)
    )
    for raw_label in raw_cluster_labels
}

population_order = sorted(
    raw_cluster_labels,
    key=lambda raw_label: (
        -cluster_sizes[int(raw_label)],
        int(raw_label),
    ),
)

raw_to_state = {
    int(raw_label): state_number
    for state_number, raw_label
    in enumerate(
        population_order,
        start=1,
    )
}

state_numbers = np.zeros(
    number_of_analysis_frames,
    dtype=int,
)

for raw_label, state_number in (
    raw_to_state.items()
):
    state_numbers[
        raw_labels == raw_label
    ] = state_number

# ============================================================
# Map HDBSCAN medoids to actual trajectory frames
# ============================================================

medoid_coordinates = np.asarray(
    cluster_model.medoids_,
    dtype=np.float64,
)

representative_index_by_raw_label: dict[
    int,
    int
] = {}

for medoid in medoid_coordinates:
    squared_distances = np.sum(
        (
            pca_scores.astype(np.float64)
            - medoid
        ) ** 2,
        axis=1,
    )

    representative_index = int(
        np.argmin(squared_distances)
    )

    corresponding_label = int(
        raw_labels[
            representative_index
        ]
    )

    if corresponding_label >= 0:
        representative_index_by_raw_label[
            corresponding_label
        ] = representative_index

for raw_label in raw_cluster_labels:
    raw_label = int(raw_label)

    if (
        raw_label
        in representative_index_by_raw_label
    ):
        continue

    cluster_indices = np.flatnonzero(
        raw_labels == raw_label
    )

    cluster_center = np.mean(
        pca_scores[
            cluster_indices
        ],
        axis=0,
    )

    distances = np.sum(
        (
            pca_scores[
                cluster_indices
            ]
            - cluster_center
        ) ** 2,
        axis=1,
    )

    representative_index_by_raw_label[
        raw_label
    ] = int(
        cluster_indices[
            np.argmin(distances)
        ]
    )

# ============================================================
# Cluster statistics
# ============================================================

summary_rows = []

representative_flags = np.zeros(
    number_of_analysis_frames,
    dtype=bool,
)

for raw_label in population_order:
    raw_label = int(raw_label)

    state_number = raw_to_state[
        raw_label
    ]

    cluster_indices = np.flatnonzero(
        raw_labels == raw_label
    )

    representative_index = (
        representative_index_by_raw_label[
            raw_label
        ]
    )

    representative_flags[
        representative_index
    ] = True

    trajectory_index = int(
        analysis_frame_indices[
            representative_index
        ]
    )

    original_frame_number = (
        trajectory_index + 1
    )

    representative_time_ps = float(
        times_ps[
            representative_index
        ]
    )

    cluster_size = len(
        cluster_indices
    )

    occupancy_percent = (
        cluster_size
        / number_of_analysis_frames
        * 100.0
    )

    summary_rows.append(
        {
            "State": state_number,
            "Raw_cluster_label": raw_label,
            "Frames": cluster_size,
            "Occupancy_percent": occupancy_percent,
            "Mean_membership_probability": float(
                np.mean(
                    membership_probabilities[
                        cluster_indices
                    ]
                )
            ),
            "First_time_ns": float(
                np.min(
                    times_ps[
                        cluster_indices
                    ]
                )
                / 1000.0
            ),
            "Last_time_ns": float(
                np.max(
                    times_ps[
                        cluster_indices
                    ]
                )
                / 1000.0
            ),
            "Representative_analysis_row_0based": (
                representative_index
            ),
            "Representative_trajectory_index_0based": (
                trajectory_index
            ),
            "Representative_frame_1based": (
                original_frame_number
            ),
            "Representative_time_ps": (
                representative_time_ps
            ),
            "Representative_time_ns": (
                representative_time_ps
                / 1000.0
            ),
            "Representative_membership_probability": float(
                membership_probabilities[
                    representative_index
                ]
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

# ============================================================
# Frame assignments
# ============================================================

assignment_data = {
    "Analysis_row_1based": np.arange(
        1,
        number_of_analysis_frames + 1,
    ),
    "Trajectory_index_0based": (
        analysis_frame_indices
    ),
    "Original_frame_1based": (
        analysis_frame_indices + 1
    ),
    "Time_ps": times_ps,
    "Time_ns": times_ps / 1000.0,
    "Fitting_RMSD_A": fitting_rmsd,
    "Raw_cluster_label": raw_labels,
    "State": state_numbers,
    "Membership_probability": (
        membership_probabilities
    ),
    "Is_noise": noise_mask,
    "Is_representative": (
        representative_flags
    ),
}

for pc_index in range(N_PCS):
    assignment_data[
        f"PC{pc_index + 1}"
    ] = pca_scores[:, pc_index]

assignment_df = pd.DataFrame(
    assignment_data
)

assignment_df.to_csv(
    ASSIGNMENT_CSV,
    index=False,
    encoding="utf-8-sig",
)

# ============================================================
# Write representative structures as PDB files
# ============================================================

print()
print("Writing representative structures")

for _, row in summary_df.iterrows():
    state_number = int(
        row["State"]
    )

    raw_label = int(
        row["Raw_cluster_label"]
    )

    trajectory_index = int(
        row[
            "Representative_trajectory_index_0based"
        ]
    )

    original_frame = int(
        row["Representative_frame_1based"]
    )

    time_ns = float(
        row["Representative_time_ns"]
    )

    universe.trajectory[
        trajectory_index
    ]

    align_current_frame()

    output_file = (
        REPRESENTATIVE_DIR
        / (
            f"State_{state_number:02d}"
            f"_cluster_{raw_label}"
            f"_frame_{original_frame}"
            f"_{time_ns:.3f}ns.pdb"
        )
    )

    output_atoms.write(
        str(output_file)
    )

    print(output_file.name)

# ============================================================
# Color settings
# ============================================================

cluster_colormap = plt.get_cmap(
    "tab20"
)

state_colors = {
    state_number: cluster_colormap(
        (state_number - 1) % 20
    )
    for state_number in range(
        1,
        number_of_clusters + 1,
    )
}

# ============================================================
# PC1-PC2 plot
# ============================================================

fig, ax = plt.subplots(
    figsize=FIGURE_SIZE
)

if np.any(noise_mask):
    ax.scatter(
        pca_scores[noise_mask, 0],
        pca_scores[noise_mask, 1],
        s=2,
        marker=".",
        color="lightgray",
        alpha=0.35,
        linewidths=0,
        rasterized=True,
        label="Noise",
    )

for state_number in range(
    1,
    number_of_clusters + 1,
):
    state_mask = (
        state_numbers == state_number
    )

    ax.scatter(
        pca_scores[state_mask, 0],
        pca_scores[state_mask, 1],
        s=3,
        marker=".",
        color=state_colors[
            state_number
        ],
        alpha=0.65,
        linewidths=0,
        rasterized=True,
        label=f"State {state_number}",
    )

for _, row in summary_df.iterrows():
    state_number = int(
        row["State"]
    )

    representative_index = int(
        row[
            "Representative_analysis_row_0based"
        ]
    )

    x_value = pca_scores[
        representative_index,
        0,
    ]

    y_value = pca_scores[
        representative_index,
        1,
    ]

    ax.scatter(
        x_value,
        y_value,
        s=50,
        marker="x",
        color="black",
        linewidths=1.2,
        zorder=10,
    )

    ax.annotate(
        f"S{state_number}",
        xy=(x_value, y_value),
        xytext=(4, 4),
        textcoords="offset points",
        fontweight="bold",
    )

ax.set_xlabel(
    f"PC1 ({explained_variance_percent[0]:.1f}%)"
)

ax.set_ylabel(
    f"PC2 ({explained_variance_percent[1]:.1f}%)"
)

ax.set_title(
    "PCA-HDBSCAN clustering",
    loc="left",
    fontweight="bold",
)

ax.tick_params(
    axis="both",
    direction="out",
)

if number_of_clusters <= 12:
    ax.legend(
        frameon=False,
        markerscale=3,
    )

fig.tight_layout()

fig.savefig(
    OUTPUT_DIR / "PCA_HDBSCAN_PC1_PC2.png",
    dpi=OUTPUT_DPI,
)

fig.savefig(
    OUTPUT_DIR / "PCA_HDBSCAN_PC1_PC2.pdf",
)

plt.close(fig)

# ============================================================
# Cluster assignment timeline
# ============================================================

fig, ax = plt.subplots(
    figsize=FIGURE_SIZE
)

if np.any(noise_mask):
    ax.scatter(
        times_ps[noise_mask] / 1000.0,
        np.zeros(
            number_of_noise_frames
        ),
        s=2,
        marker=".",
        color="lightgray",
        linewidths=0,
        rasterized=True,
    )

for state_number in range(
    1,
    number_of_clusters + 1,
):
    state_mask = (
        state_numbers == state_number
    )

    ax.scatter(
        times_ps[state_mask] / 1000.0,
        np.full(
            np.sum(state_mask),
            state_number,
        ),
        s=3,
        marker=".",
        color=state_colors[
            state_number
        ],
        linewidths=0,
        rasterized=True,
    )

ax.set_xlabel("Time (ns)")
ax.set_ylabel("State")

ax.set_yticks(
    np.arange(
        0,
        number_of_clusters + 1,
    )
)

ax.set_yticklabels(
    ["Noise"]
    + [
        str(state_number)
        for state_number in range(
            1,
            number_of_clusters + 1,
        )
    ]
)

ax.set_title(
    "Cluster assignment over time",
    loc="left",
    fontweight="bold",
)

ax.tick_params(
    axis="both",
    direction="out",
)

fig.tight_layout()

fig.savefig(
    OUTPUT_DIR / "cluster_timeline.png",
    dpi=OUTPUT_DPI,
)

fig.savefig(
    OUTPUT_DIR / "cluster_timeline.pdf",
)

plt.close(fig)

# ============================================================
# Cluster populations
# ============================================================

fig, ax = plt.subplots(
    figsize=FIGURE_SIZE
)

states = summary_df[
    "State"
].to_numpy(dtype=int)

occupancies = summary_df[
    "Occupancy_percent"
].to_numpy(dtype=float)

bar_colors = [
    state_colors[int(state)]
    for state in states
]

bars = ax.bar(
    states,
    occupancies,
    color=bar_colors,
    edgecolor="black",
    linewidth=0.6,
)

for bar, occupancy in zip(
    bars,
    occupancies,
):
    ax.text(
        bar.get_x()
        + bar.get_width() / 2.0,
        bar.get_height(),
        f"{occupancy:.1f}",
        ha="center",
        va="bottom",
    )

ax.set_xlabel("State")
ax.set_ylabel("Occupancy (%)")

ax.set_xticks(states)

ax.set_title(
    "Cluster populations",
    loc="left",
    fontweight="bold",
)

ax.tick_params(
    axis="both",
    direction="out",
)

fig.tight_layout()

fig.savefig(
    OUTPUT_DIR / "cluster_populations.png",
    dpi=OUTPUT_DPI,
)

fig.savefig(
    OUTPUT_DIR / "cluster_populations.pdf",
)

plt.close(fig)

# ============================================================
# PCA explained variance
# ============================================================

fig, ax = plt.subplots(
    figsize=FIGURE_SIZE
)

pc_numbers = np.arange(
    1,
    N_PCS + 1,
)

ax.plot(
    pc_numbers,
    cumulative_variance_percent,
    marker="o",
)

ax.set_xlabel("Principal component")
ax.set_ylabel(
    "Cumulative explained variance (%)"
)

ax.set_xticks(
    pc_numbers
)

ax.yaxis.set_major_locator(
    MaxNLocator(integer=True)
)

ax.set_ylim(
    0,
    min(
        105,
        max(
            10,
            cumulative_variance_percent[-1]
            + 5,
        ),
    ),
)

ax.set_title(
    "PCA cumulative explained variance",
    loc="left",
    fontweight="bold",
)

ax.tick_params(
    axis="both",
    direction="out",
)

fig.tight_layout()

fig.savefig(
    OUTPUT_DIR / "PCA_cumulative_variance.png",
    dpi=OUTPUT_DPI,
)

fig.savefig(
    OUTPUT_DIR / "PCA_cumulative_variance.pdf",
)

plt.close(fig)

# ============================================================
# Final summary
# ============================================================

noise_percent = (
    number_of_noise_frames
    / number_of_analysis_frames
    * 100.0
)

print()
print("Analysis completed")
print(
    f"Clusters : {number_of_clusters}"
)

print(
    f"Noise    : "
    f"{number_of_noise_frames}"
    f"/{number_of_analysis_frames} "
    f"({noise_percent:.2f}%)"
)

print()
print(
    summary_df[
        [
            "State",
            "Frames",
            "Occupancy_percent",
            "Representative_frame_1based",
            "Representative_time_ns",
        ]
    ].to_string(index=False)
)

print()
print(f"Output directory : {OUTPUT_DIR}")
print(f"Assignments      : {ASSIGNMENT_CSV}")
print(f"Cluster summary  : {SUMMARY_CSV}")
print(
    f"Representatives  : "
    f"{REPRESENTATIVE_DIR}"
)
