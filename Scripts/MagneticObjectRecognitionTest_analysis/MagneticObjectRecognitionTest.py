# -*- coding: utf-8 -*-
"""Analyze magnetic object exploration from single-animal SLEAP HDF5 tracks.

Edit the paths and experiment configuration below, then run this file.
The default analysis uses the first 60 seconds at 16 frames/s and a 100-pixel
object radius. Missing coordinates are retained without interpolation.

For each window, the script saves a skeleton-summary PDF and a JSON file of
zero-based approach-frame indices for each recording, plus a dated CSV table.
DI_window_B_A = (time_B - time_A) / (time_B + time_A); DI_Magnet is oriented
so positive values indicate more exploration of the object with magnets.
The legacy name "novel_side" refers to the magnet side in this experiment.

Dependencies: h5py, numpy, matplotlib, pandas.
Author: Li Zhang. Original version: 2025-09-02.
"""

import json
from datetime import datetime
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ----------------------------
# 1. Paths and analysis settings
# ----------------------------
subset_name = "MOR_pooled"  # Magnetic object tests.
windows = [60]  # fixed window analysis: 60 seconds

h5_folder = Path(
    r"X:\NMR_group_data\_Lab_Data\Manuscripts\MouseMagnetoreception\Revision_ScienceAdvances\Scripts\MagneticObjectRecognitionTest_analysis\h5_outputs_from_SLEAP"
)
result_folder = Path(
    rf"X:\NMR_group_data\_Lab_Data\Manuscripts\MouseMagnetoreception\Revision_ScienceAdvances\Scripts\MagneticObjectRecognitionTest_analysis\results_{subset_name}"
)

fps = 16  # Acquisition frame rate, in frames per second.
radius_px = 100  # Distance from each object center, in pixels.

# 2. Recording metadata
# Coordinates are (x, y) in image pixels. All original values are retained.
# Filenames follow the renamed format Mouse_<mouse_id>.h5.
# side_list (magnet side) and group_list follow the entry order in
# object_positions. Keep all three in the same order when editing metadata.
# Missing or unreadable files do not change the assignments for other recordings.
experiment_config = {
    'MOR_pooled': {
        'side_list': [
            "A", "B", "B", "A", "A", "B", "B", "A", "A", "B", "B", "A", "A", "B", "B", "A",
           
            "A", "B", "B", "A", "B", "A", "A", "B", "B", "A", "A", "B", "A", "B", "B",

            "A", "B", "B", "A", "A", "B","A", "B", "B", "A", "A", "B", "B", "A", "A",

            "A", "B", "B", "A", "A", "B", "B", "A", "A", "B", "B", "A",

            "A", "B", "B", "A", "A", "B", "B", "A", "A",

            "B", "A", "B", "B", "A", "B", "A", "A",

            "A", "B", "B", "A", "B", "A", "B", "B", "A",

            "B", "A", "A", "B", "B", "B", "A", "A", "B",

            "A", "B", "A", "B", "B", "B", "A", "B", "A", "A",

            "B", "A", "A", "B", "B", "A", "A", "B", "B", "A", "B", "A",
        ],

        'group_list': [
            "hxM", "hxM", "hxM", "hxM", "hxM", "hxM", "hxM", "hxM", "hxF", "hxF", "hxF", "hxF", "hxF", "hxF", "hxF", "hxF",
            
            "koM", "koM", "koM", "koM", "koM", "koM", "hxM", "hxM", "hxM", "hxM", "hxM", "hxM", "koM", "koM", "koM",

            "koF", "koF", "koF", "koF", "koF", "koF", "hxF", "hxF", "hxF", "hxF", "hxF", "hxF", "koF", "koF", "koF",

            "koF", "koF", "koF", "koF", "koF", "koF", "koM", "koM", "koM", "koM", "koM", "koM",

            "hxM_3", "hxM_3", "hxM_3", "hxM_3", "hxM_5", "hxM_5", "hxM_5", "hxM_5", "hxM_5",

            "hxF_5", "hxF_3", "hxF_3", "hxF_5", "hxF_5", "hxF_3", "hxF_3", "hxF_3",

            "hxM_3", "hxM_3", "hxM_3", "hxM_5", "hxM_5", "hxM_3", "hxM_3", "hxM_5", "hxM_5",

            "hxF_3", "hxF_5", "hxF_5", "hxF_5", "hxF_5", "hxF_3", "hxF_3", "hxF_5", "hxF_5",

            "hxM_5", "hxM_5", "hxM_3", "hxM_3", "hxM_3", "hxM_5", "hxM_5", "hxM_5", "hxM_3", "hxM_3",

            "hxF_3", "hxF_3", "hxF_3", "hxF_3", "hxF_5", "hxF_5", "hxF_5", "hxF_5", "hxF_3", "hxF_3", "hxF_5", "hxF_5",
        ],

        'object_positions': {
            "Mouse_54809.h5": {"A": (413, 513), "B": (858, 514)},
            "Mouse_54965.h5": {"A": (420, 520), "B": (862, 514)},
            "Mouse_54966.h5": {"A": (418, 526), "B": (869, 525)},
            "Mouse_54967.h5": {"A": (421, 526), "B": (860, 528)},
            "Mouse_54969.h5": {"A": (408, 532), "B": (860, 534)},
            "Mouse_54999.h5": {"A": (418, 540), "B": (860, 533)},
            "Mouse_55000.h5": {"A": (417, 553), "B": (860, 528)},
            "Mouse_55001.h5": {"A": (420, 544), "B": (858, 538)},
            "Mouse_54951.h5": {"A": (417, 544), "B": (864, 534)},
            "Mouse_54945.h5": {"A": (420, 540), "B": (864, 538)},
            "Mouse_54947.h5": {"A": (410, 546), "B": (864, 541)},
            "Mouse_54958.h5": {"A": (413, 553), "B": (865, 550)},
            "Mouse_54959.h5": {"A": (416, 553), "B": (858, 548)},
            "Mouse_54972.h5": {"A": (416, 553), "B": (857, 545)},
            "Mouse_54974.h5": {"A": (416, 548), "B": (861, 544)},
            "Mouse_54961.h5": {"A": (418, 549), "B": (874, 537)},
            "Mouse_54829.h5": {"A": (417, 550), "B": (864, 554)},
            "Mouse_54828.h5": {"A": (416, 549), "B": (860, 541)},
            "Mouse_54830.h5": {"A": (413, 560), "B": (862, 545)},
            "Mouse_54831.h5": {"A": (414, 556), "B": (866, 542)},
            "Mouse_54975.h5": {"A": (420, 558), "B": (860, 540)},
            "Mouse_55009.h5": {"A": (412, 553), "B": (862, 548)},
            "Mouse_55007.h5": {"A": (412, 569), "B": (857, 568)},
            "Mouse_55078.h5": {"A": (414, 566), "B": (856, 562)},
            "Mouse_55079.h5": {"A": (414, 569), "B": (858, 566)},
            "Mouse_55080.h5": {"A": (414, 566), "B": (860, 554)},
            "Mouse_55081.h5": {"A": (416, 568), "B": (856, 564)},
            "Mouse_55008.h5": {"A": (418, 570), "B": (862, 558)},
            "Mouse_55010.h5": {"A": (413, 572), "B": (862, 566)},
            "Mouse_55056.h5": {"A": (412, 557), "B": (857, 553)},
            "Mouse_55073.h5": {"A": (413, 558), "B": (857, 552)},
            "Mouse_54811.h5": {"A": (413, 520), "B": (868, 512)},
            "Mouse_54832.h5": {"A": (410, 522), "B": (865, 517)},
            "Mouse_54949.h5": {"A": (409, 525), "B": (870, 513)},
            "Mouse_54973.h5": {"A": (418, 514), "B": (874, 516)},
            "Mouse_54963.h5": {"A": (412, 522), "B": (868, 509)},
            "Mouse_54964.h5": {"A": (413, 510), "B": (870, 510)},
            "Mouse_55002.h5": {"A": (410, 508), "B": (864, 514)},
            "Mouse_55011.h5": {"A": (402, 512), "B": (864, 506)},
            "Mouse_55087.h5": {"A": (409, 513), "B": (866, 505)},
            "Mouse_55084.h5": {"A": (406, 518), "B": (868, 517)},
            "Mouse_55085.h5": {"A": (406, 510), "B": (868, 508)},
            "Mouse_55058.h5": {"A": (401, 509), "B": (868, 501)},
            "Mouse_55083.h5": {"A": (409, 512), "B": (866, 508)},
            "Mouse_55013.h5": {"A": (408, 516), "B": (865, 502)},
            "Mouse_55014.h5": {"A": (413, 514), "B": (870, 505)},
            "Mouse_54960.h5": {"A": (417, 514), "B": (878, 510)},
            "Mouse_55003.h5": {"A": (417, 518), "B": (878, 512)},
            "Mouse_55069.h5": {"A": (417, 525), "B": (880, 518)},
            "Mouse_55076.h5": {"A": (418, 524), "B": (876, 520)},
            "Mouse_55244.h5": {"A": (418, 522), "B": (878, 513)},
            "Mouse_55246.h5": {"A": (422, 512), "B": (881, 518)},
            "Mouse_55062.h5": {"A": (426, 513), "B": (880, 520)},
            "Mouse_55066.h5": {"A": (422, 510), "B": (878, 513)},
            "Mouse_55153.h5": {"A": (422, 520), "B": (878, 517)},
            "Mouse_55173.h5": {"A": (421, 517), "B": (880, 524)},
            "Mouse_55174.h5": {"A": (421, 524), "B": (880, 525)},
            "Mouse_55248.h5": {"A": (416, 538), "B": (876, 544)},
            "Mouse_55005.h5": {"A": (409, 537), "B": (876, 537)},
            "Mouse_55006.h5": {"A": (409, 530), "B": (874, 541)},
            "Mouse_55057.h5": {"A": (414, 536), "B": (870, 541)},
            "Mouse_55054.h5": {"A": (412, 540), "B": (869, 542)},
            "Mouse_55055.h5": {"A": (413, 540), "B": (866, 542)},
            "Mouse_55072.h5": {"A": (413, 538), "B": (876, 542)},
            "Mouse_55063.h5": {"A": (417, 538), "B": (869, 537)},
            "Mouse_55065.h5": {"A": (414, 540), "B": (872, 540)},
            "Mouse_55155.h5": {"A": (418, 533), "B": (873, 544)},
            "Mouse_55061.h5": {"A": (413, 533), "B": (869, 525)},
            "Mouse_55067.h5": {"A": (416, 530), "B": (872, 526)},
            "Mouse_55068.h5": {"A": (416, 538), "B": (873, 526)},
            "Mouse_55070.h5": {"A": (416, 537), "B": (872, 529)},
            "Mouse_55071.h5": {"A": (409, 538), "B": (872, 528)},
            "Mouse_55156.h5": {"A": (417, 537), "B": (873, 528)},
            "Mouse_55157.h5": {"A": (412, 538), "B": (870, 528)},
            "Mouse_55175.h5": {"A": (409, 540), "B": (872, 526)},
            "Mouse_55247.h5": {"A": (414, 545), "B": (874, 541)},
            "Mouse_55249.h5": {"A": (417, 537), "B": (877, 540)},
            "Mouse_55250.h5": {"A": (416, 541), "B": (874, 540)},
            "Mouse_55251.h5": {"A": (413, 545), "B": (872, 545)},
            "Mouse_55252.h5": {"A": (418, 542), "B": (873, 545)},
            "Mouse_55276.h5": {"A": (413, 538), "B": (874, 544)},
            "Mouse_55277.h5": {"A": (414, 543), "B": (874, 541)},
            "Mouse_55279.h5": {"A": (416, 545), "B": (877, 538)},
            "Mouse_55280.h5": {"A": (418, 549), "B": (874, 536)},
            "Mouse_55255.h5": {"A": (413, 545), "B": (873, 545)},
            "Mouse_55242.h5": {"A": (416, 544), "B": (873, 540)},
            "Mouse_55243.h5": {"A": (416, 542), "B": (872, 536)},
            "Mouse_55268.h5": {"A": (410, 541), "B": (872, 545)},
            "Mouse_55269.h5": {"A": (414, 542), "B": (874, 542)},
            "Mouse_55283.h5": {"A": (413, 549), "B": (872, 544)},
            "Mouse_55284.h5": {"A": (416, 546), "B": (872, 537)},
            "Mouse_55281.h5": {"A": (417, 540), "B": (874, 538)},
            "Mouse_55282.h5": {"A": (414, 542), "B": (872, 541)},
            "Mouse_55304.h5": {"A": (414, 536), "B": (872, 538)},
            "Mouse_55306.h5": {"A": (414, 541), "B": (872, 538)},
            "Mouse_55261.h5": {"A": (412, 544), "B": (870, 538)},
            "Mouse_55262.h5": {"A": (410, 542), "B": (868, 541)},
            "Mouse_55264.h5": {"A": (412, 540), "B": (870, 538)},
            "Mouse_55355.h5": {"A": (413, 542), "B": (873, 541)},
            "Mouse_55356.h5": {"A": (416, 545), "B": (870, 542)},
            "Mouse_55357.h5": {"A": (416, 545), "B": (874, 541)},
            "Mouse_55312.h5": {"A": (420, 549), "B": (870, 542)},
            "Mouse_55344.h5": {"A": (413, 546), "B": (870, 545)},
            "Mouse_55319.h5": {"A": (417, 544), "B": (868, 542)},
            "Mouse_55338.h5": {"A": (414, 545), "B": (872, 544)},
            "Mouse_55340.h5": {"A": (414, 541), "B": (870, 540)},
            "Mouse_55341.h5": {"A": (414, 545), "B": (873, 554)},
            "Mouse_55359.h5": {"A": (412, 541), "B": (870, 542)},
            "Mouse_55360.h5": {"A": (414, 544), "B": (869, 541)},
            "Mouse_55308.h5": {"A": (414, 545), "B": (868, 538)},
            "Mouse_55309.h5": {"A": (414, 540), "B": (870, 538)},
            "Mouse_55349.h5": {"A": (410, 536), "B": (877, 540)},
            "Mouse_55362.h5": {"A": (412, 545), "B": (872, 542)},
            "Mouse_55317.h5": {"A": (416, 541), "B": (876, 537)},
            "Mouse_55318.h5": {"A": (416, 546), "B": (873, 538)},

        },
    },
}

# Get current experiment configuration
side_list = experiment_config[subset_name]["side_list"]
group_list = experiment_config[subset_name]["group_list"]
object_positions = experiment_config[subset_name]["object_positions"]


# 3. SLEAP node names and plotting connections
head_nodes = ["nose", "neck", "earL", "earR"]
tail_nodes = ["tailstart", "tailmid", "tailend"]

# Skeleton connectivity (pairs of nodes)
skeleton_edges = [
    ("nose", "earL"),
    ("nose", "earR"),
    ("nose", "neck"),
    ("neck", "tailstart"),
    ("tailstart", "tailmid"),
    ("tailmid", "tailend"),
]


# ----------------------------
# 4. Tracking input and approach classification
# ----------------------------
def extract_all_nodes(h5_path):
    """Read the first tracked animal from a SLEAP analysis export.

    Expected tracks shape: (animals, 2, nodes, frames).
    Returns coordinates shaped (frames, nodes, 2), plus decoded node names.
    Missing coordinates remain NaN; no interpolation or filtering is applied.
    """
    with h5py.File(h5_path, "r") as f:
        tracks = f["tracks"][:]  # (1, 2, nodes, frames)
        node_names = [n.decode("utf-8") for n in f["node_names"][:]]
        coords = tracks[0].transpose(2, 1, 0)  # (frames, nodes, 2)
    return coords, node_names


def is_approach(frame_coords, node_names, obj_center, radius):
    """Apply the four original geometric approach criteria to one frame.

    At least three head nodes must be inside the radius; the nose and neck
    must face the object according to their distances from its center; at
    least one tail node must be outside the radius; and the nose must be
    more than 30 pixels from the center and closer than the neck.
    All comparisons are strict. NaNs follow the original NumPy comparisons.
    """
    dists = {
        node: np.linalg.norm(frame_coords[node_names.index(node)] - np.array(obj_center))
        for node in node_names
    }

    head_in_radius = len([n for n in head_nodes if dists[n] < radius]) >= 3
    head_facing_object = (dists["nose"] < dists["earL"] and dists["nose"] < dists["earR"]) and (
        dists["neck"] < dists["earL"] or dists["neck"] < dists["earR"]
    )
    tail_outside_radius = any(dists[n] > radius for n in tail_nodes)
    nose_clear_and_forward = dists["nose"] > 30 and dists["neck"] > dists["nose"]

    return head_in_radius and head_facing_object and tail_outside_radius and nose_clear_and_forward


def exploration_time_fixed_window(coords, node_names, objA, objB, radius, fps, window=300):
    """Measure approach time and contiguous episodes within the initial window.

    Analyze min(recording length, int(window * fps)) frames. Time and episode
    durations are in seconds; each qualifying frame contributes 1/fps.
    An unfinished episode is included when the window or recording ends.
    The helper retains its original 300-second default; main() explicitly
    passes the configured window (60 seconds by default).

    Returns (time_A, time_B, DI_B_minus_A, frames_A, frames_B,
             episode_durations_A, episode_durations_B).
    Frame indices are zero based. DI is NaN if neither object is explored.
    """
    tA_window = tB_window = 0
    approach_frames_A, approach_frames_B = [], []
    approach_durations_A, approach_durations_B = [], []

    max_frames = int(window * fps)
    n_frames = min(coords.shape[0], max_frames)

    # Track ongoing approach episodes
    in_approach_A = False
    in_approach_B = False
    current_len_A = 0
    current_len_B = 0

    for frame in range(n_frames):
        sec = 1 / fps
        fcoords = coords[frame]

        # --- Object A ---
        if is_approach(fcoords, node_names, objA, radius):
            tA_window += sec
            approach_frames_A.append(frame)
            current_len_A += sec
            in_approach_A = True
        elif in_approach_A:  # episode ended
            approach_durations_A.append(current_len_A)
            current_len_A = 0
            in_approach_A = False

        # --- Object B ---
        if is_approach(fcoords, node_names, objB, radius):
            tB_window += sec
            approach_frames_B.append(frame)
            current_len_B += sec
            in_approach_B = True
        elif in_approach_B:  # episode ended
            approach_durations_B.append(current_len_B)
            current_len_B = 0
            in_approach_B = False

    # Handle if still in approach at the end of window
    if in_approach_A and current_len_A > 0:
        approach_durations_A.append(current_len_A)
    if in_approach_B and current_len_B > 0:
        approach_durations_B.append(current_len_B)

    di_window = (
        (tB_window - tA_window) / (tB_window + tA_window) if (tA_window + tB_window) > 0 else np.nan
    )

    return (
        tA_window,
        tB_window,
        di_window,
        approach_frames_A,
        approach_frames_B,
        approach_durations_A,
        approach_durations_B,
    )


def get_novel_side(file_index):
    """Return the magnet side for the zero-based positional metadata index."""
    return side_list[file_index]


def draw_skeleton(ax, fcoords, node_names, color, alpha=0.3):
    """Draw the configured skeleton edges when both node names are present."""
    for n1, n2 in skeleton_edges:
        if n1 in node_names and n2 in node_names:
            i1, i2 = node_names.index(n1), node_names.index(n2)
            ax.plot(
                [fcoords[i1, 0], fcoords[i2, 0]],
                [fcoords[i1, 1], fcoords[i2, 1]],
                color=color,
                alpha=alpha,
                linewidth=1,
            )


def plot_skeleton_summary_windows(
    coords,
    node_names,
    objA,
    objB,
    novel_side,
    radius,
    di,
    approachA,
    approachB,
    title,
    window,
    savepath=None,
):
    """Overlay trajectories, approach frames, and object regions.

    For magnet-side A/B recordings, red/pink marks the magnetic object and
    blue marks the nonmagnetic object. "NO" retains the original neutral
    display. The displayed times use the global acquisition rate, fps.
    """
    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot all skeletons in light gray
    for frame in range(coords.shape[0]):
        draw_skeleton(ax, coords[frame], node_names, "lightgray", alpha=0.1)

    # Objects + ROIs
    if novel_side == "NO":
        ax.scatter(*objA, c="red", s=100, marker="o", label="object_A")
        ax.scatter(*objB, c="blue", s=100, marker="o", label="object_B")
    elif novel_side == "B":
        # Highlight approach frames
        for frame in approachA:
            draw_skeleton(ax, coords[frame], node_names, "skyblue", alpha=0.3)
        for frame in approachB:
            draw_skeleton(ax, coords[frame], node_names, "pink", alpha=0.3)
        ax.scatter(*objA, c="blue", s=100, marker="o", label="wo-Magnet")
        ax.scatter(*objB, c="red", s=100, marker="^", label="w-Magnets")
        ax.add_patch(
            plt.Circle(objA, radius, color="blue", fill=False, linestyle="--", linewidth=2)
        )
        ax.add_patch(plt.Circle(objB, radius, color="red", fill=False, linestyle="--", linewidth=2))
    elif novel_side == "A":
        for frame in approachA:
            draw_skeleton(ax, coords[frame], node_names, "pink", alpha=0.3)
        for frame in approachB:
            draw_skeleton(ax, coords[frame], node_names, "skyblue", alpha=0.3)
        ax.scatter(*objA, c="red", s=100, marker="^", label="w-Magnets")
        ax.scatter(*objB, c="blue", s=100, marker="o", label="wo-Magnet")
        ax.add_patch(plt.Circle(objA, radius, color="red", fill=False, linestyle="--", linewidth=2))
        ax.add_patch(
            plt.Circle(objB, radius, color="blue", fill=False, linestyle="--", linewidth=2)
        )

    # Fixed arena view; invert y to match the image coordinate system.
    ax.set_xlim(140, 1140)
    ax.set_ylim(12, 1012)
    ax.invert_yaxis()
    ax.set_xlabel("X (px)")
    ax.set_ylabel("Y (px)")

    # Calculate approach times for display
    tA = len(approachA) / fps
    tB = len(approachB) / fps

    # Text box
    text = f"DI={di:.2f}\nLeft={tA:.1f}s\nRight={tB:.1f}s\nfirst-{window}s"
    ax.text(
        0.02,
        0.02,
        text,
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
    )

    ax.set_title(title)
    ax.legend()

    if savepath:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")
        plt.close()  # Close figure to save memory
    else:
        plt.show()


def save_approach_frames(filename, approach_frames_A, approach_frames_B, result_folder, window):
    """Save zero-based approach indices and counts for the current window."""
    approach_data = {
        "file": filename,
        "approach_frames_A": approach_frames_A,
        "approach_frames_B": approach_frames_B,
        "total_frames_A": len(approach_frames_A),
        "total_frames_B": len(approach_frames_B),
    }

    json_path = result_folder / f"{Path(filename).stem}_approach_frames_{window}s.json"
    with open(json_path, "w") as f:
        json.dump(approach_data, f, indent=2)


# ----------------------------
# 5. Batch analysis and export
# ----------------------------
def main():
    """Analyze configured recordings and export one summary per time window."""
    result_folder.mkdir(exist_ok=True)
    print(f"Starting analysis for {subset_name} dataset")
    print(f"Processing {len(object_positions)} files")
    print(f"Fixed window analysis: {windows}")
    print("-" * 50)

    print("\n=== FIXED WINDOW ANALYSIS ===")

    for window in windows:
        print(f"\nAnalyzing {window}s window...")
        results_window = []
        # Use configuration order to retain each recording's side and group.
        # Mouse IDs need not increase in the original recording order.
        for file_index, filename in enumerate(object_positions):
            file = h5_folder / filename
            if not file.is_file():
                print(f"Skipping missing file: {filename}")
                continue

            try:
                coords, node_names = extract_all_nodes(file)
                objA, objB = object_positions[file.name]["A"], object_positions[file.name]["B"]

                # The metadata index stays fixed even if an earlier file failed.
                novel_side = get_novel_side(file_index)
                group_id = group_list[file_index]
                # Fixed window analysis
                (
                    tA_window,
                    tB_window,
                    di_window,
                    approach_frames_A,
                    approach_frames_B,
                    approach_durations_A,
                    approach_durations_B,
                ) = exploration_time_fixed_window(
                    coords, node_names, objA, objB, radius_px, fps, window
                )

                print(
                    f"{file.name}: {window}s window - A={tA_window:.1f}s, B={tB_window:.1f}s, DI={di_window:.2f}"
                )

                # Save frame indices for each configured window.
                save_approach_frames(
                    file.name, approach_frames_A, approach_frames_B, result_folder, window
                )

                # Plot and save
                savepath = result_folder / f"{file.stem}_{radius_px}px_first{window}sec.pdf"
                # Retain the original extra background frame in the plot;
                # exploration statistics still use exactly int(window * fps).
                plot_skeleton_summary_windows(
                    coords[0 : window * fps + 1],
                    node_names,
                    objA,
                    objB,
                    novel_side,
                    radius_px,
                    di_window,
                    approach_frames_A,
                    approach_frames_B,
                    file.stem,
                    window,
                    savepath,
                )

                # Calculate statistics
                mean_duration_A = np.mean(approach_durations_A) if approach_durations_A else 0
                mean_duration_B = np.mean(approach_durations_B) if approach_durations_B else 0
                num_episodes_A = len(approach_durations_A)
                num_episodes_B = len(approach_durations_B)

                # Preserve downstream column names and orient DI toward the magnet.

                results_window.append(
                    {
                        "file": file.name,
                        "group": group_id,
                        f"tA_window{window}s": round(tA_window, 2),
                        f"tB_window{window}s": round(tB_window, 2),
                        f"DI_window_B_A{window}s": round(di_window, 3),
                        "mean_approach_durations_A": round(mean_duration_A, 2),
                        "mean_approach_durations_B": round(mean_duration_B, 2),
                        "num_episodes_A": num_episodes_A,
                        "num_episodes_B": num_episodes_B,
                        "novel_side": novel_side,
                        "DI_Magnet": round(-di_window if novel_side == "A" else di_window, 3),
                        "t_Magnet": round(tA_window if novel_side == "A" else tB_window, 2),
                        "t_nonMag": round(tB_window if novel_side == "A" else tA_window, 2),
                        "dur_Magnet": round(
                            mean_duration_A if novel_side == "A" else mean_duration_B, 2
                        ),
                        "dur_nonMag": round(
                            mean_duration_B if novel_side == "A" else mean_duration_A, 2
                        ),
                        "epi_Magnet": round(
                            num_episodes_A if novel_side == "A" else num_episodes_B, 2
                        ),
                        "epi_nonMag": round(
                            num_episodes_B if novel_side == "A" else num_episodes_A, 2
                        ),
                    }
                )

            except Exception as e:
                print(f"Error processing {file}: {e}")

        # Sort once after collecting rows; equal-group ordering remains stable.
        if results_window:
            results_window.sort(key=lambda x: x["group"], reverse=True)
            results_df = pd.DataFrame(results_window)
            date_str = datetime.now().strftime("%Y%m%d")
            out_csv = result_folder / f"summary_{subset_name}_{window}s_window_{date_str}.csv"
            results_df.to_csv(out_csv, index=False)
            print(f"Results saved to {out_csv}")

    print(f"\n Analysis complete! All results saved to {result_folder}")


if __name__ == "__main__":
    main()
