# -*- coding: utf-8 -*-
"""
Bilateral Concordance Analysis
Compares left and right hemisphere results to identify robust bilateral changes
@author: zhang

"""


import numpy as np
import matplotlib.pyplot as plt
from skimage import io
import csv
import pandas as pd

# =============================================================================
# Configuration - Update these paths
# =============================================================================

WORKING_PATH    = r"\resource_folder_1"
ATLAS_PATH      = WORKING_PATH + r"\Wholebrain_New_Annotation_Map_T.tif"  
PVALUE_MAP_PATH = WORKING_PATH + r"\pvaluemap_L_T.tif"
OUTPUT_PATH     = WORKING_PATH + r"\volcano_plot_L.pdf"
CSV_PATH        = WORKING_PATH + r"\volcano_plot_L_summary_stats.csv"
XLSX_PATH       = WORKING_PATH + r"\volcano_plot_L_summary_stats.xlsx"

# Hemisphere selection: 'left', 'right', or 'both'
HEMISPHERE = 'left'  # Set to 'left' or 'right' for single hemisphere analysis
# Note: For single hemisphere analysis, use corresponding p-value map (L or R)
# Atlas will be automatically split at the midline

# Labeling control: Set how many regions to label on each side
MAX_LABELS_DECREASED = 45  # Maximum labels for decreased regions (negative x-axis)
MAX_LABELS_INCREASED = 30  # Maximum labels for increased regions (positive x-axis)

# Additional labeling control: Label top N regions by absolute voxel count
MAX_LABELS_BY_COUNT_DECREASED = 10  # Top N decreased regions by voxel count
MAX_LABELS_BY_COUNT_INCREASED = 10  # Top N increased regions by voxel count

# Read the Excel file
df = pd.read_excel(WORKING_PATH + r'\New_Annotation_Table.xlsx')
# Create dictionary mapping id (first column) to acronym (third column)
REGION_NAMES = dict(zip(df['id'], df['acronym']))
# Create dictionary mapping id to full name
REGION_FULL_NAMES = dict(zip(df['id'], df['name']))

# P-value map format: 2-channel 32bit TIFF stack
# Channel 0: decreased values
# Channel 1: increased values
# Non-zero values indicate significance (p < 0.05)

# =============================================================================
# Functions
# =============================================================================

def load_tif_stack(filepath):
    """Load a multi-page TIFF stack."""
    return io.imread(filepath)

def split_hemisphere(volume, hemisphere='both'):
    """
    Split volume to get only left or right hemisphere.
    Assumes the atlas is symmetric with midline at the center of the width dimension.
    
    Parameters:
        volume: 3D or 4D numpy array (slices, [channels,] height, width)
        hemisphere: 'left', 'right', or 'both'
    
    Returns:
        volume (modified or original)
    """
    if hemisphere == 'both':
        return volume
    
    # Determine the shape and find the width dimension (last dimension)
    if len(volume.shape) == 3:
        # Shape: (slices, height, width)
        slices, height, width = volume.shape
        midline = width // 2
        
        if hemisphere == 'left':
            # Keep left half (assuming left is first half)
            return volume[:, :, :midline]
        elif hemisphere == 'right':
            # Keep right half
            return volume[:, :, midline:]
        else:
            raise ValueError(f"Invalid hemisphere option: {hemisphere}. Use 'left', 'right', or 'both'")
    
    elif len(volume.shape) == 4:
        # Shape: (slices, channels, height, width)
        slices, channels, height, width = volume.shape
        midline = width // 2
        
        if hemisphere == 'left':
            return volume[:, :, :, :midline]
        elif hemisphere == 'right':
            return volume[:, :, :, midline:]
        else:
            raise ValueError(f"Invalid hemisphere option: {hemisphere}. Use 'left', 'right', or 'both'")
    
    else:
        raise ValueError(f"Unexpected volume shape: {volume.shape}")

def identify_pvalue_voxels(pvalue_map):
    """
    Identify increased and decreased voxels from the p-value map.
    P-value map is a 2-channel 32bit TIFF stack.
    Shape: (slices, 2, height, width) e.g., (254, 2, 577, 375)
    Channel 0: decreased values
    Channel 1: increased values
    Non-zero values indicate significance.
    
    Returns two binary masks.
    """
    # pvalue_map shape: (slices, 2, height, width)
    if len(pvalue_map.shape) == 4 and pvalue_map.shape[1] == 2:
        decreased_mask = pvalue_map[:, 0, :, :] != 0
        increased_mask = pvalue_map[:, 1, :, :] != 0
    else:
        raise ValueError(f"Unexpected p-value map shape: {pvalue_map.shape}. "
                        f"Expected shape: (slices, 2, height, width)")
    
    return increased_mask, decreased_mask

def calculate_region_statistics(atlas, increased_mask, decreased_mask, region_names, region_full_names):
    """
    Calculate statistics for each brain region.
    
    Returns:
        results: dict with region names as keys and statistics as values
    """
    results = {}
    
    # Get total brain size (non-zero voxels in atlas)
    total_brain_voxels = np.sum(atlas > 0)
    
    # Get unique region labels
    unique_labels = np.unique(atlas)
    unique_labels = unique_labels[unique_labels > 0]  # Exclude background
    
    for label in unique_labels:
        if label not in region_names:
            continue
            
        region_name = region_names[label]
        region_full_name = region_full_names.get(label, "")  # Get full name, default to empty string
        region_mask = (atlas == label)
        
        # Count voxels in this region
        region_size = np.sum(region_mask)
        
        # Count significant voxels in this region
        increased_in_region = np.sum(region_mask & increased_mask)
        decreased_in_region = np.sum(region_mask & decreased_mask)
        
        # Calculate percentages relative to THIS REGION'S size (MODIFIED)
        increased_pct = (increased_in_region / region_size) * 100 if region_size > 0 else 0
        decreased_pct = (decreased_in_region / region_size) * 100 if region_size > 0 else 0
        
        # Total significant voxels
        total_significant = increased_in_region + decreased_in_region
        
        results[region_name] = {
            'full_name': region_full_name,
            'region_size': region_size,
            'increased_count': increased_in_region,
            'decreased_count': decreased_in_region,
            'increased_pct': increased_pct,
            'decreased_pct': decreased_pct,
            'total_significant': total_significant
        }
    
    return results, total_brain_voxels

def create_volcano_plot(results, output_path, max_labels_decreased=15, max_labels_increased=15,
                       max_labels_by_count_decreased=5, max_labels_by_count_increased=5):
    """
    Create a volcano-style plot similar to the example.
    X-axis: percentage for decreased (negative) and increased (positive) - RELATIVE TO REGION SIZE
    Y-axis: significant voxel counts (log scale)
    Regions with both increased and decreased are plotted twice.
    Uses adjustText library to avoid label overlaps.
    
    Parameters:
        results: dictionary of region statistics
        output_path: path to save the plot
        max_labels_decreased: maximum number of labels for decreased regions (negative x-axis)
        max_labels_increased: maximum number of labels for increased regions (positive x-axis)
        max_labels_by_count_decreased: top N decreased regions by absolute voxel count to label
        max_labels_by_count_increased: top N increased regions by absolute voxel count to label
    """
    # Set font to Arial
    plt.rcParams['font.family'] = 'Arial'
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Prepare data for plotting
    regions = []
    x_values = []
    y_values = []
    colors_list = []
    
    for region_name, stats in results.items():
        # Plot decreased voxels (negative x-axis)
        if stats['decreased_count'] > 0:
            regions.append(region_name)
            x_values.append(-stats['decreased_pct'])  # Negative for decreased
            y_values.append(stats['decreased_count'])
            colors_list.append('#E91E63')  # Pink for decreased
        
        # Plot increased voxels (positive x-axis)
        if stats['increased_count'] > 0:
            regions.append(region_name)
            x_values.append(stats['increased_pct'])  # Positive for increased
            y_values.append(stats['increased_count'])
            colors_list.append('#00BCD4')  # Cyan for increased
    
    # Convert to numpy arrays
    x_values = np.array(x_values)
    y_values = np.array(y_values)
    
    # Plot all points in gray first
    ax.scatter(x_values, y_values, c='lightgray', s=30, alpha=0.3, zorder=1)
    
    # Plot colored points on top
    ax.scatter(x_values, y_values, c=colors_list, s=30, alpha=0.6, zorder=2)
    
    # Set log scale for y-axis and symlog for x-axis
    ax.set_yscale('log')
    ax.set_xscale('symlog', linthresh=0.1)  # linthresh: linear range near zero
    
    # define x-axis ticks and labels
    xtick_values = [-100, -50, -20, -10, -5, -1, 0, 1, 5, 10, 20, 50, 100]
    xtick_labels = ['100', '50', '20', '10', '5', '1', '0', '1', '5', '10', '20', '50', '100']
    ax.set_xticks(xtick_values)
    ax.set_xticklabels(xtick_labels)

    # Select which regions to label based on importance - SEPARATELY for each side
    # Separate indices for decreased (negative x) and increased (positive x)
    decreased_indices = np.where(x_values < 0)[0]
    increased_indices = np.where(x_values > 0)[0]
    
    indices_to_label = []
    
    # Process decreased regions (negative x-axis)
    if len(decreased_indices) > 0:
        decreased_y = y_values[decreased_indices]
        decreased_x_abs = np.abs(x_values[decreased_indices])
        
        # Get indices sorted by count (y-axis)
        indices_by_count_local = np.argsort(decreased_y)[::-1]
        top_by_count = decreased_indices[indices_by_count_local[:max_labels_by_count_decreased]]
        
        # Get indices sorted by percentage (x-axis)
        indices_by_pct_local = np.argsort(decreased_x_abs)[::-1]
        top_by_pct = decreased_indices[indices_by_pct_local[:max_labels_decreased]]
        
        # Combine and get unique indices
        selected_decreased = np.unique(np.concatenate([top_by_count, top_by_pct]))
        
        print(f"Labeling {len(selected_decreased)} decreased regions "
              f"(top {max_labels_by_count_decreased} by count + top {max_labels_decreased} by %, "
              f"out of {len(decreased_indices)} total)")
        
        indices_to_label.extend(selected_decreased)
    
    # Process increased regions (positive x-axis)
    if len(increased_indices) > 0:
        increased_y = y_values[increased_indices]
        increased_x_abs = np.abs(x_values[increased_indices])
        
        # Get indices sorted by count (y-axis)
        indices_by_count_local = np.argsort(increased_y)[::-1]
        top_by_count = increased_indices[indices_by_count_local[:max_labels_by_count_increased]]
        
        # Get indices sorted by percentage (x-axis)
        indices_by_pct_local = np.argsort(increased_x_abs)[::-1]
        top_by_pct = increased_indices[indices_by_pct_local[:max_labels_increased]]
        
        # Combine and get unique indices
        selected_increased = np.unique(np.concatenate([top_by_count, top_by_pct]))
        
        print(f"Labeling {len(selected_increased)} increased regions "
              f"(top {max_labels_by_count_increased} by count + top {max_labels_increased} by %, "
              f"out of {len(increased_indices)} total)")
        
        indices_to_label.extend(selected_increased)
    
    indices_to_label = np.array(indices_to_label)

    # Try to use adjustText if available, otherwise fall back to manual spacing
    try:
        from adjustText import adjust_text
        
        # Create text annotations ONLY for selected regions
        texts = []
        for idx in indices_to_label:
            txt = ax.text(x_values[idx], y_values[idx], regions[idx],
                         fontsize=8,
                         family='Arial'
                         # bbox=dict(boxstyle='round,pad=0.2', facecolor='white', 
                         #          edgecolor='none', alpha=0.7)
                         )
            texts.append(txt)
        
        # Adjust text positions to avoid overlaps
        adjust_text(texts, 
                   x=x_values[indices_to_label], 
                   y=y_values[indices_to_label],
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.6),
                   expand_points=(1.5, 1.5),
                   force_points=(0.5, 0.5),
                   force_text=(3, 3),
                   ax=ax)
        
        print(f"Successfully labeled regions using adjustText")
    
    except ImportError:
        print("Warning: adjustText not available. Using manual label spacing.")
        print("Install with: pip install adjustText")
        
        # Fallback: manual annotation
        for idx in indices_to_label:
            # Vary offset based on position to reduce overlap
            offset_x = 25 if x_values[idx] > 0 else -25
            offset_y = 15 + (idx % 3) * 10  # Stagger vertically
            
            ax.annotate(regions[idx], 
                       xy=(x_values[idx], y_values[idx]),
                       xytext=(offset_x, offset_y),
                       textcoords='offset points',
                       fontsize=7,
                       family='Arial',
                       ha='left' if x_values[idx] > 0 else 'right',
                       va='bottom',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                edgecolor='lightgray', alpha=0.8),
                       arrowprops=dict(arrowstyle='->', 
                                     connectionstyle='arc3,rad=0.2',
                                     color='gray', 
                                     linewidth=0.8,
                                     alpha=0.7))
    
    # Adjust x-axis range to fit data with some padding
    if len(x_values) > 0:
        x_abs_max = max(abs(x_values.min()), abs(x_values.max()))
        ax.set_xlim(-x_abs_max * 3, x_abs_max * 3)
    
    ax.set_ylim(0.5, max(y_values) * 2)
    
    # Labels and title
    ax.set_xlabel('% of region voxels that are significant\n' + 
                  'Decreased (in CMF)' + ' ' * 20 + 'Increased (in CMF)',
                  fontsize=13, fontweight='normal', family='Arial')
    ax.set_ylabel('significant voxel counts', fontsize=13, fontweight='normal', family='Arial')
    ax.set_title('Brain Region P-value Volcano Plot', fontsize=15, fontweight='bold', family='Arial', pad=15)
    
    # Set tick label font sizes
    ax.tick_params(axis='both', which='major', labelsize=11)
    
    # Add vertical line at x=0
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    
    # Grid - enable for both major ticks
    ax.grid(True, alpha=0.3, which='major', linestyle='-', linewidth=0.5)
    
    # Create custom legend - move to lower left
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E91E63', label='Decreased (in CMF)'),
        Patch(facecolor='#00BCD4', label='Increased (in CMF)')
    ]
    ax.legend(handles=legend_elements, loc='lower left', fontsize=11,
             framealpha=0.9, edgecolor='gray', prop={'family': 'Arial'})
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_path}")
    
    return fig, ax

def print_summary(results, total_brain_voxels, csv_path="example_summary_stats.csv", xlsx_path="example_summary_stats.xlsx"):
    """Print summary statistics and save table as CSV/XLSX."""
    
    print("\n" + "="*90)
    print("SUMMARY STATISTICS")
    print("="*90)
    print(f"Total brain voxels: {total_brain_voxels:,}")
    print(f"\nNumber of regions analyzed: {len(results)}")
    print("\nAll regions by total significant voxels:")
    print("-"*90)

    # Sort by total significant voxels
    sorted_regions = sorted(results.items(),
                            key=lambda x: x[1]['total_significant'],
                            reverse=True)

    # Print header
    print(f"{'Region':<10} {'Region Size':<12} {'Total Sig':<12} {'Increased':<12} "
          f"{'Decreased':<12} {'Inc %':<10} {'Dec %':<10}")
    print("-"*90)

    # Prepare table rows for saving
    rows = []

    for region, stats in sorted_regions:
        print(f"{region:<10} {stats['region_size']:<12} {stats['total_significant']:<12} "
              f"{stats['increased_count']:<12} {stats['decreased_count']:<12} "
              f"{stats['increased_pct']:>8.2f}%  {stats['decreased_pct']:>8.2f}%")

        rows.append({
            "Region": region,
            "Full Name": stats["full_name"],
            "Region Size": stats["region_size"],
            "Total Significant Voxels": stats["total_significant"],
            "Increased Count": stats["increased_count"],
            "Decreased Count": stats["decreased_count"],
            "Increased % (of region)": stats["increased_pct"],
            "Decreased % (of region)": stats["decreased_pct"],
        })

    # --- Save CSV ---
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # --- Save XLSX ---
    df = pd.DataFrame(rows)
    df.to_excel(xlsx_path, index=False)

    print(f"\nSaved CSV to: {csv_path}")
    print(f"Saved XLSX to: {xlsx_path}")

# =============================================================================
# Main Analysis
# =============================================================================

def main():
    print("Loading data...")
    print(f"Hemisphere setting: {HEMISPHERE}")
    
    # Load atlas and p-value map
    atlas = load_tif_stack(ATLAS_PATH)
    pvalue_map = load_tif_stack(PVALUE_MAP_PATH)
    
    print(f"Original atlas shape: {atlas.shape}")
    print(f"Original p-value map shape: {pvalue_map.shape}")
    
    # Split hemispheres if needed
    atlas = split_hemisphere(atlas, HEMISPHERE)
    pvalue_map = split_hemisphere(pvalue_map, HEMISPHERE)
    
    print(f"Processing atlas shape: {atlas.shape}")
    print(f"Processing p-value map shape: {pvalue_map.shape}")
    
    # Identify increased and decreased voxels
    print("\nIdentifying significant voxels...")
    increased_mask, decreased_mask = identify_pvalue_voxels(pvalue_map)
    
    print(f"Total increased voxels: {np.sum(increased_mask):,}")
    print(f"Total decreased voxels: {np.sum(decreased_mask):,}")
    
    # Calculate statistics for each region
    print("\nCalculating region statistics...")
    results, total_brain_voxels = calculate_region_statistics(
        atlas, increased_mask, decreased_mask, REGION_NAMES, REGION_FULL_NAMES
    )
    
    # Print summary
    print_summary(results, total_brain_voxels, csv_path=CSV_PATH, xlsx_path=XLSX_PATH)
    
    # Create volcano plot
    print("\nCreating volcano plot...")
    create_volcano_plot(results, OUTPUT_PATH, 
                       max_labels_decreased=MAX_LABELS_DECREASED,
                       max_labels_increased=MAX_LABELS_INCREASED,
                       max_labels_by_count_decreased=MAX_LABELS_BY_COUNT_DECREASED,
                       max_labels_by_count_increased=MAX_LABELS_BY_COUNT_INCREASED)
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()