# -*- coding: utf-8 -*-
"""
Bilateral Concordance Analysis
Compares left and right hemisphere results to identify robust bilateral changes
@author: zhang

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# =============================================================================
# Configuration
# =============================================================================

WORKING_PATH = r"\resource_folder_2"
LEFT_FILE = WORKING_PATH + r"\volcano_plot_L_summary_stats.xlsx"
RIGHT_FILE = WORKING_PATH + r"\volcano_plot_R_summary_stats.xlsx"
OUTPUT_PLOT = WORKING_PATH + r"\bilateral_concordance_plot_100.pdf"
OUTPUT_TABLE = WORKING_PATH + r"\bilateral_concordance_summary.xlsx"

# Concordance threshold - only plot regions present in both hemispheres
MIN_SIGNIFICANCE_THRESHOLD = 0.1  # Minimum % to consider (filter out very small changes)

# Labeling control
MAX_LABELS_DECREASED = 100  # Maximum number of labels for decreased regions (negative x-axis)
MAX_LABELS_INCREASED = 100  # Maximum number of labels for increased regions (positive x-axis)

# =============================================================================
# Functions
# =============================================================================

def load_hemisphere_data(filepath):
    """Load hemisphere data from Excel file."""
    df = pd.read_excel(filepath)
    return df

def calculate_bilateral_metrics(left_df, right_df):
    """
    Calculate concordance metrics for regions present in both hemispheres.
    
    Returns:
        bilateral_df: DataFrame with concordance metrics
    """
    results = []
    
    # Merge on Region name
    merged = pd.merge(left_df, right_df, on='Region', suffixes=('_L', '_R'), how='inner')
    
    print(f"Found {len(merged)} regions present in both hemispheres")
    
    for idx, row in merged.iterrows():
        region = row['Region']
        full_name = row.get('Full Name_L', row.get('Full Name_R', ''))
        
        # Get percentages for each hemisphere
        left_inc_pct = row['Increased % (of region)_L']
        left_dec_pct = row['Decreased % (of region)_L']
        right_inc_pct = row['Increased % (of region)_R']
        right_dec_pct = row['Decreased % (of region)_R']
        
        # Get counts
        left_inc_count = row['Increased Count_L']
        left_dec_count = row['Decreased Count_L']
        right_inc_count = row['Increased Count_R']
        right_dec_count = row['Decreased Count_R']
        
        # Determine dominant direction in each hemisphere
        left_net_pct = left_inc_pct - left_dec_pct
        right_net_pct = right_inc_pct - right_dec_pct
        
        # Calculate average effect size (net percentage)
        avg_net_pct = (left_net_pct + right_net_pct) / 2
        
        # Calculate concordance score (correlation-like measure)
        # If both hemispheres change in same direction, concordance is high
        # Concordance = 1 - normalized_difference
        max_abs = max(abs(left_net_pct), abs(right_net_pct))
        if max_abs > 0:
            concordance = 1 - abs(left_net_pct - right_net_pct) / (2 * max_abs)
        else:
            concordance = 0
        
        # Alternative: Pearson-like concordance (both positive or both negative)
        if left_net_pct * right_net_pct > 0:  # Same sign
            concordance_direction = 1
        else:  # Opposite signs
            concordance_direction = -1
        
        # Total significant voxels (combined both hemispheres)
        total_significant = (left_inc_count + left_dec_count + 
                           right_inc_count + right_dec_count)
        
        # Determine overall direction
        if avg_net_pct > 0:
            direction = 'Increased'
            color = '#00BCD4'  # Cyan
        else:
            direction = 'Decreased'
            color = '#E91E63'  # Pink
        
        results.append({
            'Region': region,
            'Full Name': full_name,
            'Left_Net_%': left_net_pct,
            'Right_Net_%': right_net_pct,
            'Average_Net_%': avg_net_pct,
            'Concordance': concordance,
            'Concordance_Direction': concordance_direction,
            'Total_Significant_Voxels': total_significant,
            'Direction': direction,
            'Color': color,
            'Left_Inc_%': left_inc_pct,
            'Left_Dec_%': left_dec_pct,
            'Right_Inc_%': right_inc_pct,
            'Right_Dec_%': right_dec_pct,
        })
    
    bilateral_df = pd.DataFrame(results)
    return bilateral_df

def create_bilateral_concordance_plot(bilateral_df, output_path, max_labels_decreased=15, max_labels_increased=15):
    """
    Create a concordance plot showing bilateral consistency.
    
    X-axis: Average net % change (across both hemispheres)
    Y-axis: Concordance score (how similar L and R are)
    
    Parameters:
        bilateral_df: DataFrame with bilateral metrics
        output_path: path to save the plot
        max_labels_decreased: maximum number of labels for decreased regions (negative x-axis)
        max_labels_increased: maximum number of labels for increased regions (positive x-axis)
    """
    # Set font to Arial
    plt.rcParams['font.family'] = 'Arial'
 
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Filter out very small changes
    filtered_df = bilateral_df[
        (np.abs(bilateral_df['Average_Net_%']) > MIN_SIGNIFICANCE_THRESHOLD)
    ].copy()
    
    print(f"Plotting {len(filtered_df)} regions (filtered by threshold)")
    
    if len(filtered_df) == 0:
        print("Warning: No regions pass the significance threshold!")
        return None, None
    
    # Prepare data
    x_values = filtered_df['Average_Net_%'].values
    y_values = filtered_df['Concordance'].values
    colors = filtered_df['Color'].values
    sizes = filtered_df['Total_Significant_Voxels'].values
    regions = filtered_df['Region'].values
    
    # Normalize sizes for plotting (scale between 30 and 500)
    if len(sizes) > 0 and sizes.max() > 0:
        sizes_normalized = 30 + (sizes - sizes.min()) / (sizes.max() - sizes.min() + 1e-10) * 470
    else:
        sizes_normalized = np.ones_like(sizes) * 100
    
    # Plot all points with gray background
    ax.scatter(x_values, y_values, s=sizes_normalized, 
               c='lightgray', alpha=0.3, zorder=1, edgecolors='none')
    
    # Plot colored points on top
    ax.scatter(x_values, y_values, s=sizes_normalized, 
               c=colors, alpha=0.6, zorder=2, edgecolors='white', linewidths=0.5)
    
    # Set axis scales
    ax.set_ylim(-0.1, 1.1)
    
    # Get x-axis range for symmetric display
    # x_abs_max = max(abs(x_values.min()), abs(x_values.max()))
    # ax.set_xlim(-x_abs_max * 1.1, x_abs_max * 1.1)
    
    # Set x-axis ticks to show positive values on both sides
    x_tick_positions = ax.get_xticks()
    x_tick_labels = [f'{abs(x):.0f}' for x in x_tick_positions]
    ax.set_xticklabels(x_tick_labels)
    
    # Add reference lines
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.3, linewidth=1)
    
    # Select regions to label - SEPARATELY for decreased and increased
    indices_to_label = []
    
    # Separate into decreased (negative x) and increased (positive x)
    decreased_mask = filtered_df['Average_Net_%'] < 0
    increased_mask = filtered_df['Average_Net_%'] > 0
    
    decreased_df = filtered_df[decreased_mask].copy()
    increased_df = filtered_df[increased_mask].copy()
    
    # Process decreased regions (negative x-axis)
    if len(decreased_df) > 0:
        # Score based on concordance * abs(effect size)
        decreased_df['label_score'] = decreased_df['Concordance'] * np.abs(decreased_df['Average_Net_%'])
        
        if len(decreased_df) > max_labels_decreased:
            top_decreased = decreased_df.nlargest(max_labels_decreased, 'label_score')
            print(f"Labeling top {len(top_decreased)} decreased regions (out of {len(decreased_df)} total)")
        else:
            top_decreased = decreased_df
            print(f"Labeling all {len(decreased_df)} decreased regions")
        
        # Get indices in the original filtered_df
        selected_decreased_indices = filtered_df.index.get_indexer(top_decreased.index)
        indices_to_label.extend(selected_decreased_indices)
    
    # Process increased regions (positive x-axis)
    if len(increased_df) > 0:
        # Score based on concordance * abs(effect size)
        increased_df['label_score'] = increased_df['Concordance'] * np.abs(increased_df['Average_Net_%'])
        
        if len(increased_df) > max_labels_increased:
            top_increased = increased_df.nlargest(max_labels_increased, 'label_score')
            print(f"Labeling top {len(top_increased)} increased regions (out of {len(increased_df)} total)")
        else:
            top_increased = increased_df
            print(f"Labeling all {len(increased_df)} increased regions")
        
        # Get indices in the original filtered_df
        selected_increased_indices = filtered_df.index.get_indexer(top_increased.index)
        indices_to_label.extend(selected_increased_indices)
    
    indices_to_label = np.array(indices_to_label)
    
    # Try to use adjustText for better label placement
    try:
        from adjustText import adjust_text
        
        texts = []
        for idx in indices_to_label:
            txt = ax.text(x_values[idx], y_values[idx], regions[idx],
                         fontsize=3)
            texts.append(txt)
        
        adjust_text(texts,
                   x=x_values[indices_to_label],
                   y=y_values[indices_to_label],
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.6),
                   expand_points=(1.5, 1.5),
                   force_points=(0.5, 0.5),
                   force_text=(0.5, 0.5),
                   ax=ax)
        
        print("Successfully labeled regions using adjustText")
        
    except ImportError:
        print("Warning: adjustText not available. Using manual labels.")
        
        for idx in indices_to_label:
            offset_x = 2 if x_values[idx] > 0 else -2
            offset_y = 0.03
            
            ax.annotate(regions[idx],
                       xy=(x_values[idx], y_values[idx]),
                       xytext=(offset_x, offset_y),
                       textcoords='offset points',
                       fontsize=8,
                       ha='left' if x_values[idx] > 0 else 'right',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                edgecolor='lightgray', alpha=0.8),
                       arrowprops=dict(arrowstyle='->', color='gray',
                                     linewidth=0.8, alpha=0.7))
    
    # Labels and formatting
    ax.set_xlabel('Average Net % Change (Left + Right Hemispheres)\n' +
                  'Decreased (in CMF) ← → Increased (in CMF)',
                  fontsize=14, fontweight='bold')
    ax.set_ylabel('Bilateral Concordance Score\n(1.0 = Perfect Agreement between L & R)',
                  fontsize=14, fontweight='bold')
    ax.set_title('Bilateral Concordance Analysis: Regions Changed in Both Hemispheres',
                 fontsize=16, fontweight='bold', pad=20)
    
    # Add grid
    ax.grid(True, alpha=0.2, which='major', linestyle='-', linewidth=0.5)
    
    # Legend - both in lower left
    legend_elements = [
        Patch(facecolor='#E91E63', label='Decreased (in CMF)'),
        Patch(facecolor='#00BCD4', label='Increased (in CMF)'),
    ]
    
    # Add size legend with more diverse sizes
    size_legend_elements = [
        plt.scatter([], [], s=30, c='gray', alpha=0.6, label='Low'),
        plt.scatter([], [], s=150, c='gray', alpha=0.6, label='Medium'),
        plt.scatter([], [], s=300, c='gray', alpha=0.6, label='High'),
    ]
    
    # Place first legend (direction) at lower left
    legend1 = ax.legend(handles=legend_elements, loc='lower left', 
                       bbox_to_anchor=(0, 0),
                       title='Direction', fontsize=10, title_fontsize=10,
                       framealpha=0.9, edgecolor='gray', prop={'family': 'Arial'})
    
    # Place second legend (size) slightly to the right of first legend
    legend2 = ax.legend(handles=size_legend_elements, loc='lower left',
                       bbox_to_anchor=(0.25, 0),
                       title='Total voxel count', fontsize=10, title_fontsize=10,
                       framealpha=0.9, edgecolor='gray', prop={'family': 'Arial'})
    
    # Add the first legend back (since creating legend2 removes legend1)
    ax.add_artist(legend1)
    
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to {output_path}")
    
    return fig, ax

def save_bilateral_summary(bilateral_df, output_path):
    """Save bilateral concordance summary to Excel."""
    
    # Sort by concordance score (descending)
    sorted_df = bilateral_df.sort_values('Concordance', ascending=False)
    
    # Add effect size column (absolute value of average net %)
    sorted_df['Effect_Size'] = np.abs(sorted_df['Average_Net_%'])
    
    # Add concordance × effect size score
    sorted_df['Concordance_x_Effect'] = sorted_df['Concordance'] * sorted_df['Effect_Size']
    
    # Select relevant columns for output
    output_df = sorted_df[[
        'Region',
        'Full Name',
        'Direction',
        'Average_Net_%',
        'Effect_Size',
        'Concordance',
        'Concordance_x_Effect',
        'Left_Net_%',
        'Right_Net_%',
        'Total_Significant_Voxels',
        'Left_Inc_%',
        'Left_Dec_%',
        'Right_Inc_%',
        'Right_Dec_%',
    ]].copy()
    
    # Rename columns for clarity
    output_df.columns = [
        'Region',
        'Full Name',
        'Direction',
        'Average Net % Change',
        'Effect Size (Absolute %)',
        'Bilateral Concordance',
        'Concordance × Effect Size',
        'Left Net %',
        'Right Net %',
        'Total Significant Voxels',
        'Left Increased %',
        'Left Decreased %',
        'Right Increased %',
        'Right Decreased %',
    ]
    
    # Save to Excel
    output_df.to_excel(output_path, index=False)
    print(f"\nBilateral concordance summary saved to {output_path}")
    
    # Print top 10 most concordant regions
    print("\n" + "="*110)
    print("TOP 10 MOST CONCORDANT BILATERAL REGIONS")
    print("="*110)
    print(f"{'Region':<10} {'Direction':<12} {'Effect Size':<12} {'Concordance':<12} {'Score':<12} {'Total Voxels':<12}")
    print("-"*110)
    
    for idx, row in output_df.head(10).iterrows():
        print(f"{row['Region']:<10} {row['Direction']:<12} "
              f"{row['Effect Size (Absolute %)']:>10.2f}% "
              f"{row['Bilateral Concordance']:>10.3f} "
              f"{row['Concordance × Effect Size']:>10.2f} "
              f"{row['Total Significant Voxels']:>12}")

# =============================================================================
# Main Analysis
# =============================================================================

def main():
    print("Loading left and right hemisphere data...")
    
    # Load data
    left_df = load_hemisphere_data(LEFT_FILE)
    right_df = load_hemisphere_data(RIGHT_FILE)
    
    print(f"Left hemisphere: {len(left_df)} regions")
    print(f"Right hemisphere: {len(right_df)} regions")
    
    # Calculate bilateral metrics
    print("\nCalculating bilateral concordance metrics...")
    bilateral_df = calculate_bilateral_metrics(left_df, right_df)
    
    # Create concordance plot
    print("\nCreating bilateral concordance plot...")
    create_bilateral_concordance_plot(bilateral_df, OUTPUT_PLOT, 
                                     max_labels_decreased=MAX_LABELS_DECREASED,
                                     max_labels_increased=MAX_LABELS_INCREASED)
    
    # Save summary table
    print("\nSaving bilateral concordance summary...")
    save_bilateral_summary(bilateral_df, OUTPUT_TABLE)
    
    print("\n" + "="*100)
    print("ANALYSIS COMPLETE!")
    print("="*100)
    print(f"\nOutputs generated:")
    print(f"  - Plot: {OUTPUT_PLOT}")
    print(f"  - Table: {OUTPUT_TABLE}")

if __name__ == "__main__":
    main()