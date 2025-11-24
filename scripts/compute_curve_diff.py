"""
Script to compute curve differences for benchmarking.
Input: merged_df files from propr wilcoxon benchmark (cluster)
"""
#%%
# ------ IMPORT LIBRARIES AND VARIABLES ------

import pandas as pd
import os
import numpy as np
import glob
import sys
from pathlib import Path


# Add project root to sys.path so "config" can be imported

from propr_cluster.config.paths import MERGED_DF_DIR, RESULTS_DIR

 #%%
# ----- DEFINE FUNCTIONS -----
def build_cumulative_curve(sorted_data_frame):
    """Return cumulative count of significant genes (as a numpy array)."""
    significant_as_integer = sorted_data_frame["significant"].astype(int)
    cumulative_significant = significant_as_integer.cumsum().to_numpy()
    return cumulative_significant


def area_between_curves(curve_one, curve_two):
    """Return area between two cumulative curves (absolute vertical difference)."""
    if curve_one.shape != curve_two.shape:
        raise ValueError("Curves must have the same length to compare them")
    x_values = np.arange(1, curve_one.shape[0] + 1)
    vertical_difference = np.abs(curve_one - curve_two)
    area = np.trapz(vertical_difference, x_values)
    return float(area)


def process_single_file(file_path):
    """Compute all requested metrics for a single merged_all file."""
    data_frame = pd.read_csv(file_path)

    number_of_rows = data_frame.shape[0]

    # Targeted gene from file name: merged_all_<GENE>.csv
    file_name = os.path.basename(file_path)
    if file_name.startswith("merged_all_") and file_name.endswith(".csv"):
        targeted_gene = file_name[len("merged_all_") : -len(".csv")]
    else:
        targeted_gene = file_name

    # Build all curves from the same original data_frame

    # Black: Wilcoxon ground truth (WGT) – sorted by p_value ascending
    sorted_by_wilcoxon = data_frame.sort_values("p_value", ascending=True)
    curve_wgt = build_cumulative_curve(sorted_by_wilcoxon)
    number_de_wilcoxon = int(curve_wgt[-1])

    # Yellow: propd_single ground truth (PGT) – sorted by ranked_weighted_connectivity
    sorted_by_propd_single = data_frame.sort_values(
        "ranked_weighted_connectivity", ascending=True
    )
    curve_pgt = build_cumulative_curve(sorted_by_propd_single)

    # Blue: propd_combined – sorted by ranked_theta_combined
    sorted_by_blue = data_frame.sort_values(
        "ranked_theta_combined", ascending=True
    )
    curve_blue = build_cumulative_curve(sorted_by_blue)

    # Green: enrichment absolute – sorted by es_abs_rank
    sorted_by_green = data_frame.sort_values("es_abs_rank", ascending=True)
    curve_green = build_cumulative_curve(sorted_by_green)

    # Brown: enrichment positive – sorted by es_pos descending
    sorted_by_brown = data_frame.sort_values("es_pos", ascending=False)
    curve_brown = build_cumulative_curve(sorted_by_brown)

    # Purple: enrichment sign – sorted by es_rank
    sorted_by_purple = data_frame.sort_values("es_rank", ascending=True)
    curve_purple = build_cumulative_curve(sorted_by_purple)

    # Pink: enrichment negative – sorted by es_neg ascending
    sorted_by_pink = data_frame.sort_values("es_neg", ascending=True)
    curve_pink = build_cumulative_curve(sorted_by_pink)

    # Sanity check: all curves should have the same length
    expected_length = number_of_rows
    for name, curve in [
        ("WGT", curve_wgt),
        ("PGT", curve_pgt),
        ("blue", curve_blue),
        ("green", curve_green),
        ("brown", curve_brown),
        ("purple", curve_purple),
        ("pink", curve_pink),
    ]:
        if curve.shape[0] != expected_length:
            raise ValueError(f"Curve {name} in {file_name} has wrong length")

    # Areas versus Wilcoxon ground truth (black)
    blue_vs_wgt = area_between_curves(curve_blue, curve_wgt)
    green_vs_wgt = area_between_curves(curve_green, curve_wgt)
    brown_vs_wgt = area_between_curves(curve_brown, curve_wgt)
    purple_vs_wgt = area_between_curves(curve_purple, curve_wgt)
    pink_vs_wgt = area_between_curves(curve_pink, curve_wgt)

    # Area between Wilcoxon ground truth and propd_single ground truth
    wgt_vs_pgt = area_between_curves(curve_wgt, curve_pgt)

    # Areas versus propd_single ground truth (yellow)
    blue_vs_pgt = area_between_curves(curve_blue, curve_pgt)
    green_vs_pgt = area_between_curves(curve_green, curve_pgt)
    brown_vs_pgt = area_between_curves(curve_brown, curve_pgt)
    purple_vs_pgt = area_between_curves(curve_purple, curve_pgt)
    pink_vs_pgt = area_between_curves(curve_pink, curve_pgt)

    result_row = {
        "Targeted_gene": targeted_gene,
        "num_DE_Wilcoxon": number_de_wilcoxon,
        "blue_vs_WGT": blue_vs_wgt,
        "green_vs_WGT": green_vs_wgt,
        "brown_vs_WGT": brown_vs_wgt,
        "purple_vs_WGT": purple_vs_wgt,
        "pink_vs_WGT": pink_vs_wgt,
        "WGT_vs_PGT": wgt_vs_pgt,
        "blue_vs_PGT": blue_vs_pgt,
        "green_vs_PGT": green_vs_pgt,
        "brown_vs_PGT": brown_vs_pgt,
        "purple_vs_PGT": purple_vs_pgt,
        "pink_vs_PGT": pink_vs_pgt,
    }
    return result_row


# ---- MAIN SCRIPT EXECUTION -----

file_paths = list(MERGED_DF.glob("merged_all_*.csv"))
output_summary_path = os.path.join(RESULTS, "curve_differences_summary.csv")

if not file_paths:
    raise RuntimeError(f"No files matched pattern {MERGED_DF / 'merged_all_*.csv'}")

summary_rows = []
for f in file_paths:
    
    df = pd.read_csv(f)
    print("Loaded:", f)
    summary_row = process_single_file(f)
    summary_rows.append(summary_row)

summary_data_frame = pd.DataFrame(summary_rows)

# ---- SAVE OUTPUT TO CSV -----
summary_data_frame.to_csv(output_summary_path, index=False)
print(f"Wrote summary to {output_summary_path}")

