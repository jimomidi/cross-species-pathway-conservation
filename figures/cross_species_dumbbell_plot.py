"""
Dumbbell plot comparing per-pathway fold enrichment between two species,
one row per conserved pathway, with an optional database-abbreviation
column.

Input: the three "conserved shared pathways" CSVs produced by
``pathway_enrichment/compare_enrichment_results.py``.

Usage
-----
    python cross_species_dumbbell_plot.py \\
        --input-dir results/pathway_conservation/tables \\
        --output-dir results/figures/cross_species_dumbbell_plot
"""

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import DATABASE_COLORS, SPECIES_A, SPECIES_B  # noqa: E402
from plotting_utils import apply_publication_style, load_conserved_pathways, save_figure  # noqa: E402

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
USE_LOG2_SCALE = True  # recommended: fold-enrichment ranges often differ substantially

# "ratio", "difference", "mean_fold", or f"{SPECIES_A.lower()}_fold"
SORT_METHOD = "ratio"

# None = plot all pathways; an integer keeps only the top N (by SORT_METHOD).
TOP_N = None

PATHWAY_WRAP_WIDTH = 42
FIGURE_WIDTH, MIN_FIGURE_HEIGHT, HEIGHT_PER_PATHWAY = 11.5, 7.0, 0.46
POINT_SIZE, POINT_EDGE_WIDTH, LINE_WIDTH, POINT_ALPHA = 95, 1.1, 2.0, 0.92
SHOW_VALUES = True
SHOW_DATABASE_COLUMN = True

SPECIES_A_COLOR = "#76B7E5"
SPECIES_B_COLOR = "#D8A0C8"


def prepare_data(combined_df, species_a: str, species_b: str):
    a_col, b_col = f"{species_a} Fold Enrichment", f"{species_b} Fold Enrichment"

    combined_df["Fold Difference"] = combined_df[a_col] - combined_df[b_col]
    combined_df["Absolute Fold Difference"] = combined_df["Fold Difference"].abs()
    combined_df[f"Fold Enrichment Ratio {species_a}/{species_b}"] = combined_df[a_col] / combined_df[b_col]
    combined_df[f"Log2 Fold Ratio {species_a}/{species_b}"] = np.log2(
        combined_df[f"Fold Enrichment Ratio {species_a}/{species_b}"]
    )

    sort_columns = {
        "ratio": f"Fold Enrichment Ratio {species_a}/{species_b}",
        "difference": "Absolute Fold Difference",
        f"{species_a.lower()}_fold": a_col,
        "mean_fold": "Mean Fold Enrichment",
    }
    if SORT_METHOD not in sort_columns:
        raise ValueError(f"SORT_METHOD must be one of: {list(sort_columns)}")
    combined_df = combined_df.sort_values(by=sort_columns[SORT_METHOD], ascending=True)

    if TOP_N is not None:
        if not isinstance(TOP_N, int) or TOP_N <= 0:
            raise ValueError("TOP_N must be None or a positive integer.")
        combined_df = combined_df.tail(TOP_N)

    combined_df = combined_df.reset_index(drop=True)
    combined_df["Pathway Wrapped"] = combined_df["Pathway"].apply(
        lambda p: "\n".join(textwrap.wrap(str(p), width=PATHWAY_WRAP_WIDTH))
    )

    if USE_LOG2_SCALE:
        combined_df[f"{species_a} Plot Value"] = np.log2(combined_df[a_col])
        combined_df[f"{species_b} Plot Value"] = np.log2(combined_df[b_col])
        x_axis_label = "Fold Enrichment (log$_2$ scale)"
    else:
        combined_df[f"{species_a} Plot Value"] = combined_df[a_col]
        combined_df[f"{species_b} Plot Value"] = combined_df[b_col]
        x_axis_label = "Fold Enrichment"

    return combined_df, x_axis_label


def make_dumbbell_plot(combined_df, species_a: str, species_b: str, x_axis_label: str):
    a_col, b_col = f"{species_a} Fold Enrichment", f"{species_b} Fold Enrichment"
    a_plot, b_plot = f"{species_a} Plot Value", f"{species_b} Plot Value"

    number_of_pathways = len(combined_df)
    figure_height = max(MIN_FIGURE_HEIGHT, number_of_pathways * HEIGHT_PER_PATHWAY + 2.0)
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, figure_height))

    y_positions = np.arange(number_of_pathways)

    for position, (_, row) in zip(y_positions, combined_df.iterrows()):
        ax.plot([row[b_plot], row[a_plot]], [position, position], color="#A6A6A6", linewidth=LINE_WIDTH, alpha=0.78, solid_capstyle="round", zorder=1)

    ax.scatter(combined_df[b_plot], y_positions, s=POINT_SIZE, marker="o", c=SPECIES_B_COLOR, edgecolors="black", linewidths=POINT_EDGE_WIDTH, alpha=POINT_ALPHA, label=species_b, zorder=3)
    ax.scatter(combined_df[a_plot], y_positions, s=POINT_SIZE, marker="o", c=SPECIES_A_COLOR, edgecolors="black", linewidths=POINT_EDGE_WIDTH, alpha=POINT_ALPHA, label=species_a, zorder=4)

    if SHOW_VALUES:
        all_plot_values = np.concatenate([combined_df[a_plot].to_numpy(), combined_df[b_plot].to_numpy()])
        value_offset = max((np.nanmax(all_plot_values) - np.nanmin(all_plot_values)) * 0.018, 0.04)

        for position, (_, row) in zip(y_positions, combined_df.iterrows()):
            ax.text(row[b_plot] - value_offset, position, f"{row[b_col]:.2f}", ha="right", va="center", fontsize=8.5, zorder=5)
            ax.text(row[a_plot] + value_offset, position, f"{row[a_col]:.2f}", ha="left", va="center", fontsize=8.5, zorder=5)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(combined_df["Pathway Wrapped"], fontsize=10)

    if SHOW_DATABASE_COLUMN:
        x_min = min(combined_df[b_plot].min(), combined_df[a_plot].min())
        x_max = max(combined_df[b_plot].max(), combined_df[a_plot].max())
        database_x_position = x_min - max((x_max - x_min) * 0.20, 0.70)

        for position, (_, row) in zip(y_positions, combined_df.iterrows()):
            database_name = row["Database"]
            ax.text(
                database_x_position, position, database_name, ha="center", va="center", fontsize=8.8, fontweight="bold",
                color=DATABASE_COLORS[database_name],
                bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": DATABASE_COLORS[database_name], "linewidth": 0.9, "alpha": 0.95},
                clip_on=False,
            )

    ax.axvline(x=0 if USE_LOG2_SCALE else 1, linestyle=(0, (6, 4)), linewidth=1.2, color="black", alpha=0.65, zorder=0)

    ax.set_xlabel(x_axis_label, fontweight="bold", labelpad=12)
    ax.set_ylabel("")
    ax.set_title("Cross-Species Comparison of Conserved Pathway Enrichment", fontweight="bold", pad=18)

    ax.grid(axis="x", linestyle="--", linewidth=0.65, alpha=0.35, zorder=0)
    ax.grid(axis="y", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.1)
    ax.spines["bottom"].set_linewidth(1.1)
    ax.tick_params(axis="x", which="major", width=1.0, length=5)
    ax.tick_params(axis="y", which="major", length=0, pad=8)

    main_legend = ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=SPECIES_A_COLOR, markeredgecolor="black", markersize=9, label=species_a),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=SPECIES_B_COLOR, markeredgecolor="black", markersize=9, label=species_b),
        ],
        loc="lower right", frameon=True, edgecolor="black", fancybox=False,
    )
    ax.add_artist(main_legend)

    if SHOW_DATABASE_COLUMN:
        database_handles = [
            Line2D([0], [0], marker="s", color="none", markerfacecolor="white", markeredgecolor=DATABASE_COLORS[name], markeredgewidth=1.4, markersize=8, label=name)
            for name in DATABASE_COLORS
        ]
        ax.legend(handles=database_handles, loc="upper right", frameon=True, edgecolor="black", fancybox=False, title="Database")

    all_plot_values = np.concatenate([combined_df[a_plot].to_numpy(), combined_df[b_plot].to_numpy()])
    plot_range = np.nanmax(all_plot_values) - np.nanmin(all_plot_values)
    ax.set_xlim(np.nanmin(all_plot_values) - max(plot_range * 0.10, 0.35), np.nanmax(all_plot_values) + max(plot_range * 0.16, 0.55))

    plt.subplots_adjust(left=0.39 if SHOW_DATABASE_COLUMN else 0.34, right=0.97, top=0.93, bottom=0.08)

    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("results/pathway_conservation/tables"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures/cross_species_dumbbell_plot"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    apply_publication_style()
    plt.rcParams["font.size"] = 11
    plt.rcParams["axes.labelsize"] = 14
    plt.rcParams["xtick.labelsize"] = 10.5
    plt.rcParams["ytick.labelsize"] = 10.5

    combined_df = load_conserved_pathways(args.input_dir, SPECIES_A, SPECIES_B)
    combined_df, x_axis_label = prepare_data(combined_df, SPECIES_A, SPECIES_B)

    fig = make_dumbbell_plot(combined_df, SPECIES_A, SPECIES_B, x_axis_label)

    data_output_path = args.output_dir / "conserved_pathways_used_for_dumbbell_plot.csv"
    columns_to_save = [
        "Database", "Pathway", f"{SPECIES_A} Fold Enrichment", f"{SPECIES_B} Fold Enrichment",
        "Fold Difference", f"Fold Enrichment Ratio {SPECIES_A}/{SPECIES_B}",
        f"Log2 Fold Ratio {SPECIES_A}/{SPECIES_B}", "Mean Fold Enrichment",
    ]
    combined_df[columns_to_save].to_csv(data_output_path, index=False)

    written = save_figure(fig, args.output_dir / "cross_species_pathway_enrichment_dumbbell", formats=("png", "tiff", "pdf", "svg"))
    plt.show()

    print("=" * 75)
    print("DUMBBELL PLOT COMPLETED")
    print("=" * 75)
    print(f"\nTotal pathways plotted: {len(combined_df)}")
    print("\nPathways by database:")
    print(combined_df["Database"].value_counts().to_string())
    print(f"\nGenerated files:\n" + "\n".join(str(p) for p in written) + f"\n{data_output_path}")


if __name__ == "__main__":
    main()
