"""
Two complementary views of conserved pathway enrichment: a fold-enrichment
heatmap (one row per pathway, one column per species) and a lollipop plot
of the log2 fold-enrichment ratio between species.

Input: the three "conserved shared pathways" CSVs produced by
``pathway_enrichment/compare_enrichment_results.py``.

Usage
-----
    python heatmap_lollipop_plot.py \\
        --input-dir results/pathway_conservation/tables \\
        --output-dir results/figures/heatmap_lollipop_plot
"""

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import DATABASE_COLORS, SPECIES_A, SPECIES_B  # noqa: E402
from plotting_utils import apply_publication_style, load_conserved_pathways, save_figure  # noqa: E402

# None = use all pathways; an integer keeps the top N by fold-enrichment ratio.
TOP_N_PATHWAYS = None

# "ratio", "<species_a>_fold" (lowercased SPECIES_A), "mean_fold", or "database"
SORT_METHOD = "ratio"
PATHWAY_WRAP_WIDTH = 42
SHOW_HEATMAP_VALUES = True
SHOW_LOLLIPOP_VALUES = True
SHOW_DATABASE_LABELS = True

HEATMAP_CMAP = LinearSegmentedColormap.from_list("custom_enrichment", ["#F6F6F6", "#DCEAF4", "#A8CEE5", "#6AAED6", "#2F6F9F"])


def prepare_data(combined_df, species_a: str, species_b: str):
    a_col, b_col = f"{species_a} Fold Enrichment", f"{species_b} Fold Enrichment"

    combined_df["Mean Fold Enrichment"] = combined_df[[a_col, b_col]].mean(axis=1)
    combined_df[f"{species_a}/{species_b} Fold Ratio"] = combined_df[a_col] / combined_df[b_col]
    combined_df[f"Log2 {species_a}/{species_b} Fold Ratio"] = np.log2(combined_df[f"{species_a}/{species_b} Fold Ratio"])
    combined_df[f"{species_a} log2 Fold Enrichment"] = np.log2(combined_df[a_col])
    combined_df[f"{species_b} log2 Fold Enrichment"] = np.log2(combined_df[b_col])

    if TOP_N_PATHWAYS is not None:
        if not isinstance(TOP_N_PATHWAYS, int) or TOP_N_PATHWAYS <= 0:
            raise ValueError("TOP_N_PATHWAYS must be None or a positive integer.")
        combined_df = combined_df.sort_values(f"{species_a}/{species_b} Fold Ratio", ascending=False).head(TOP_N_PATHWAYS).copy()

    database_order = {name: i for i, name in enumerate(DATABASE_COLORS)}
    combined_df["Database Order"] = combined_df["Database"].map(database_order)

    sort_columns = {
        "ratio": (f"Log2 {species_a}/{species_b} Fold Ratio", False),
        f"{species_a.lower()}_fold": (a_col, False),
        "mean_fold": ("Mean Fold Enrichment", False),
    }
    if SORT_METHOD in sort_columns:
        column, ascending = sort_columns[SORT_METHOD]
        combined_df = combined_df.sort_values(column, ascending=ascending)
    elif SORT_METHOD == "database":
        combined_df = combined_df.sort_values(["Database Order", a_col], ascending=[True, False])
    else:
        raise ValueError(f"SORT_METHOD must be one of: {list(sort_columns) + ['database']}")

    combined_df = combined_df.reset_index(drop=True)
    combined_df["Pathway Wrapped"] = combined_df["Pathway"].apply(lambda p: "\n".join(textwrap.wrap(str(p), width=PATHWAY_WRAP_WIDTH)))

    return combined_df


def make_heatmap(plot_df, species_a: str, species_b: str):
    a_col, b_col = f"{species_a} Fold Enrichment", f"{species_b} Fold Enrichment"
    a_log, b_log = f"{species_a} log2 Fold Enrichment", f"{species_b} log2 Fold Enrichment"

    n = len(plot_df)
    fig, ax = plt.subplots(figsize=(8.8, max(6.5, n * 0.43 + 2.0)))

    heatmap_matrix = plot_df[[a_log, b_log]].to_numpy()
    image = ax.imshow(heatmap_matrix, aspect="auto", cmap=HEATMAP_CMAP, interpolation="nearest")

    ax.set_xticks([0, 1])
    ax.set_xticklabels([species_a, species_b], fontsize=12, fontweight="bold")
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(plot_df["Pathway Wrapped"], fontsize=9.5)

    if SHOW_HEATMAP_VALUES:
        color_threshold = (np.nanmin(heatmap_matrix) + np.nanmax(heatmap_matrix)) / 2
        for row_index in range(n):
            raw_values = [plot_df.loc[row_index, a_col], plot_df.loc[row_index, b_col]]
            for column_index in range(2):
                text_color = "white" if heatmap_matrix[row_index, column_index] > color_threshold else "black"
                ax.text(column_index, row_index, f"{raw_values[column_index]:.2f}", ha="center", va="center", fontsize=8.7, color=text_color, fontweight="bold")

    if SHOW_DATABASE_LABELS:
        for row_index, database_name in enumerate(plot_df["Database"]):
            ax.scatter(-0.62, row_index, s=70, marker="s", color=DATABASE_COLORS[database_name], edgecolor="black", linewidth=0.7, clip_on=False, zorder=4)

    ax.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=2.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="y", length=0, pad=12)

    colorbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.04)
    colorbar.set_label("Fold Enrichment (log$_2$)", fontsize=11, fontweight="bold")

    database_handles = [Patch(facecolor=DATABASE_COLORS[d], edgecolor="black", label=d) for d in DATABASE_COLORS]
    ax.legend(handles=database_handles, title="Database", loc="upper left", bbox_to_anchor=(1.16, 1.0), frameon=True, edgecolor="black", fancybox=False)

    ax.set_title("Cross-Species Conserved Pathway Enrichment", fontsize=16, fontweight="bold", pad=15)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.subplots_adjust(left=0.47, right=0.78, top=0.94, bottom=0.06)
    return fig


def make_lollipop(plot_df, species_a: str, species_b: str):
    ratio_col = f"{species_a}/{species_b} Fold Ratio"
    log_ratio_col = f"Log2 {species_a}/{species_b} Fold Ratio"

    lollipop_df = plot_df.sort_values(log_ratio_col, ascending=True).reset_index(drop=True)
    n = len(lollipop_df)
    fig, ax = plt.subplots(figsize=(11.0, max(6.5, n * 0.43 + 2.0)))
    y_positions = np.arange(n)

    for row_index, row in lollipop_df.iterrows():
        value = row[log_ratio_col]
        ax.hlines(y=row_index, xmin=min(0, value), xmax=max(0, value), color=DATABASE_COLORS[row["Database"]], linewidth=2.2, alpha=0.70, zorder=1)

    for database_name in DATABASE_COLORS:
        subset = lollipop_df[lollipop_df["Database"] == database_name]
        if subset.empty:
            continue
        point_sizes = 70 + 20 * np.log2(subset["Mean Fold Enrichment"] + 1)
        ax.scatter(subset[log_ratio_col], subset.index, s=point_sizes, color=DATABASE_COLORS[database_name], edgecolor="black", linewidth=1.0, alpha=0.90, label=database_name, zorder=3)

    ax.axvline(x=0, color="black", linestyle=(0, (6, 4)), linewidth=1.3, alpha=0.75, zorder=0)

    if SHOW_LOLLIPOP_VALUES:
        x_values = lollipop_df[log_ratio_col].to_numpy()
        text_offset = max((np.nanmax(x_values) - np.nanmin(x_values)) * 0.025, 0.05)
        for row_index, row in lollipop_df.iterrows():
            value, ratio = row[log_ratio_col], row[ratio_col]
            text_x = value + text_offset if value >= 0 else value - text_offset
            ax.text(text_x, row_index, f"{ratio:.1f}×", ha="left" if value >= 0 else "right", va="center", fontsize=8.5)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(lollipop_df["Pathway Wrapped"], fontsize=9.5)
    ax.tick_params(axis="y", length=0, pad=8)

    ax.set_xlabel(f"Relative Enrichment, log$_2$({species_a}/{species_b})", fontweight="bold", labelpad=10)
    ax.set_ylabel("")
    ax.set_title("Relative Enrichment of Conserved Pathways", fontweight="bold", pad=15)

    ax.text(0.01, 1.01, f"Higher in {species_b.lower()}", transform=ax.transAxes, ha="left", va="bottom", fontsize=10, fontstyle="italic")
    ax.text(0.99, 1.01, f"Higher in {species_a.lower()}", transform=ax.transAxes, ha="right", va="bottom", fontsize=10, fontstyle="italic")

    ax.grid(axis="x", linestyle="--", linewidth=0.65, alpha=0.30)
    ax.grid(axis="y", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.1)

    legend_handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=DATABASE_COLORS[d], markeredgecolor="black", markersize=9, label=d) for d in DATABASE_COLORS]
    ax.legend(handles=legend_handles, title="Database", loc="lower right", frameon=True, edgecolor="black", fancybox=False)

    minimum_x, maximum_x = min(0, lollipop_df[log_ratio_col].min()), max(0, lollipop_df[log_ratio_col].max())
    x_span = maximum_x - minimum_x
    ax.set_xlim(minimum_x - max(0.15, x_span * 0.08), maximum_x + max(0.40, x_span * 0.18))

    plt.subplots_adjust(left=0.40, right=0.97, top=0.94, bottom=0.08)
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("results/pathway_conservation/tables"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures/heatmap_lollipop_plot"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    apply_publication_style()
    plt.rcParams["font.size"] = 11
    plt.rcParams["axes.titlesize"] = 17
    plt.rcParams["axes.labelsize"] = 14
    plt.rcParams["xtick.labelsize"] = 11
    plt.rcParams["ytick.labelsize"] = 10

    combined_df = load_conserved_pathways(args.input_dir, SPECIES_A, SPECIES_B)
    plot_df = prepare_data(combined_df, SPECIES_A, SPECIES_B)

    data_output_path = args.output_dir / "pathway_data_used_for_heatmap_and_lollipop.csv"
    columns_to_save = [
        "Database", "Pathway", f"{SPECIES_A} Fold Enrichment", f"{SPECIES_B} Fold Enrichment",
        f"{SPECIES_A} log2 Fold Enrichment", f"{SPECIES_B} log2 Fold Enrichment", "Mean Fold Enrichment",
        f"{SPECIES_A}/{SPECIES_B} Fold Ratio", f"Log2 {SPECIES_A}/{SPECIES_B} Fold Ratio",
    ]
    plot_df[columns_to_save].to_csv(data_output_path, index=False)

    heatmap_fig = make_heatmap(plot_df, SPECIES_A, SPECIES_B)
    heatmap_written = save_figure(heatmap_fig, args.output_dir / "cross_species_pathway_heatmap", formats=("png", "tiff", "pdf", "svg"))
    plt.show()

    lollipop_fig = make_lollipop(plot_df, SPECIES_A, SPECIES_B)
    lollipop_written = save_figure(lollipop_fig, args.output_dir / "cross_species_pathway_lollipop", formats=("png", "tiff", "pdf", "svg"))
    plt.show()

    print("=" * 80)
    print("HEATMAP AND LOLLIPOP PLOTS COMPLETED")
    print("=" * 80)
    print(f"\nTotal pathways plotted: {len(plot_df)}")
    print("\nPathways by database:")
    print(plot_df["Database"].value_counts().to_string())
    print("\nGenerated files:")
    for p in heatmap_written + lollipop_written:
        print(p)
    print(data_output_path)


if __name__ == "__main__":
    main()
