"""
Scatter plot comparing fold enrichment of conserved pathways between two
species, colored by source database (GO BP / KEGG / Reactome), with an
identity line and Pearson/Spearman correlation annotation.

Input: the three "conserved shared pathways" CSVs produced by
``pathway_enrichment/compare_enrichment_results.py``.

Usage
-----
    python cross_species_scatter_plot.py \\
        --input-dir results/pathway_conservation/tables \\
        --output-dir results/figures/cross_species_scatter_plot
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.stats import pearsonr, spearmanr

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import DATABASE_COLORS, DATABASE_MARKERS, SPECIES_A, SPECIES_B  # noqa: E402
from plotting_utils import apply_publication_style, load_conserved_pathways, save_figure  # noqa: E402

# ---------------------------------------------------------------------------
# Figure settings
# ---------------------------------------------------------------------------
USE_LOG_SCALE = False
SHOW_LABELS = True

# Pathway names to label on the plot (case-insensitive match against the
# "Pathway" column). Empty by default -- add the pathway names you want
# labeled for your own results, or set SHOW_LABELS = False to disable
# labeling entirely.
LABEL_PATHWAYS: list[str] = []

FIGURE_WIDTH, FIGURE_HEIGHT = 9.5, 8.0
POINT_SIZE, POINT_ALPHA, EDGE_WIDTH = 105, 0.82, 1.1


def make_scatter_plot(combined_df, species_a: str, species_b: str):
    a_col, b_col = f"{species_a} Fold Enrichment", f"{species_b} Fold Enrichment"

    a_values = combined_df[a_col].to_numpy()
    b_values = combined_df[b_col].to_numpy()

    if len(combined_df) >= 3:
        pearson_r, pearson_p = pearsonr(b_values, a_values)
        spearman_rho, spearman_p = spearmanr(b_values, a_values)
    else:
        pearson_r = pearson_p = spearman_rho = spearman_p = np.nan

    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))

    for database_name in DATABASE_COLORS:
        subset = combined_df[combined_df["Database"] == database_name]
        if subset.empty:
            continue
        ax.scatter(
            subset[b_col], subset[a_col], s=POINT_SIZE, marker=DATABASE_MARKERS[database_name],
            c=DATABASE_COLORS[database_name], edgecolors="black", linewidths=EDGE_WIDTH, alpha=POINT_ALPHA,
            label=f"{database_name} (n={len(subset)})", zorder=3,
        )

    all_values = np.concatenate([a_values, b_values])
    minimum_value, maximum_value = np.nanmin(all_values), np.nanmax(all_values)

    if USE_LOG_SCALE:
        lower_limit, upper_limit = max(minimum_value / 1.35, 0.1), maximum_value * 1.35
        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
    else:
        lower_limit, upper_limit = max(0, minimum_value - 0.05 * maximum_value), maximum_value * 1.10

    ax.set_xlim(lower_limit, upper_limit)
    ax.set_ylim(lower_limit, upper_limit)

    identity_values = np.linspace(lower_limit, upper_limit, 500)
    ax.plot(identity_values, identity_values, linestyle=(0, (8, 5)), linewidth=1.7, color="black", alpha=0.75, label="Equal enrichment", zorder=2)

    if SHOW_LABELS:
        normalized_labels = {p.strip().lower() for p in LABEL_PATHWAYS}
        labeled_rows = combined_df[combined_df["Pathway"].astype(str).str.strip().str.lower().isin(normalized_labels)]
        offsets = [(7, 7), (7, -13), (-8, 9), (-8, -14), (10, 2), (-10, 2)]

        for i, (_, row) in enumerate(labeled_rows.iterrows()):
            offset = offsets[i % len(offsets)]
            ax.annotate(
                str(row["Pathway"]), xy=(row[b_col], row[a_col]), xytext=offset, textcoords="offset points",
                fontsize=9.5, ha="left" if offset[0] >= 0 else "right", va="bottom" if offset[1] >= 0 else "top",
                arrowprops={"arrowstyle": "-", "linewidth": 0.65, "color": "gray", "alpha": 0.75}, zorder=4,
            )

    ax.set_xlabel(f"{species_b} Fold Enrichment", fontweight="bold", labelpad=10)
    ax.set_ylabel(f"{species_a} Fold Enrichment", fontweight="bold", labelpad=10)
    ax.set_title("Cross-Species Comparison of Conserved Pathway Enrichment", fontweight="bold", pad=15)

    if len(combined_df) >= 3:
        correlation_text = (
            f"Pearson r = {pearson_r:.3f} (P = {pearson_p:.2e})\n"
            f"Spearman ρ = {spearman_rho:.3f} (P = {spearman_p:.2e})"
        )
        ax.text(
            0.04, 0.96, correlation_text, transform=ax.transAxes, ha="left", va="top", fontsize=11,
            bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "edgecolor": "black", "linewidth": 0.8, "alpha": 0.90},
            zorder=5,
        )

    if USE_LOG_SCALE:
        formatter = ScalarFormatter()
        formatter.set_scientific(False)
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)

    ax.grid(True, which="major", linestyle="--", linewidth=0.65, alpha=0.35, zorder=1)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.45, alpha=0.20, zorder=1)
    ax.legend(loc="lower right", frameon=True, edgecolor="black", fancybox=False)

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    ax.tick_params(axis="both", which="major", width=1.1, length=6)
    ax.tick_params(axis="both", which="minor", width=0.8, length=3)

    plt.tight_layout()

    stats = {"pearson_r": pearson_r, "pearson_p": pearson_p, "spearman_rho": spearman_rho, "spearman_p": spearman_p}
    return fig, stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("results/pathway_conservation/tables"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures/cross_species_scatter_plot"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    apply_publication_style()

    combined_df = load_conserved_pathways(args.input_dir, SPECIES_A, SPECIES_B)

    fig, stats = make_scatter_plot(combined_df, SPECIES_A, SPECIES_B)

    combined_df.to_csv(args.output_dir / "combined_conserved_pathways_used_for_scatter.csv", index=False)
    written = save_figure(fig, args.output_dir / "cross_species_pathway_fold_enrichment_scatter")
    plt.show()

    print("=" * 70)
    print("SCATTER PLOT COMPLETED")
    print("=" * 70)
    print(f"\nTotal pathways plotted: {len(combined_df)}")
    print("\nPathways by database:")
    print(combined_df["Database"].value_counts())
    print(f"\nPearson r = {stats['pearson_r']:.4f}, P-value = {stats['pearson_p']:.4e}" if np.isfinite(stats["pearson_r"]) else "\nPearson correlation not calculable")
    print(f"Spearman rho = {stats['spearman_rho']:.4f}, P-value = {stats['spearman_p']:.4e}" if np.isfinite(stats["spearman_rho"]) else "Spearman correlation not calculable")
    print(f"\nGenerated files:\n" + "\n".join(str(p) for p in written))


if __name__ == "__main__":
    main()
