"""
Conceptual summary figure: conserved pathways grouped into broad biological
themes using transparent keyword-based classification, with representative
pathways and mean cross-species fold enrichment shown per theme.

Input: the three "conserved shared pathways" CSVs produced by
``pathway_enrichment/compare_enrichment_results.py`` (must include
``Dog Count`` -> ``<species> Count`` and ``<species> P-Value`` columns).

Usage
-----
    python conceptual_theme_summary.py \\
        --input-dir results/pathway_conservation/tables \\
        --output-dir results/figures/conceptual_theme_summary
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import SPECIES_A, SPECIES_B  # noqa: E402
from plotting_utils import apply_publication_style, load_conserved_pathways, save_figure  # noqa: E402

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
TOP_N_PATHWAYS = None  # None = use all pathways; an integer selects the top N before theme classification.
PATHWAY_SELECTION_METHOD = "mean_fold"  # "mean_fold", "<species_a>_fold", "ratio", or "min_pvalue"
MAX_PATHWAYS_PER_THEME = 10

CONCEPT_FIGURE_WIDTH = 20.0
ROW_HEIGHT = 2.65
FIGURE_VERTICAL_PADDING = 3.6
TITLE_FONT_SIZE, COLUMN_HEADER_FONT_SIZE, THEME_FONT_SIZE = 25, 17, 15
PATHWAY_FONT_SIZE, ENRICHMENT_FONT_SIZE, COUNT_FONT_SIZE, NOTE_FONT_SIZE = 12.2, 13.2, 11.2, 11.5
BOX_ROUNDING = 0.030

# Broad biological themes and the keywords used to classify pathway names
# into each one. This is a transparent, editable heuristic, provided only
# as a generic starting template -- replace the theme names and keywords
# with categories and terms that match the biology of your own pathway
# set before relying on this figure.
THEME_KEYWORDS = {
    "Theme A\n(edit me)": ["cytoskeleton", "cell movement"],
    "Theme B\n(edit me)": ["signaling", "receptor"],
    "Theme C\n(edit me)": ["dna repair", "chromatin"],
    "Theme D\n(edit me)": ["transcription", "rna"],
    "Theme E\n(edit me)": ["immune", "differentiation"],
}
OTHER_THEME = "Other conserved\nprocesses"
THEME_ORDER = list(THEME_KEYWORDS) + [OTHER_THEME]

THEME_COLORS = {
    "Theme A\n(edit me)": "#CFE8F6",
    "Theme B\n(edit me)": "#D5EEDC",
    "Theme C\n(edit me)": "#E8D9F1",
    "Theme D\n(edit me)": "#F7E3C8",
    "Theme E\n(edit me)": "#F3D7DE",
    OTHER_THEME: "#E4E7EB",
}
THEME_EDGE_COLORS = {
    "Theme A\n(edit me)": "#2E78A6",
    "Theme B\n(edit me)": "#31845A",
    "Theme C\n(edit me)": "#76508F",
    "Theme D\n(edit me)": "#B06E23",
    "Theme E\n(edit me)": "#A34B63",
    OTHER_THEME: "#65717D",
}

SPECIES_A_COLOR, SPECIES_B_COLOR = "#126A9C", "#A23B72"
PATHWAY_BOX_FACE, PATHWAY_BOX_EDGE = "#FFFFFF", "#6D7783"
ENRICHMENT_BOX_FACE, ENRICHMENT_BOX_EDGE = "#F7F9FB", "#394B59"
ARROW_COLOR_1, ARROW_COLOR_2 = "#417A9B", "#7A4E78"
HEADER_COLOR, NOTE_COLOR = "#1F2D3D", "#4F5964"


def assign_biological_theme(pathway_name: str) -> str:
    """Assign a pathway to a broad biological theme via keyword matching."""
    name = str(pathway_name).lower()
    for theme, keywords in THEME_KEYWORDS.items():
        if any(keyword in name for keyword in keywords):
            return theme
    return OTHER_THEME


def select_top_pathways(plot_df, species_a: str, species_b: str):
    if TOP_N_PATHWAYS is None:
        return plot_df

    if not isinstance(TOP_N_PATHWAYS, int) or TOP_N_PATHWAYS <= 0:
        raise ValueError("TOP_N_PATHWAYS must be None or a positive integer.")

    selection_columns = {
        "mean_fold": ("Mean Fold Enrichment", False),
        f"{species_a.lower()}_fold": (f"{species_a} Fold Enrichment", False),
        "ratio": (f"{species_a}/{species_b} Fold Ratio", False),
        "min_pvalue": ("Minimum P-Value", True),
    }
    if PATHWAY_SELECTION_METHOD not in selection_columns:
        raise ValueError(f"PATHWAY_SELECTION_METHOD must be one of: {list(selection_columns)}")

    column, ascending = selection_columns[PATHWAY_SELECTION_METHOD]
    return plot_df.sort_values(column, ascending=ascending).head(TOP_N_PATHWAYS).copy()


def build_theme_records(concept_df, species_a: str, species_b: str):
    records = []
    for theme in THEME_ORDER:
        subset = concept_df[concept_df["Biological Theme"] == theme]
        if subset.empty:
            continue
        subset = subset.sort_values("Mean Fold Enrichment", ascending=False)
        records.append(
            {
                "Theme": theme,
                "Pathways": subset.head(MAX_PATHWAYS_PER_THEME),
                "Total Pathways": len(subset),
                f"Mean {species_a} Fold": subset[f"{species_a} Fold Enrichment"].mean(),
                f"Mean {species_b} Fold": subset[f"{species_b} Fold Enrichment"].mean(),
            }
        )
    return records


def draw_conceptual_summary(concept_records, species_a: str, species_b: str):
    n_themes = len(concept_records)
    fig, ax = plt.subplots(figsize=(CONCEPT_FIGURE_WIDTH, n_themes * ROW_HEIGHT + FIGURE_VERTICAL_PADDING), facecolor="white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n_themes + 1.35)
    ax.axis("off")

    header_y = n_themes + 0.78
    for x, text in [(0.145, "Conserved biological theme"), (0.535, "Representative shared pathways"), (0.875, "Cross-species enrichment")]:
        ax.text(x, header_y, text, ha="center", va="center", fontsize=COLUMN_HEADER_FONT_SIZE, fontweight="bold", color=HEADER_COLOR)

    for index, record in enumerate(concept_records):
        y = n_themes - index - 0.02
        theme = record["Theme"]

        theme_x, theme_width, theme_height = 0.018, 0.255, 1.30
        pathway_x, pathway_width, pathway_height = 0.325, 0.445, 1.55
        enrichment_x, enrichment_width, enrichment_height = 0.825, 0.155, 1.18

        ax.add_patch(FancyBboxPatch((theme_x, y - theme_height / 2), theme_width, theme_height, boxstyle=f"round,pad=0.025,rounding_size={BOX_ROUNDING}", facecolor=THEME_COLORS[theme], edgecolor=THEME_EDGE_COLORS[theme], linewidth=2.0, zorder=2))
        ax.text(theme_x + theme_width / 2, y, theme, ha="center", va="center", fontsize=THEME_FONT_SIZE, fontweight="bold", color="#17212B", linespacing=1.25, zorder=3)

        pathway_text = "\n".join(f"• {row['Pathway']} [{row['Database']}]" for _, row in record["Pathways"].iterrows())
        ax.add_patch(FancyBboxPatch((pathway_x, y - pathway_height / 2), pathway_width, pathway_height, boxstyle=f"round,pad=0.022,rounding_size={BOX_ROUNDING}", facecolor=PATHWAY_BOX_FACE, edgecolor=PATHWAY_BOX_EDGE, linewidth=1.6, zorder=2))
        ax.text(pathway_x + 0.018, y, pathway_text, ha="left", va="center", fontsize=PATHWAY_FONT_SIZE, color="#20262D", linespacing=1.32, zorder=3)

        ax.add_patch(FancyBboxPatch((enrichment_x, y - enrichment_height / 2), enrichment_width, enrichment_height, boxstyle=f"round,pad=0.024,rounding_size={BOX_ROUNDING}", facecolor=ENRICHMENT_BOX_FACE, edgecolor=ENRICHMENT_BOX_EDGE, linewidth=1.8, zorder=2))
        ax.text(enrichment_x + enrichment_width / 2, y + 0.27, f"{species_a}: {record[f'Mean {species_a} Fold']:.2f}", ha="center", va="center", fontsize=ENRICHMENT_FONT_SIZE, fontweight="bold", color=SPECIES_A_COLOR, zorder=3)
        ax.text(enrichment_x + enrichment_width / 2, y - 0.04, f"{species_b}: {record[f'Mean {species_b} Fold']:.2f}", ha="center", va="center", fontsize=ENRICHMENT_FONT_SIZE, fontweight="bold", color=SPECIES_B_COLOR, zorder=3)
        ax.text(enrichment_x + enrichment_width / 2, y - 0.37, f"n = {record['Total Pathways']} pathways", ha="center", va="center", fontsize=COUNT_FONT_SIZE, fontweight="bold", color="#3F4852", zorder=3)

        ax.add_patch(FancyArrowPatch((theme_x + theme_width + 0.008, y), (pathway_x - 0.010, y), arrowstyle="-|>", mutation_scale=24, linewidth=2.6, color=ARROW_COLOR_1, shrinkA=0, shrinkB=0, zorder=4))
        ax.add_patch(FancyArrowPatch((pathway_x + pathway_width + 0.008, y), (enrichment_x - 0.010, y), arrowstyle="-|>", mutation_scale=24, linewidth=2.6, color=ARROW_COLOR_2, shrinkA=0, shrinkB=0, zorder=4))

    ax.set_title("Conceptual Summary of Conserved Molecular Programs", fontsize=TITLE_FONT_SIZE, fontweight="bold", color=HEADER_COLOR, pad=28)
    ax.text(
        0.5, 0.10,
        "Pathways were grouped into broad biological themes using predefined keyword-based classification. "
        "Values represent mean fold enrichment within each theme.",
        ha="center", va="center", fontsize=NOTE_FONT_SIZE, fontstyle="italic", color=NOTE_COLOR,
    )

    plt.subplots_adjust(left=0.025, right=0.985, top=0.925, bottom=0.055)
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("results/pathway_conservation/tables"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures/conceptual_theme_summary"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    apply_publication_style()
    plt.rcParams["axes.titlesize"] = TITLE_FONT_SIZE

    combined_df = load_conserved_pathways(args.input_dir, SPECIES_A, SPECIES_B)

    required = [f"{SPECIES_A} Count", f"{SPECIES_B} Count", f"{SPECIES_A} P-Value", f"{SPECIES_B} P-Value"]
    missing = [c for c in required if c not in combined_df.columns]
    if missing:
        raise KeyError(f"Missing columns required for this figure: {missing}")

    combined_df[f"{SPECIES_A}/{SPECIES_B} Fold Ratio"] = combined_df[f"{SPECIES_A} Fold Enrichment"] / combined_df[f"{SPECIES_B} Fold Enrichment"]
    combined_df["Minimum P-Value"] = combined_df[[f"{SPECIES_A} P-Value", f"{SPECIES_B} P-Value"]].min(axis=1)

    plot_df = select_top_pathways(combined_df, SPECIES_A, SPECIES_B).reset_index(drop=True)

    concept_df = plot_df.copy()
    concept_df["Biological Theme"] = concept_df["Pathway"].apply(assign_biological_theme)

    concept_records = build_theme_records(concept_df, SPECIES_A, SPECIES_B)
    if not concept_records:
        raise ValueError("No pathways remained for the conceptual figure.")

    fig = draw_conceptual_summary(concept_records, SPECIES_A, SPECIES_B)
    written = save_figure(fig, args.output_dir / "conceptual_pathway_summary", formats=("png", "tiff", "pdf", "svg"))
    plt.show()

    table_path = args.output_dir / "conceptual_pathway_theme_classification.csv"
    concept_df[
        ["Database", "Pathway", "Biological Theme", f"{SPECIES_A} Fold Enrichment", f"{SPECIES_B} Fold Enrichment",
         "Mean Fold Enrichment", f"{SPECIES_A} Count", f"{SPECIES_B} Count", f"{SPECIES_A} P-Value", f"{SPECIES_B} P-Value"]
    ].to_csv(table_path, index=False)

    print("=" * 85)
    print("CONCEPTUAL PATHWAY FIGURE COMPLETED")
    print("=" * 85)
    print(f"\nTotal pathways used: {len(plot_df)}")
    print("\nConceptual themes:")
    print(concept_df["Biological Theme"].value_counts().to_string())
    print(f"\nGenerated files:\n" + "\n".join(str(p) for p in written) + f"\n{table_path}")


if __name__ == "__main__":
    main()
