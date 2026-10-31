"""
Standalone legend figure describing interaction-evidence categories
(in the style of a protein-protein interaction network legend, e.g.
STRING-database-style edge colors).

This figure has no data dependency -- it is a fixed, editable reference
legend. Edit ``LEGEND_SECTIONS`` below to change the categories, colors,
or labels.

Usage
-----
    python evidence_legend.py --output-dir results/figures/evidence_legend
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Legend content: (section title, [(line color, label), ...])
# ---------------------------------------------------------------------------
LEGEND_SECTIONS = [
    ("Known Interactions", [
        ("#00c8ff", "From Curated Databases"),
        ("#ff00ff", "Experimentally Determined"),
    ]),
    ("Predicted Interactions", [
        ("#00a000", "Gene Neighborhood"),
        ("#ff0000", "Gene Fusions"),
        ("#0000cc", "Gene co-Occurrence"),
    ]),
    ("Other", [
        ("#a0b000", "Text Mining"),
        ("#000000", "Co-Expression"),
        ("#9090ff", "Protein Homology"),
    ]),
]

FONT_FAMILY_CANDIDATES = ["Times New Roman", "Liberation Serif", "Times"]
LINE_GAP = 0.09
SECTION_GAP = 0.05


def draw_legend(sections: list[tuple[str, list[tuple[str, str]]]]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(3.2, 4.5), dpi=300)
    ax.axis("off")

    y = 1.0
    for title, items in sections:
        ax.text(0, y, title, fontsize=11, fontweight="bold", va="top", ha="left")
        y -= LINE_GAP

        for color, label in items:
            x_start, x_end = 0.02, 0.16
            ax.add_line(Line2D([x_start, x_end], [y, y], color=color, linewidth=2.5))
            ax.plot(x_start, y, marker="o", markersize=7, markerfacecolor="white", markeredgecolor="black", markeredgewidth=0.8)
            ax.plot(x_end, y, marker="o", markersize=7, markerfacecolor="white", markeredgecolor="black", markeredgewidth=0.8)
            ax.text(x_end + 0.03, y, label, fontsize=10, va="center", ha="left")
            y -= LINE_GAP

        y -= SECTION_GAP

    ax.set_xlim(0, 1)
    ax.set_ylim(y, 1.05)
    fig.tight_layout()
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures/evidence_legend"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = FONT_FAMILY_CANDIDATES

    fig = draw_legend(LEGEND_SECTIONS)

    output_path = args.output_dir / "interaction_evidence_legend.png"
    fig.savefig(output_path, dpi=300, transparent=True, bbox_inches="tight")
    plt.show()

    print(f"Legend saved to: {output_path}")


if __name__ == "__main__":
    main()
