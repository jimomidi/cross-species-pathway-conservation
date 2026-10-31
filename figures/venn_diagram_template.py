"""
Blank two-circle Venn diagram template (overlapping-circle outline only,
no counts or labels). Useful as a background layer for manually annotated
figures, or as a starting point for a proportional/area-scaled Venn diagram.

This figure has no data dependency; it draws two fixed-size overlapping
circles in a chosen color pair.

Usage
-----
    python venn_diagram_template.py --output-dir results/figures/venn_template \\
        --color-1 "#A7D8F0" --color-2 "#F8C8DC"
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

CIRCLE_RADIUS = 1.5
CIRCLE_OFFSET = 0.8


def draw_two_circle_template(color_1: str, color_2: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(6, 6))

    circle_1 = Circle((-CIRCLE_OFFSET, 0), CIRCLE_RADIUS, facecolor=color_1, edgecolor="black", alpha=0.8, linewidth=2)
    circle_2 = Circle((CIRCLE_OFFSET, 0), CIRCLE_RADIUS, facecolor=color_2, edgecolor="black", alpha=0.8, linewidth=2)
    ax.add_patch(circle_1)
    ax.add_patch(circle_2)

    ax.set_xlim(-3, 3)
    ax.set_ylim(-2.2, 2.2)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.tight_layout()
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures/venn_template"))
    parser.add_argument("--color-1", default="#A7D8F0", help="Fill color for the left circle.")
    parser.add_argument("--color-2", default="#F8C8DC", help="Fill color for the right circle.")
    parser.add_argument("--output-name", default="venn_diagram_template", help="Output file stem.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fig = draw_two_circle_template(args.color_1, args.color_2)

    output_path = args.output_dir / f"{args.output_name}.png"
    fig.savefig(output_path, dpi=300, transparent=True)
    plt.show()

    print(f"Venn diagram template saved to: {output_path}")


if __name__ == "__main__":
    main()
