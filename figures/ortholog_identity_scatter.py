"""
Scatter plot comparing reciprocal ortholog sequence-identity percentages
for one-to-one orthologs between the two species.

Input: the "one-to-one orthologs" CSV produced by
``data_preparation/map_orthologs.py`` (columns
"Human sequence identity to dog (%)" and "Dog sequence identity to
human (%)" -- despite the column names, these represent the reciprocal
identity of species A relative to species B and vice versa; see
``config.py`` if you have renamed the species).

Usage
-----
    python ortholog_identity_scatter.py \\
        --input results/ortholog_mapping/one_to_one_orthologs.csv \\
        --output-dir results/figures/ortholog_identity_scatter
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import FONT_FAMILY  # noqa: E402
from plotting_utils import save_figure  # noqa: E402

IDENTITY_COLUMN_X = "Human sequence identity to dog (%)"
IDENTITY_COLUMN_Y = "Dog sequence identity to human (%)"

POINT_SIZE = 85
POINT_COLOR = "#A78BFA"
POINT_ALPHA = 0.75


def make_scatter_plot(df: pd.DataFrame, x_column: str, y_column: str) -> plt.Figure:
    x, y = df[x_column], df[y_column]

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x, y, s=POINT_SIZE, color=POINT_COLOR, alpha=POINT_ALPHA, edgecolors="black", linewidth=1.5)

    min_val = min(x.min(), y.min())
    max_val = max(x.max(), y.max())
    ax.plot([min_val, max_val], [min_val, max_val], color="black", linewidth=0.8, dashes=(10, 6), label="y = x")

    padding = 2
    ax.set_xlim(min_val - padding, max_val + padding)
    ax.set_ylim(min_val - padding, max_val + padding)

    ax.set_xlabel(x_column, fontsize=15, fontweight="bold", labelpad=10)
    ax.set_ylabel(y_column, fontsize=15, fontweight="bold", labelpad=10)
    ax.set_title("Ortholog Sequence Identity", fontsize=16, fontweight="bold", pad=8)
    ax.tick_params(labelsize=12)

    fig.tight_layout()
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="One-to-one orthologs CSV.")
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures/ortholog_identity_scatter"))
    parser.add_argument("--x-column", default=IDENTITY_COLUMN_X)
    parser.add_argument("--y-column", default=IDENTITY_COLUMN_Y)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams["font.family"] = FONT_FAMILY

    df = pd.read_csv(args.input)
    for column in (args.x_column, args.y_column):
        if column not in df.columns:
            raise KeyError(f"Column '{column}' was not found in {args.input}. Available: {df.columns.tolist()}")
    df[args.x_column] = pd.to_numeric(df[args.x_column], errors="coerce")
    df[args.y_column] = pd.to_numeric(df[args.y_column], errors="coerce")
    df = df.dropna(subset=[args.x_column, args.y_column])

    fig = make_scatter_plot(df, args.x_column, args.y_column)
    written = save_figure(fig, args.output_dir / "ortholog_sequence_identity_scatter", formats=("png", "pdf", "svg"))
    plt.show()

    print(f"Ortholog pairs plotted: {len(df)}")
    print("Generated files:\n" + "\n".join(str(p) for p in written))


if __name__ == "__main__":
    main()
