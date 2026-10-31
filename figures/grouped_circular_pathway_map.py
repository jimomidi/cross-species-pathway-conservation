"""
Circular network map of conserved pathways: the same gene-overlap network
as ``pathway_network_map.py``, but laid out as nodes grouped by database
around a single circle (a "chord diagram"-style view), with leader-line
labels and per-database arc labels.

Input: the three "conserved shared pathways" CSVs produced by
``pathway_enrichment/compare_enrichment_results.py``.

Usage
-----
    python grouped_circular_pathway_map.py \\
        --input-dir results/pathway_conservation/tables \\
        --output-dir results/figures/grouped_circular_pathway_map
"""

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import DATABASE_COLORS, DATABASE_MARKERS, SPECIES_A, SPECIES_B  # noqa: E402
from plotting_utils import apply_publication_style, load_conserved_pathways, save_figure, scale_values  # noqa: E402
from pathway_network_map import build_pathway_network  # noqa: E402

# ---------------------------------------------------------------------------
# Layout settings
# ---------------------------------------------------------------------------
CIRCLE_RADIUS = 6.3
START_ANGLE_DEGREES = 90
DATABASE_GROUP_GAP_DEGREES = 9
SORT_WITHIN_DATABASE = "mean_fold"  # "mean_fold", "degree", or "pathway"

LABEL_WRAP_WIDTH = 28
LABEL_FONT_SIZE = 8.7
LABEL_RADIUS_MULTIPLIER = 1.34
LEADER_LINE_RADIUS_MULTIPLIER = 1.25

MIN_NODE_SIZE, MAX_NODE_SIZE = 450, 1750
MIN_EDGE_WIDTH, MAX_EDGE_WIDTH = 0.45, 3.2
EDGE_COLOR, EDGE_ALPHA = "#969696", 0.28
NODE_EDGE_WIDTH = 1.25


def create_grouped_circular_layout(graph, database_order, radius, start_angle_degrees, group_gap_degrees, sort_method="mean_fold"):
    """Place nodes around a circle, grouped by database, in a chosen order within each group."""
    sort_keys = {
        "mean_fold": lambda n: graph.nodes[n]["mean_fold"],
        "degree": lambda n: graph.degree(n),
        "pathway": lambda n: graph.nodes[n]["pathway"].lower(),
    }
    if sort_method not in sort_keys:
        raise ValueError("sort_method must be 'mean_fold', 'degree', or 'pathway'.")

    reverse = sort_method in ("mean_fold", "degree")
    grouped_nodes = {
        db: sorted([n for n in graph.nodes if graph.nodes[n]["database"] == db], key=sort_keys[sort_method], reverse=reverse)
        for db in database_order
    }

    total_nodes = sum(len(nodes) for nodes in grouped_nodes.values())
    if total_nodes == 0:
        raise ValueError("No nodes are available for circular positioning.")

    nonempty_groups = [db for db in database_order if grouped_nodes[db]]
    available_degrees = 360.0 - group_gap_degrees * len(nonempty_groups)
    angle_per_node = available_degrees / total_nodes

    positions, node_angles = {}, {}
    current_angle = float(start_angle_degrees)

    for database_name in nonempty_groups:
        for node in grouped_nodes[database_name]:
            angle_radians = np.deg2rad(current_angle)
            positions[node] = np.array([radius * np.cos(angle_radians), radius * np.sin(angle_radians)])
            node_angles[node] = current_angle % 360
            current_angle -= angle_per_node
        current_angle -= group_gap_degrees

    return positions, node_angles, grouped_nodes


def draw_circular_map(graph, database_order):
    mean_fold_values = np.array([graph.nodes[n]["mean_fold"] for n in graph.nodes])
    node_sizes = dict(zip(graph.nodes, scale_values(np.log2(mean_fold_values + 1), MIN_NODE_SIZE, MAX_NODE_SIZE)))

    edges = list(graph.edges)
    edge_widths = scale_values([graph.edges[e]["jaccard"] for e in edges], MIN_EDGE_WIDTH, MAX_EDGE_WIDTH) if edges else []

    pos, node_angles, grouped_nodes = create_grouped_circular_layout(
        graph, database_order, CIRCLE_RADIUS, START_ANGLE_DEGREES, DATABASE_GROUP_GAP_DEGREES, SORT_WITHIN_DATABASE
    )

    fig, ax = plt.subplots(figsize=(18, 18))

    if edges:
        nx.draw_networkx_edges(graph, pos, ax=ax, edgelist=edges, width=edge_widths, edge_color=EDGE_COLOR, alpha=EDGE_ALPHA)

    for database_name in database_order:
        nodes = grouped_nodes.get(database_name, [])
        if not nodes:
            continue
        nx.draw_networkx_nodes(
            graph, pos, ax=ax, nodelist=nodes, node_size=[node_sizes[n] for n in nodes],
            node_color=DATABASE_COLORS[database_name], node_shape=DATABASE_MARKERS[database_name],
            edgecolors="black", linewidths=NODE_EDGE_WIDTH, alpha=0.92, label=database_name,
        )

    # Leader lines + labels, pointing radially outward from each node.
    for node in graph.nodes:
        node_x, node_y = pos[node]
        vector_length = np.hypot(node_x, node_y)
        unit_x, unit_y = (0, 1) if np.isclose(vector_length, 0) else (node_x / vector_length, node_y / vector_length)

        label_x, label_y = unit_x * CIRCLE_RADIUS * LABEL_RADIUS_MULTIPLIER, unit_y * CIRCLE_RADIUS * LABEL_RADIUS_MULTIPLIER
        leader_x, leader_y = unit_x * CIRCLE_RADIUS * LEADER_LINE_RADIUS_MULTIPLIER, unit_y * CIRCLE_RADIUS * LEADER_LINE_RADIUS_MULTIPLIER

        ax.plot([node_x, leader_x], [node_y, leader_y], color="#8B8B8B", linewidth=0.75, alpha=0.65, zorder=1)

        ha = "left" if label_x > 0.25 else ("right" if label_x < -0.25 else "center")
        va = "bottom" if label_y > 0.25 else ("top" if label_y < -0.25 else "center")
        wrapped = "\n".join(textwrap.wrap(graph.nodes[node]["pathway"], width=LABEL_WRAP_WIDTH))

        ax.text(label_x, label_y, wrapped, ha=ha, va=va, fontsize=LABEL_FONT_SIZE, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.18}, zorder=5)

    # Database arc labels, placed at the circular mean angle of each group.
    for database_name in database_order:
        nodes = grouped_nodes.get(database_name, [])
        if not nodes:
            continue
        angles = np.deg2rad(np.array([node_angles[n] for n in nodes]))
        mean_angle = np.arctan2(np.mean(np.sin(angles)), np.mean(np.cos(angles)))
        label_radius = CIRCLE_RADIUS * 1.63

        ax.text(
            label_radius * np.cos(mean_angle), label_radius * np.sin(mean_angle), database_name, ha="center", va="center",
            fontsize=12, fontweight="bold", color=DATABASE_COLORS[database_name],
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": DATABASE_COLORS[database_name], "linewidth": 1.1, "alpha": 0.95},
            zorder=6,
        )

    ax.set_title("Cross-Species Enrichment Map of Conserved Pathways", fontweight="bold", pad=24)

    database_legend = [
        Line2D([0], [0], marker=DATABASE_MARKERS[d], color="none", markerfacecolor=DATABASE_COLORS[d], markeredgecolor="black", markeredgewidth=1.0, markersize=11, label=d)
        for d in database_order
    ]
    first_legend = ax.legend(handles=database_legend, title="Database", loc="upper left", frameon=True, edgecolor="black", fancybox=False)
    ax.add_artist(first_legend)

    size_legend_values = [np.nanmin(mean_fold_values), np.nanmedian(mean_fold_values), np.nanmax(mean_fold_values)]
    size_legend_sizes = scale_values(np.log2(np.array(size_legend_values) + 1), MIN_NODE_SIZE, MAX_NODE_SIZE)
    size_legend = [ax.scatter([], [], s=s, facecolor="white", edgecolor="black", linewidth=1.0, label=f"{v:.1f}") for v, s in zip(size_legend_values, size_legend_sizes)]
    ax.legend(handles=size_legend, title="Mean fold enrichment", loc="lower left", frameon=True, edgecolor="black", fancybox=False, labelspacing=1.4)

    summary_text = f"Pathways: {graph.number_of_nodes()}\nConnections: {graph.number_of_edges()}"
    ax.text(0.99, 0.01, summary_text, transform=ax.transAxes, ha="right", va="bottom", fontsize=9.5, bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "edgecolor": "black", "linewidth": 0.8, "alpha": 0.92})

    outer_limit = CIRCLE_RADIUS * 1.85
    ax.set_xlim(-outer_limit, outer_limit)
    ax.set_ylim(-outer_limit, outer_limit)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    plt.tight_layout()

    return fig, pos, node_angles


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("results/pathway_conservation/tables"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures/grouped_circular_pathway_map"))
    parser.add_argument("--gene-set-source", choices=["species_a", "species_b", "combined"], default="species_b")
    parser.add_argument("--min-jaccard", type=float, default=0.01)
    parser.add_argument("--min-shared-genes", type=int, default=1)
    parser.add_argument("--remove-isolated-nodes", action="store_true", default=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    apply_publication_style()
    plt.rcParams["font.size"] = 11
    plt.rcParams["axes.titlesize"] = 18

    combined_df = load_conserved_pathways(args.input_dir, SPECIES_A, SPECIES_B)
    database_order = list(DATABASE_COLORS)

    graph, edge_records = build_pathway_network(
        combined_df, SPECIES_A, SPECIES_B, args.gene_set_source,
        args.min_jaccard, args.min_shared_genes, require_both=True, remove_isolates=args.remove_isolated_nodes,
    )

    fig, pos, node_angles = draw_circular_map(graph, database_order)

    node_records = [
        {
            "Node ID": n, "Pathway": graph.nodes[n]["pathway"], "Database": graph.nodes[n]["database"],
            "Mean Fold Enrichment": graph.nodes[n]["mean_fold"], "Gene Count": graph.nodes[n]["gene_count"],
            "Degree": graph.nodes[n]["degree"], "Betweenness Centrality": graph.nodes[n]["betweenness"],
            "Closeness Centrality": graph.nodes[n]["closeness"], "Connected Component": graph.nodes[n]["component"],
            "Circular Angle Degrees": node_angles[n], "Circular X": float(pos[n][0]), "Circular Y": float(pos[n][1]),
            "Genes": "; ".join(sorted(graph.nodes[n]["genes"])),
        }
        for n in graph.nodes
    ]

    node_path = args.output_dir / "pathway_network_nodes_circular.csv"
    edge_path = args.output_dir / "pathway_network_edges_circular.csv"
    graphml_path = args.output_dir / "grouped_circular_pathway_map.graphml"

    pd.DataFrame(node_records).to_csv(node_path, index=False)
    pd.DataFrame(edge_records).to_csv(edge_path, index=False)

    graphml_graph = graph.copy()
    for node in graphml_graph.nodes:
        graphml_graph.nodes[node]["genes"] = "; ".join(sorted(graphml_graph.nodes[node]["genes"]))
        graphml_graph.nodes[node]["x"], graphml_graph.nodes[node]["y"] = float(pos[node][0]), float(pos[node][1])
        graphml_graph.nodes[node]["circular_angle_degrees"] = float(node_angles[node])
    for a, b in graphml_graph.edges:
        graphml_graph.edges[a, b]["shared_genes"] = "; ".join(graphml_graph.edges[a, b]["shared_genes"])
    nx.write_graphml(graphml_graph, graphml_path)

    written = save_figure(fig, args.output_dir / "grouped_circular_pathway_map", formats=("png", "tiff", "pdf", "svg"))
    plt.show()

    print("=" * 85)
    print("CIRCULAR CONSERVED PATHWAY ENRICHMENT MAP COMPLETED")
    print("=" * 85)
    print(f"Nodes: {graph.number_of_nodes()} | Edges: {graph.number_of_edges()}")
    print(f"\nGenerated files:\n" + "\n".join(str(p) for p in written) + f"\n{node_path}\n{edge_path}\n{graphml_path}")


if __name__ == "__main__":
    main()
