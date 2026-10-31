"""
Network map of conserved pathways, connected by shared-gene overlap
(Jaccard similarity of the gene sets underlying each pathway).

Nodes are pathways (colored/shaped by source database, sized by mean fold
enrichment); edges connect pathway pairs that share at least
``--min-shared-genes`` genes with Jaccard similarity at or above
``--min-jaccard`` (default 0.01, with isolated nodes removed by default).
An optional pathway allowlist (``--pathway-allowlist``, one pathway name
per line) restricts the plot to a curated subset -- for example, when
building a smaller "main text" figure from a larger supplementary network
(a higher ``--min-jaccard`` such as 0.08 is typically paired with an
allowlist for that kind of curated figure).

Input: the three "conserved shared pathways" CSVs produced by
``pathway_enrichment/compare_enrichment_results.py`` (must include
``<species> User Ids`` columns, i.e. the DAVID gene-identifier list per
term).

Usage
-----
    python pathway_network_map.py \\
        --input-dir results/pathway_conservation/tables \\
        --output-dir results/figures/pathway_network_map

    python pathway_network_map.py \\
        --input-dir results/pathway_conservation/tables \\
        --output-dir results/figures/pathway_network_map_main_text \\
        --pathway-allowlist data/main_text_pathways.txt \\
        --min-jaccard 0.08
"""

from __future__ import annotations

import argparse
import sys
import textwrap
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import DATABASE_COLORS, DATABASE_MARKERS, SPECIES_A, SPECIES_B  # noqa: E402
from plotting_utils import apply_publication_style, load_conserved_pathways, parse_gene_ids, save_figure, scale_values  # noqa: E402

# ---------------------------------------------------------------------------
# Layout / appearance defaults
# ---------------------------------------------------------------------------
LAYOUT_SEED = 42
LAYOUT_K = 1.8
LAYOUT_ITERATIONS = 1000
LABEL_WRAP_WIDTH = 30
EDGE_COLOR = "#8A8A8A"
MIN_NODE_SIZE, MAX_NODE_SIZE = 850, 3300
MIN_EDGE_WIDTH, MAX_EDGE_WIDTH = 0.8, 5.5


def build_pathway_network(
    combined_df,
    species_a: str,
    species_b: str,
    gene_set_source: str,
    min_jaccard: float,
    min_shared_genes: int,
    require_both: bool,
    remove_isolates: bool,
) -> tuple[nx.Graph, list[dict]]:
    """Build a gene-overlap network among pathways.

    ``gene_set_source`` selects which species' "User Ids" column defines
    each pathway's gene set: "species_a", "species_b", or "combined"
    (union of both). Using the reference/ortholog-annotated species
    (typically species_b/human) is recommended, since all pathways are
    then expressed in the same gene-identifier namespace.
    """
    a_genes_col, b_genes_col = f"{species_a} User Ids", f"{species_b} User Ids"
    combined_df["A Gene Set"] = combined_df[a_genes_col].apply(parse_gene_ids)
    combined_df["B Gene Set"] = combined_df[b_genes_col].apply(parse_gene_ids)

    if gene_set_source == "species_a":
        combined_df["Network Gene Set"] = combined_df["A Gene Set"]
    elif gene_set_source == "species_b":
        combined_df["Network Gene Set"] = combined_df["B Gene Set"]
    elif gene_set_source == "combined":
        combined_df["Network Gene Set"] = combined_df.apply(lambda r: r["A Gene Set"] | r["B Gene Set"], axis=1)
    else:
        raise ValueError("gene_set_source must be 'species_a', 'species_b', or 'combined'.")

    combined_df = combined_df[combined_df["Network Gene Set"].apply(len).gt(0)].copy()
    combined_df = combined_df.drop_duplicates(subset=["Database", "Pathway"]).reset_index(drop=True)
    combined_df["Node ID"] = [f"P{i + 1}" for i in range(len(combined_df))]

    graph = nx.Graph()
    for _, row in combined_df.iterrows():
        graph.add_node(
            row["Node ID"], pathway=row["Pathway"], database=row["Database"],
            mean_fold=row["Mean Fold Enrichment"], gene_count=len(row["Network Gene Set"]), genes=row["Network Gene Set"],
        )

    edge_records = []
    for (_, row_a), (_, row_b) in combinations(combined_df.iterrows(), 2):
        intersection = row_a["Network Gene Set"] & row_b["Network Gene Set"]
        union = row_a["Network Gene Set"] | row_b["Network Gene Set"]
        shared_count = len(intersection)
        jaccard = shared_count / len(union) if union else 0.0

        retain = (jaccard >= min_jaccard and shared_count >= min_shared_genes) if require_both else (
            jaccard >= min_jaccard or shared_count >= min_shared_genes
        )

        if retain:
            graph.add_edge(row_a["Node ID"], row_b["Node ID"], jaccard=jaccard, shared_gene_count=shared_count, shared_genes=sorted(intersection))
            edge_records.append(
                {
                    "Node 1": row_a["Node ID"], "Pathway 1": row_a["Pathway"], "Database 1": row_a["Database"],
                    "Node 2": row_b["Node ID"], "Pathway 2": row_b["Pathway"], "Database 2": row_b["Database"],
                    "Shared Gene Count": shared_count, "Jaccard Similarity": jaccard, "Shared Genes": "; ".join(sorted(intersection)),
                }
            )

    if remove_isolates:
        graph.remove_nodes_from(list(nx.isolates(graph)))

    if graph.number_of_nodes() == 0:
        raise ValueError("The graph has no nodes after filtering; relax the thresholds.")

    for node, component_id in ((n, i) for i, comp in enumerate(nx.connected_components(graph), start=1) for n in comp):
        graph.nodes[node]["component"] = component_id

    degree = dict(graph.degree())
    nx.set_node_attributes(graph, degree, "degree")
    nx.set_node_attributes(graph, nx.betweenness_centrality(graph), "betweenness")
    nx.set_node_attributes(graph, nx.closeness_centrality(graph), "closeness")

    return graph, edge_records


def draw_network(graph: nx.Graph, title: str):
    mean_fold_values = np.array([graph.nodes[n]["mean_fold"] for n in graph.nodes])
    node_sizes = dict(zip(graph.nodes, scale_values(np.log2(mean_fold_values + 1), MIN_NODE_SIZE, MAX_NODE_SIZE)))

    edges = list(graph.edges)
    edge_widths = scale_values([graph.edges[e]["jaccard"] for e in edges], MIN_EDGE_WIDTH, MAX_EDGE_WIDTH) if edges else []

    pos = nx.spring_layout(graph, seed=LAYOUT_SEED, k=LAYOUT_K, iterations=LAYOUT_ITERATIONS, weight="jaccard") if edges else nx.circular_layout(graph)

    fig, ax = plt.subplots(figsize=(16, 13))

    if edges:
        nx.draw_networkx_edges(graph, pos, ax=ax, edgelist=edges, width=edge_widths, edge_color=EDGE_COLOR, alpha=0.48)

    for database_name in DATABASE_COLORS:
        nodes = [n for n in graph.nodes if graph.nodes[n]["database"] == database_name]
        if not nodes:
            continue
        nx.draw_networkx_nodes(
            graph, pos, ax=ax, nodelist=nodes, node_size=[node_sizes[n] for n in nodes],
            node_color=DATABASE_COLORS[database_name], node_shape=DATABASE_MARKERS[database_name],
            edgecolors="black", linewidths=1.4, alpha=0.92, label=database_name,
        )

    labels = {n: "\n".join(textwrap.wrap(graph.nodes[n]["pathway"], width=LABEL_WRAP_WIDTH)) for n in graph.nodes}
    nx.draw_networkx_labels(graph, pos, labels=labels, ax=ax, font_size=9.2, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.20})

    ax.set_title(title, fontweight="bold", pad=20)

    database_legend = [
        Line2D([0], [0], marker=DATABASE_MARKERS[d], color="none", markerfacecolor=DATABASE_COLORS[d], markeredgecolor="black", markersize=11, label=d)
        for d in DATABASE_COLORS
    ]
    first_legend = ax.legend(handles=database_legend, title="Database", loc="upper left", frameon=True, edgecolor="black", fancybox=False)
    ax.add_artist(first_legend)

    size_legend_values = [np.nanmin(mean_fold_values), np.nanmedian(mean_fold_values), np.nanmax(mean_fold_values)]
    size_legend_sizes = scale_values(np.log2(np.array(size_legend_values) + 1), MIN_NODE_SIZE, MAX_NODE_SIZE)
    size_legend = [
        plt.scatter([], [], s=size, facecolor="white", edgecolor="black", linewidth=1.1, label=f"{value:.1f}")
        for value, size in zip(size_legend_values, size_legend_sizes)
    ]
    ax.legend(handles=size_legend, title="Mean fold enrichment", loc="lower left", frameon=True, edgecolor="black", fancybox=False, labelspacing=1.5)

    summary_text = f"Pathways: {graph.number_of_nodes()}\nConnections: {graph.number_of_edges()}"
    ax.text(0.99, 0.01, summary_text, transform=ax.transAxes, ha="right", va="bottom", fontsize=10, bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "edgecolor": "black", "linewidth": 0.8, "alpha": 0.90})

    ax.axis("off")
    plt.tight_layout()
    return fig, pos


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", type=Path, default=Path("results/pathway_conservation/tables"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures/pathway_network_map"))
    parser.add_argument("--gene-set-source", choices=["species_a", "species_b", "combined"], default="species_b")
    parser.add_argument("--min-jaccard", type=float, default=0.01)
    parser.add_argument("--min-shared-genes", type=int, default=1)
    parser.add_argument("--require-both-criteria", action="store_true", default=True)
    parser.add_argument("--remove-isolated-nodes", action="store_true", default=True)
    parser.add_argument(
        "--pathway-allowlist", type=Path, default=None,
        help="Optional text file, one pathway name per line, to restrict the network to a curated subset.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    apply_publication_style()
    plt.rcParams["font.size"] = 11
    plt.rcParams["axes.titlesize"] = 18

    combined_df = load_conserved_pathways(args.input_dir, SPECIES_A, SPECIES_B)

    if args.pathway_allowlist is not None:
        allowlist = {line.strip().lower() for line in args.pathway_allowlist.read_text().splitlines() if line.strip()}
        combined_df["_normalized"] = combined_df["Pathway"].str.strip().str.lower()
        missing = allowlist - set(combined_df["_normalized"])
        if missing:
            print("Requested pathways not found in the input tables:")
            for pathway in sorted(missing):
                print(f"  - {pathway}")
        combined_df = combined_df[combined_df["_normalized"].isin(allowlist)].drop(columns="_normalized")
        if combined_df.empty:
            raise ValueError("None of the allowlisted pathway names matched the input tables.")

    graph, edge_records = build_pathway_network(
        combined_df, SPECIES_A, SPECIES_B, args.gene_set_source,
        args.min_jaccard, args.min_shared_genes, args.require_both_criteria, args.remove_isolated_nodes,
    )

    title = "Conserved Pathway Enrichment Map" if args.pathway_allowlist is None else "Conserved Pathway Enrichment Map (Curated Subset)"
    fig, _ = draw_network(graph, title)

    node_df_records = [
        {
            "Node ID": n, "Pathway": graph.nodes[n]["pathway"], "Database": graph.nodes[n]["database"],
            "Mean Fold Enrichment": graph.nodes[n]["mean_fold"], "Gene Count": graph.nodes[n]["gene_count"],
            "Degree": graph.nodes[n]["degree"], "Betweenness Centrality": graph.nodes[n]["betweenness"],
            "Closeness Centrality": graph.nodes[n]["closeness"], "Connected Component": graph.nodes[n]["component"],
            "Genes": "; ".join(sorted(graph.nodes[n]["genes"])),
        }
        for n in graph.nodes
    ]

    node_path = args.output_dir / "pathway_network_nodes.csv"
    edge_path = args.output_dir / "pathway_network_edges.csv"
    graphml_path = args.output_dir / "pathway_network_map.graphml"

    pd.DataFrame(node_df_records).to_csv(node_path, index=False)
    pd.DataFrame(edge_records).to_csv(edge_path, index=False)

    graphml_graph = graph.copy()
    for node in graphml_graph.nodes:
        graphml_graph.nodes[node]["genes"] = "; ".join(sorted(graphml_graph.nodes[node]["genes"]))
    for node_a, node_b in graphml_graph.edges:
        graphml_graph.edges[node_a, node_b]["shared_genes"] = "; ".join(graphml_graph.edges[node_a, node_b]["shared_genes"])
    nx.write_graphml(graphml_graph, graphml_path)

    written = save_figure(fig, args.output_dir / "pathway_network_map", formats=("png", "tiff", "pdf", "svg"))
    plt.show()

    print("=" * 80)
    print("PATHWAY NETWORK MAP COMPLETED")
    print("=" * 80)
    print(f"Nodes: {graph.number_of_nodes()} | Edges: {graph.number_of_edges()} | Connected components: {nx.number_connected_components(graph)}")
    print(f"\nGenerated files:\n" + "\n".join(str(p) for p in written) + f"\n{node_path}\n{edge_path}\n{graphml_path}")


if __name__ == "__main__":
    main()
