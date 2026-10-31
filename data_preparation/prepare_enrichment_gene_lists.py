"""
Prepare gene lists for functional-enrichment tools (e.g. DAVID) and check
their overlap with mapped orthologs.

This script provides three related, independently runnable subcommands:

  overlap
      Intersect a mapped one-to-one ortholog table (see map_orthologs.py)
      with a mutated-gene table from a second species, to find genes with
      evidence of involvement in both species.

  unique-list
      Collapse a per-transcript mutation table into a single deduplicated
      gene list, suitable for submission to DAVID or a similar enrichment
      tool.

  filtered-list
      Same as unique-list, but additionally requires each gene to be
      mutated in at least a minimum number of samples, and reports the
      per-gene mutation frequency.

Usage
-----
    python prepare_enrichment_gene_lists.py overlap \\
        --mapped-orthologs results/ortholog_mapping/one_to_one_human_orthologs.csv \\
        --mutated-genes data/human_mutated_genes.csv \\
        --output results/gene_lists/ortholog_mutation_overlap.csv

    python prepare_enrichment_gene_lists.py unique-list \\
        --mutation-table data/query_gene_mutations.csv \\
        --output-prefix results/gene_lists/unique_genes_for_enrichment

    python prepare_enrichment_gene_lists.py filtered-list \\
        --mutation-table data/query_gene_mutations.csv \\
        --min-mutated-samples 5 \\
        --output-prefix results/gene_lists/filtered_genes_for_enrichment
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

# Transcript-level suffixes (e.g. "A1CF_ENST00000282641") are collapsed to
# the base gene symbol using this pattern.
TRANSCRIPT_SUFFIX_PATTERN = r"_ENST\d+.*$"


def find_ortholog_mutation_overlap(mapped_orthologs: pd.DataFrame, mutated_genes: pd.DataFrame) -> pd.DataFrame:
    """Find genes shared between a mapped-ortholog table and a mutated-gene table.

    ``mapped_orthologs`` must contain a ``Human ortholog gene`` column (as
    produced by ``map_orthologs.py``); ``mutated_genes`` must contain a
    ``Gene name`` column.
    """
    mapped = mapped_orthologs.copy()
    mutated = mutated_genes.copy()

    mapped["Comparison gene"] = mapped["Human ortholog gene"].astype("string").str.strip().str.upper()
    mutated["Comparison gene"] = mutated["Gene name"].astype("string").str.strip().str.upper()

    mapped_clean = mapped[mapped["Comparison gene"].notna() & (mapped["Comparison gene"] != "")]
    mutated_clean = mutated[mutated["Comparison gene"].notna() & (mutated["Comparison gene"] != "")]

    overlap_genes = sorted(set(mapped_clean["Comparison gene"]) & set(mutated_clean["Comparison gene"]))

    overlap_df = mapped_clean[mapped_clean["Comparison gene"].isin(overlap_genes)].copy()
    overlap_df["Found in comparison mutation dataset"] = "Yes"
    overlap_df = overlap_df.drop_duplicates(subset=["Comparison gene"])

    column_order = ["Comparison gene"] + [c for c in overlap_df.columns if c != "Comparison gene"]
    return overlap_df[column_order]


def build_unique_gene_list(mutation_table: pd.DataFrame) -> pd.DataFrame:
    """Collapse a per-transcript mutation table to one row per gene symbol."""
    df = mutation_table.copy()
    df.columns = df.columns.str.strip()

    df["Clean_Gene_Name"] = (
        df["Gene name"].astype(str).str.strip().str.replace(TRANSCRIPT_SUFFIX_PATTERN, "", regex=True)
    )

    unique_genes = df[["Clean_Gene_Name"]].dropna().drop_duplicates().reset_index(drop=True)
    unique_genes.columns = ["Gene name"]
    return unique_genes


def build_filtered_gene_list(mutation_table: pd.DataFrame, min_mutated_samples: int) -> pd.DataFrame:
    """Collapse to gene level, compute mutation frequency, and filter by recurrence."""
    df = mutation_table.copy()
    df.columns = df.columns.str.strip()

    df["Gene name"] = df["Gene name"].astype(str).str.strip()
    df["Mutated samples"] = pd.to_numeric(df["Mutated samples"], errors="coerce")

    df["Gene"] = df["Gene name"].str.replace(TRANSCRIPT_SUFFIX_PATTERN, "", regex=True)
    df = df[df["Gene"].notna() & (df["Gene"] != "") & (~df["Gene"].str.fullmatch(r"\d+"))].copy()

    # A gene's transcripts typically report the same set of mutated samples,
    # so "max" is used to collapse to one row per gene rather than summing.
    gene_summary = df.groupby("Gene", as_index=False).agg(
        Mutated_samples=("Mutated samples", "max"),
        Samples_tested=("Samples tested", "max"),
    )

    gene_summary["Mutation_frequency_percent"] = (
        gene_summary["Mutated_samples"] / gene_summary["Samples_tested"] * 100
    )

    filtered = gene_summary[gene_summary["Mutated_samples"] >= min_mutated_samples].copy()
    filtered = filtered.sort_values(
        by=["Mutated_samples", "Mutation_frequency_percent"], ascending=False
    ).reset_index(drop=True)

    return filtered


def _write_gene_list(genes: pd.Series, output_prefix: Path) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    genes.to_frame(name="Gene name").to_csv(output_prefix.with_suffix(".csv"), index=False)
    genes.to_csv(output_prefix.with_suffix(".txt"), index=False, header=False)


def run_overlap(args: argparse.Namespace) -> None:
    mapped_orthologs = pd.read_csv(args.mapped_orthologs, dtype=str)
    mutated_genes = pd.read_csv(args.mutated_genes, dtype=str)

    overlap_df = find_ortholog_mutation_overlap(mapped_orthologs, mutated_genes)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    overlap_df.to_csv(args.output, index=False, encoding="utf-8-sig")

    print("================ OVERLAP SUMMARY ================")
    print(f"Unique mapped orthologs:         {mapped_orthologs['Human ortholog gene'].nunique()}")
    print(f"Unique comparison-dataset genes: {mutated_genes['Gene name'].nunique()}")
    print(f"Shared genes:                    {len(overlap_df)}")
    print(f"\nOutput saved to:\n{args.output}")


def run_unique_list(args: argparse.Namespace) -> None:
    mutation_table = pd.read_csv(args.mutation_table, sep=None, engine="python")
    unique_genes = build_unique_gene_list(mutation_table)
    _write_gene_list(unique_genes["Gene name"], args.output_prefix)

    print(f"Number of unique genes: {len(unique_genes)}")
    print(f"CSV saved to: {args.output_prefix.with_suffix('.csv')}")
    print(f"TXT saved to: {args.output_prefix.with_suffix('.txt')}")


def run_filtered_list(args: argparse.Namespace) -> None:
    mutation_table = pd.read_csv(args.mutation_table, sep=None, engine="python")
    filtered = build_filtered_gene_list(mutation_table, args.min_mutated_samples)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(args.output_prefix.with_suffix(".csv"), index=False)
    _write_gene_list(filtered["Gene"], args.output_prefix)

    print(f"Genes mutated in at least {args.min_mutated_samples} sample(s): {len(filtered)}")
    print(f"Gene list saved to: {args.output_prefix.with_suffix('.txt')}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest="command", required=True)

    overlap_parser = subparsers.add_parser("overlap", help="Intersect mapped orthologs with a mutated-gene list.")
    overlap_parser.add_argument("--mapped-orthologs", type=Path, required=True)
    overlap_parser.add_argument("--mutated-genes", type=Path, required=True)
    overlap_parser.add_argument("--output", type=Path, required=True)
    overlap_parser.set_defaults(func=run_overlap)

    unique_parser = subparsers.add_parser("unique-list", help="Export a deduplicated gene list for enrichment tools.")
    unique_parser.add_argument("--mutation-table", type=Path, required=True)
    unique_parser.add_argument("--output-prefix", type=Path, required=True)
    unique_parser.set_defaults(func=run_unique_list)

    filtered_parser = subparsers.add_parser(
        "filtered-list", help="Export a recurrence-filtered gene list for enrichment tools."
    )
    filtered_parser.add_argument("--mutation-table", type=Path, required=True)
    filtered_parser.add_argument("--min-mutated-samples", type=int, default=5)
    filtered_parser.add_argument("--output-prefix", type=Path, required=True)
    filtered_parser.set_defaults(func=run_filtered_list)

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
