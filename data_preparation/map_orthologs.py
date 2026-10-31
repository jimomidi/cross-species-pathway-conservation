"""
Map a query-species gene list to human orthologs using an Ensembl BioMart
export, and summarize the resulting mapping.

This script expects two inputs:

1. A single-column gene list for the query species (e.g. genes carrying a
   mutation of interest), produced upstream of this pipeline.
2. An Ensembl BioMart "Homologues" export (converted to CSV by
   ``convert_biomart_export.py``) containing, at minimum, the columns:
   ``Gene name``, ``Gene stable ID``, ``Human gene stable ID``,
   ``Human gene name``, ``Human homology type``,
   ``%id. target Human gene identical to query gene``,
   ``%id. query gene identical to target Human gene``, and
   ``Human orthology confidence [0 low, 1 high]``.
   These are the standard column names produced by BioMart when exporting
   human-orthology attributes, so they are kept as-is regardless of the
   query species label configured in ``config.py``.

For each query gene, the script reports whether a human ortholog was found,
the orthology type (one-to-one, one-to-many, many-to-many) and confidence,
and sequence-identity statistics. It also writes a one-to-one-only subset,
which is typically the most appropriate set for downstream pathway
enrichment comparisons.

Usage
-----
    python map_orthologs.py \\
        --gene-list data/query_gene_list.csv \\
        --mapping data/biomart/mart_export.csv \\
        --output-dir results/ortholog_mapping
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import SPECIES_A  # noqa: E402

MAPPING_COLUMNS = [
    "Query gene key",
    "Gene name",
    "Gene stable ID",
    "Human gene stable ID",
    "Human gene name",
    "Human homology type",
    "%id. target Human gene identical to query gene",
    "%id. query gene identical to target Human gene",
    "Human orthology confidence [0 low, 1 high]",
]

MAPPING_RENAME = {
    "Gene name": "Query gene name",
    "Gene stable ID": "Query Ensembl gene ID",
    "Human gene stable ID": "Human Ensembl gene ID",
    "Human gene name": "Human ortholog gene",
    "Human homology type": "Orthology type",
    "%id. target Human gene identical to query gene": "Human sequence identity to query (%)",
    "%id. query gene identical to target Human gene": "Query sequence identity to human (%)",
    "Human orthology confidence [0 low, 1 high]": "Orthology confidence",
}

RESULT_COLUMN_ORDER = [
    "Original query gene",
    "Query gene name",
    "Query Ensembl gene ID",
    "Human ortholog gene",
    "Human Ensembl gene ID",
    "Orthology type",
    "Orthology confidence",
    "Human sequence identity to query (%)",
    "Query sequence identity to human (%)",
    "Mapping status",
]


def load_gene_list(gene_list_path: Path) -> pd.DataFrame:
    """Read a single-column gene list and standardize it for matching."""
    genes = pd.read_csv(gene_list_path, dtype=str)

    if genes.shape[1] != 1:
        raise ValueError(
            f"{gene_list_path.name} should contain exactly one column, "
            f"but {genes.shape[1]} columns were found."
        )

    genes = genes.rename(columns={genes.columns[0]: "Original query gene"})
    genes["Query gene key"] = genes["Original query gene"].astype(str).str.strip().str.upper()

    return genes


def load_mapping_table(mapping_path: Path) -> pd.DataFrame:
    """Read and clean a BioMart export prepared by convert_biomart_export.py."""
    mapping = pd.read_csv(mapping_path, dtype=str)

    missing = [c for c in MAPPING_RENAME if c not in mapping.columns]
    if missing:
        raise KeyError(
            f"The following expected BioMart columns are missing from "
            f"{mapping_path.name}: {missing}\n\n"
            f"Available columns: {mapping.columns.tolist()}"
        )

    mapping = mapping[
        mapping["Gene name"].notna() & (mapping["Gene name"].str.strip() != "")
    ].copy()

    mapping["Query gene key"] = mapping["Gene name"].astype(str).str.strip().str.upper()
    mapping = mapping[MAPPING_COLUMNS].drop_duplicates()
    mapping = mapping.rename(columns=MAPPING_RENAME)

    return mapping


def map_to_human_orthologs(genes: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    """Left-join the query gene list to the ortholog mapping table."""
    result = genes.merge(mapping, on="Query gene key", how="left")

    result["Mapping status"] = result["Human ortholog gene"].apply(
        lambda value: "Mapped" if pd.notna(value) and str(value).strip() != "" else "Unmapped"
    )

    result = result[RESULT_COLUMN_ORDER].drop_duplicates()

    return result


def summarize_mapping(result: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Compute mapping-rate, orthology-type, confidence, and identity summaries."""
    mapped = result[result["Mapping status"] == "Mapped"].copy()
    unmapped = result[result["Mapping status"] == "Unmapped"].copy()

    total_genes = result["Original query gene"].nunique()
    mapped_genes = mapped["Original query gene"].nunique()
    unmapped_genes = unmapped["Original query gene"].nunique()
    mapping_rate = (mapped_genes / total_genes * 100) if total_genes else 0.0

    general_summary = pd.DataFrame(
        [
            {"Metric": "Total unique query genes", "Value": total_genes},
            {"Metric": "Mapped unique query genes", "Value": mapped_genes},
            {"Metric": "Unmapped unique query genes", "Value": unmapped_genes},
            {"Metric": "Mapping rate (%)", "Value": round(mapping_rate, 2)},
            {"Metric": "Unique human ortholog genes", "Value": mapped["Human ortholog gene"].nunique()},
        ]
    )

    orthology_rows = []
    for orthology_type in ["ortholog_one2one", "ortholog_one2many", "ortholog_many2many"]:
        subset = mapped[mapped["Orthology type"] == orthology_type]
        unique_query = subset["Original query gene"].nunique()
        orthology_rows.append(
            {
                "Category": orthology_type,
                "Unique query genes": unique_query,
                "Unique human orthologs": subset["Human ortholog gene"].nunique(),
                "Mapping rows": len(subset),
                "Percent of mapped query genes": round(
                    unique_query / mapped_genes * 100 if mapped_genes else 0.0, 2
                ),
            }
        )
    orthology_summary = pd.DataFrame(orthology_rows)

    for column in ["Human sequence identity to query (%)", "Query sequence identity to human (%)", "Orthology confidence"]:
        mapped[column] = pd.to_numeric(mapped[column], errors="coerce")

    confidence_summary = (
        mapped.dropna(subset=["Orthology confidence"])
        .groupby("Orthology confidence")
        .agg(
            Unique_query_genes=("Original query gene", "nunique"),
            Unique_human_genes=("Human ortholog gene", "nunique"),
            Mapping_rows=("Original query gene", "size"),
        )
        .reset_index()
    )
    confidence_summary["Confidence category"] = (
        confidence_summary["Orthology confidence"].map({0: "Low confidence", 1: "High confidence"}).fillna("Other")
    )
    confidence_summary["Percent of mapped query genes"] = (
        confidence_summary["Unique_query_genes"] / mapped_genes * 100 if mapped_genes else 0.0
    ).round(2) if mapped_genes else confidence_summary["Unique_query_genes"]

    identity_records = []
    for group_name, group_df in [
        ("All mapped orthologs", mapped),
        ("One-to-one orthologs", mapped[mapped["Orthology type"] == "ortholog_one2one"]),
    ]:
        for column in ["Human sequence identity to query (%)", "Query sequence identity to human (%)"]:
            values = group_df[column].dropna()
            identity_records.append(
                {
                    "Mapping group": group_name,
                    "Identity measure": column,
                    "N": values.count(),
                    "Mean (%)": round(values.mean(), 2) if len(values) else None,
                    "Median (%)": round(values.median(), 2) if len(values) else None,
                    "Standard deviation": round(values.std(), 2) if len(values) > 1 else None,
                    "Minimum (%)": round(values.min(), 2) if len(values) else None,
                    "Maximum (%)": round(values.max(), 2) if len(values) else None,
                    "Q1 (%)": round(values.quantile(0.25), 2) if len(values) else None,
                    "Q3 (%)": round(values.quantile(0.75), 2) if len(values) else None,
                }
            )
    identity_summary = pd.DataFrame(identity_records)

    changed_names = mapped[
        mapped["Original query gene"].notna()
        & mapped["Human ortholog gene"].notna()
        & (mapped["Original query gene"].str.upper() != mapped["Human ortholog gene"].str.upper())
    ][
        [
            "Original query gene",
            "Human ortholog gene",
            "Query Ensembl gene ID",
            "Human Ensembl gene ID",
            "Orthology type",
            "Orthology confidence",
            "Human sequence identity to query (%)",
            "Query sequence identity to human (%)",
        ]
    ].drop_duplicates()

    unmapped_supplementary = (
        unmapped[["Original query gene", "Query Ensembl gene ID"]]
        .drop_duplicates(subset=["Original query gene"])
        .sort_values("Original query gene")
    )

    one_to_one = (
        mapped[mapped["Orthology type"] == "ortholog_one2one"]
        .drop_duplicates(subset=["Original query gene", "Human ortholog gene"])
    )

    return {
        "general_summary": general_summary,
        "orthology_summary": orthology_summary,
        "confidence_summary": confidence_summary,
        "identity_summary": identity_summary,
        "changed_names": changed_names,
        "unmapped_supplementary": unmapped_supplementary,
        "one_to_one": one_to_one,
        "mapped": mapped,
        "unmapped": unmapped,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gene-list",
        type=Path,
        default=Path("data/query_gene_list.csv"),
        help=f"Single-column CSV of {SPECIES_A} genes to map (default: %(default)s).",
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=Path("data/biomart/mart_export.csv"),
        help="BioMart export CSV, see convert_biomart_export.py (default: %(default)s).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/ortholog_mapping"),
        help="Directory for output tables (default: %(default)s).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    genes = load_gene_list(args.gene_list)
    mapping = load_mapping_table(args.mapping)

    result = map_to_human_orthologs(genes, mapping)
    result.to_csv(output_dir / "gene_list_human_orthologs.csv", index=False, encoding="utf-8-sig")

    summaries = summarize_mapping(result)

    summaries["one_to_one"].to_csv(output_dir / "one_to_one_human_orthologs.csv", index=False, encoding="utf-8-sig")
    summaries["unmapped_supplementary"].to_csv(output_dir / "unmapped_query_genes.csv", index=False, encoding="utf-8-sig")
    summaries["changed_names"].to_csv(output_dir / "query_human_changed_gene_names.csv", index=False, encoding="utf-8-sig")

    combined_summary = pd.concat(
        [
            summaries["general_summary"],
            pd.DataFrame([{"Metric": "", "Value": ""}]),
            summaries["orthology_summary"].rename(columns={"Category": "Metric", "Unique query genes": "Value"})[
                ["Metric", "Value"]
            ],
        ],
        ignore_index=True,
    )
    combined_summary.to_csv(output_dir / "ortholog_mapping_summary.csv", index=False, encoding="utf-8-sig")
    summaries["confidence_summary"].to_csv(output_dir / "orthology_confidence_summary.csv", index=False, encoding="utf-8-sig")
    summaries["identity_summary"].to_csv(output_dir / "ortholog_identity_statistics.csv", index=False, encoding="utf-8-sig")

    print("================ ORTHOLOG MAPPING REPORT ================")
    print(combined_summary.to_string(index=False))
    print("\nIdentity statistics:")
    print(summaries["identity_summary"].to_string(index=False))
    print(f"\nGenes with different {SPECIES_A}/human symbols: {summaries['changed_names']['Original query gene'].nunique()}")
    print(f"\nOutput files written to:\n{output_dir}")


if __name__ == "__main__":
    main()
