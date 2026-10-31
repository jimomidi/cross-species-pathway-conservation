"""
Test whether the pathway overlap identified by ``compare_enrichment_results.py``
is greater than expected by chance, using a real pathway "universe" size (N)
retrieved from public databases rather than an arbitrary assumption.

Two alternative universe definitions are provided as subcommands, since
they answer slightly different questions and were both explored during the
analysis this toolkit is based on:

  combined-universe
      N = the number of pathways annotated for EITHER species (KEGG and
      Reactome only; a species-specific union is not available for GO
      Biological Process without downloading a full annotation file, so
      GO BP is skipped in this mode). Tests with Fisher's exact test.

  reference-universe
      N = the number of pathways annotated for the reference species only
      (GO BP, KEGG, and Reactome). This framing is appropriate when the
      query species' genes were mapped to reference-species orthologs and
      both gene lists were analyzed against reference-species annotation,
      so a single shared background applies. Tests with both Fisher's
      exact test and the hypergeometric test. GO Biological Process
      requires a local, pre-downloaded GOA GAF file for the reference
      species (see --goa-file).

Both subcommands read the combined summary table produced by
``compare_enrichment_results.py`` (``Combined_Pathway_Conservation_Summary.csv``)
for the observed shared / species-specific pathway counts.

Usage
-----
    python pathway_universe_enrichment_tests.py combined-universe \\
        --summary results/pathway_conservation/tables/Combined_Pathway_Conservation_Summary.csv \\
        --output results/pathway_conservation/tables/combined_universe_enrichment_tests.csv

    python pathway_universe_enrichment_tests.py reference-universe \\
        --summary results/pathway_conservation/tables/Combined_Pathway_Conservation_Summary.csv \\
        --goa-file data/reference_annotation/goa_human.gaf.gz \\
        --output results/pathway_conservation/tables/reference_universe_enrichment_tests.csv
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd
import requests
from scipy.stats import fisher_exact, hypergeom

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import KEGG_ORGANISM_CODE, NCBI_TAXONOMY_ID, SPECIES_A, SPECIES_B  # noqa: E402

DATABASE_DISPLAY_NAMES = {"GO_BP": "GO Biological Process", "KEGG": "KEGG", "REACTOME": "Reactome"}


# ---------------------------------------------------------------------------
# Universe retrieval
# ---------------------------------------------------------------------------
def get_kegg_pathway_set(organism_code: str) -> set[str]:
    """Return KEGG pathway IDs currently annotated for a KEGG organism code."""
    response = requests.get(f"https://rest.kegg.jp/list/pathway/{organism_code}", timeout=30)
    response.raise_for_status()
    return {line.split("\t")[0].replace("path:", "") for line in response.text.strip().splitlines() if line}


def _collect_reactome_pathway_ids(node, pathway_ids: set[str]) -> None:
    if isinstance(node, dict):
        if node.get("type") == "Pathway" and "stId" in node:
            pathway_ids.add(node["stId"])
        for child in node.get("children", []):
            _collect_reactome_pathway_ids(child, pathway_ids)
    elif isinstance(node, list):
        for item in node:
            _collect_reactome_pathway_ids(item, pathway_ids)


def get_reactome_pathway_set(taxon_id: str) -> set[str]:
    """Return curated Reactome pathway stable IDs for an NCBI taxonomy ID."""
    response = requests.get(
        f"https://reactome.org/ContentService/data/eventsHierarchy/{taxon_id}",
        headers={"Accept": "application/json"},
        timeout=60,
    )
    response.raise_for_status()
    pathway_ids: set[str] = set()
    _collect_reactome_pathway_ids(response.json(), pathway_ids)
    return pathway_ids


def get_go_bp_terms_from_gaf(gaf_gz_path: Path) -> set[str]:
    """Read a gzip-compressed GOA GAF file and return unique Biological Process GO IDs.

    GAF column 5 is the GO ID and column 9 is the aspect ("P" = Biological
    Process, "F" = Molecular Function, "C" = Cellular Component).
    """
    if not gaf_gz_path.exists():
        raise FileNotFoundError(f"GOA annotation file was not found:\n{gaf_gz_path}")

    go_bp_terms: set[str] = set()
    with gzip.open(gaf_gz_path, mode="rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[8].strip() == "P" and fields[4].strip():
                go_bp_terms.add(fields[4].strip())

    return go_bp_terms


# ---------------------------------------------------------------------------
# Summary-table helpers
# ---------------------------------------------------------------------------
def get_summary_value(summary_row: pd.Series, possible_columns: list[str]) -> int:
    """Return the first matching column's value, tolerating minor naming differences."""
    for column in possible_columns:
        if column in summary_row.index:
            return int(summary_row[column])
    raise KeyError(f"None of these columns were found: {possible_columns}\nAvailable: {list(summary_row.index)}")


def _database_row(summary: pd.DataFrame, database_display_name: str) -> pd.Series:
    matches = summary[summary["Database"] == database_display_name]
    if matches.empty:
        raise KeyError(f"No row for database '{database_display_name}' in the summary table.")
    return matches.iloc[0]


# ---------------------------------------------------------------------------
# Combined-universe mode (KEGG, Reactome)
# ---------------------------------------------------------------------------
def run_combined_universe(summary: pd.DataFrame) -> pd.DataFrame:
    results = []

    for database_name in ["KEGG", "REACTOME"]:
        display_name = DATABASE_DISPLAY_NAMES[database_name]
        row = _database_row(summary, display_name)

        a = get_summary_value(row, ["Shared Conserved Pathways", "Shared Conserved Terms"])
        b = get_summary_value(row, [f"{SPECIES_A}-Specific Pathways", f"{SPECIES_A}-Specific Terms", "Dog-Specific Terms", "Dog-Specific Pathways"])
        c = get_summary_value(row, [f"{SPECIES_B}-Specific Pathways", f"{SPECIES_B}-Specific Terms", "Human-Specific Terms", "Human-Specific Pathways"])

        if database_name == "KEGG":
            organism_set_a = get_kegg_pathway_set(KEGG_ORGANISM_CODE[SPECIES_A])
            organism_set_b = get_kegg_pathway_set(KEGG_ORGANISM_CODE[SPECIES_B])
        else:
            organism_set_a = get_reactome_pathway_set(NCBI_TAXONOMY_ID[SPECIES_A])
            organism_set_b = get_reactome_pathway_set(NCBI_TAXONOMY_ID[SPECIES_B])

        universe_n = len(organism_set_a | organism_set_b)
        d = universe_n - a - b - c

        if d < 0:
            print(f"{display_name}: universe (N={universe_n}) smaller than a+b+c; skipping.")
            continue

        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

        results.append(
            {
                "Database": display_name,
                "Universe Definition": f"Union of {SPECIES_A} and {SPECIES_B} annotated pathways",
                "Universe Size (N)": universe_n,
                "Shared (a)": a,
                f"{SPECIES_A}-only (b)": b,
                f"{SPECIES_B}-only (c)": c,
                "Neither (d)": d,
                "Odds Ratio": odds_ratio,
                "Fisher Exact P-Value": p_value,
            }
        )

    print(
        f"\nNote: GO Biological Process is intentionally excluded from the combined-universe test. "
        f"A species-specific GO BP universe requires the full annotation set for each species and "
        f"is not available through a simple pathway-listing API; see the reference-universe subcommand."
    )

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Reference-universe mode (GO BP, KEGG, Reactome)
# ---------------------------------------------------------------------------
def run_reference_universe(summary: pd.DataFrame, goa_file: Path | None) -> pd.DataFrame:
    universes: dict[str, int] = {
        "KEGG": len(get_kegg_pathway_set(KEGG_ORGANISM_CODE[SPECIES_B])),
        "REACTOME": len(get_reactome_pathway_set(NCBI_TAXONOMY_ID[SPECIES_B])),
    }

    if goa_file is not None:
        universes["GO_BP"] = len(get_go_bp_terms_from_gaf(goa_file))
    else:
        print("No --goa-file provided; GO Biological Process will be skipped in this mode.")

    print("Reference-species annotation universes:")
    for database_name, size in universes.items():
        print(f"  {DATABASE_DISPLAY_NAMES[database_name]}: {size}")

    results = []
    for database_name, universe_n in universes.items():
        display_name = DATABASE_DISPLAY_NAMES[database_name]
        row = _database_row(summary, display_name)

        a_total = get_summary_value(row, [f"{SPECIES_A} Retained Pathways", f"{SPECIES_A} Retained Terms"])
        b_total = get_summary_value(row, [f"{SPECIES_B} Retained Pathways", f"{SPECIES_B} Retained Terms"])
        shared = get_summary_value(row, ["Shared Conserved Pathways", "Shared Conserved Terms"])

        a_only = a_total - shared
        b_only = b_total - shared
        neither = universe_n - shared - a_only - b_only

        if neither < 0:
            raise ValueError(
                f"{display_name}: reference universe N={universe_n} is smaller than the observed "
                f"union ({shared + a_only + b_only}). This may indicate a mismatch between the "
                f"DAVID annotation version and the reference universe retrieval date."
            )

        odds_ratio, fisher_p = fisher_exact([[shared, a_only], [b_only, neither]], alternative="greater")
        hypergeom_p = hypergeom.sf(shared - 1, universe_n, b_total, a_total)

        results.append(
            {
                "Database": display_name,
                f"{SPECIES_B} Annotation Universe": universe_n,
                f"{SPECIES_A} Pathways": a_total,
                f"{SPECIES_B} Pathways": b_total,
                "Shared Conserved Pathways": shared,
                f"{SPECIES_A}-Specific Pathways": a_only,
                f"{SPECIES_B}-Specific Pathways": b_only,
                "Neither": neither,
                f"{SPECIES_A} Pathway Conservation (%)": 100 * shared / a_total if a_total else float("nan"),
                "Odds Ratio": odds_ratio,
                "Fisher Exact P-Value": fisher_p,
                "Hypergeometric P-Value": hypergeom_p,
            }
        )

    return pd.DataFrame(results)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest="command", required=True)

    combined_parser = subparsers.add_parser("combined-universe", help="Union-of-both-species universe (KEGG, Reactome).")
    combined_parser.add_argument("--summary", type=Path, required=True)
    combined_parser.add_argument("--output", type=Path, required=True)

    reference_parser = subparsers.add_parser("reference-universe", help="Reference-species-only universe (GO BP, KEGG, Reactome).")
    reference_parser.add_argument("--summary", type=Path, required=True)
    reference_parser.add_argument("--goa-file", type=Path, default=None, help="Local gzip-compressed GOA GAF file for the reference species.")
    reference_parser.add_argument("--output", type=Path, required=True)

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = pd.read_csv(args.summary)

    if args.command == "combined-universe":
        results = run_combined_universe(summary)
    else:
        results = run_reference_universe(summary, args.goa_file)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output, index=False)

    print(f"\n{results.to_string(index=False)}")
    print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()
