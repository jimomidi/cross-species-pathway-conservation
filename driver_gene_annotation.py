"""
Cross-reference a candidate gene list against two curated cancer-driver-gene
resources: the COSMIC Cancer Gene Census (CGC) and IntOGen.

For each gene in the input list, this script reports whether it is
catalogued in the COSMIC CGC and/or identified as a driver by IntOGen,
combines both sources into a consensus classification, and (optionally)
flags CGC entries associated with a configurable set of tissue/disease
keywords (see ``TISSUE_KEYWORDS`` below) as an annotation aid only -- it is
not evidence that a gene is specific to any particular condition.

Expected inputs
----------------
--genes
    A CSV with a single ``Gene`` column of gene symbols to annotate.
--cosmic-cgc
    The COSMIC Cancer Gene Census export (tab-separated), with at least the
    columns: GENE_SYMBOL, TIER, COSMIC_GENE_ID, ROLE_IN_CANCER,
    MUTATION_TYPES, TUMOUR_TYPES_SOMATIC, TISSUE_TYPE, MOLECULAR_GENETICS.
--intogen-compendium
    The IntOGen "Compendium_Cancer_Genes" driver table (tab-separated),
    with at least the columns: SYMBOL, IS_DRIVER, COHORT, CANCER_TYPE,
    ROLE, METHODS, QVALUE_COMBINATION, SAMPLES, %_SAMPLES_COHORT.
--intogen-unfiltered
    IntOGen's unfiltered per-cohort driver-call table (tab-separated), with
    at least the columns: SYMBOL, DRIVER, COHORT, TIER, NUM_PAPERS, FILTER.

Usage
-----
    python driver_gene_annotation.py \\
        --genes data/candidate_genes.csv \\
        --cosmic-cgc data/driver_catalogs/cosmic_cancer_gene_census.tsv \\
        --intogen-compendium data/driver_catalogs/intogen_compendium_cancer_genes.tsv \\
        --intogen-unfiltered data/driver_catalogs/intogen_unfiltered_drivers.tsv \\
        --output-dir results/driver_gene_annotation
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

# Keywords used only to flag COSMIC CGC entries whose reported tumour
# type(s) mention a tissue or disease category of interest, as a
# convenience annotation. Left empty by default; populate it with terms
# relevant to your own study (e.g. specific tumour-type names) to enable
# the flag. It does not affect the CGC/IntOGen driver classification
# itself -- only the "CGC_Tissue_Keyword_Match" column.
TISSUE_KEYWORDS: list[str] = []


def clean_gene_symbol(value) -> str | float:
    if pd.isna(value):
        return np.nan
    return str(value).strip().upper()


def join_unique(series: pd.Series, separator: str = "; ") -> str:
    """Join unique, non-empty values while preserving first-seen order."""
    values: list[str] = []
    for value in series:
        if pd.isna(value):
            continue
        value = str(value).strip()
        if value and value not in values:
            values.append(value)
    return separator.join(values)


def numeric_min(series: pd.Series) -> float:
    numeric = pd.to_numeric(series, errors="coerce")
    return numeric.min() if numeric.notna().any() else np.nan


def numeric_max(series: pd.Series) -> float:
    numeric = pd.to_numeric(series, errors="coerce")
    return numeric.max() if numeric.notna().any() else np.nan


def contains_tissue_keyword(text) -> bool:
    """Conservative annotation flag; not proof of tissue-specific driver status."""
    if pd.isna(text):
        return False
    text = str(text).lower()
    return any(keyword in text for keyword in TISSUE_KEYWORDS)


def summarize_cosmic_cgc(genes: list[str], cosmic: pd.DataFrame) -> pd.DataFrame:
    records = []
    for gene in genes:
        subset = cosmic[cosmic["GENE_SYMBOL"] == gene]

        if subset.empty:
            records.append(
                {
                    "Gene": gene, "CGC_status": "No", "CGC_Tier": np.nan, "CGC_COSMIC_Gene_ID": "",
                    "CGC_Role_in_Cancer": "", "CGC_Mutation_Types": "", "CGC_Tumour_Types_Somatic": "",
                    "CGC_Tissue_Type": "", "CGC_Molecular_Genetics": "", "CGC_Tissue_Keyword_Match": "No",
                }
            )
            continue

        tumour_types = join_unique(subset["TUMOUR_TYPES_SOMATIC"])
        tissue_types = join_unique(subset["TISSUE_TYPE"])

        keyword_flag = contains_tissue_keyword(tumour_types)

        records.append(
            {
                "Gene": gene,
                "CGC_status": "Yes",
                "CGC_Tier": join_unique(subset["TIER"].astype(str)),
                "CGC_COSMIC_Gene_ID": join_unique(subset["COSMIC_GENE_ID"]),
                "CGC_Role_in_Cancer": join_unique(subset["ROLE_IN_CANCER"]),
                "CGC_Mutation_Types": join_unique(subset["MUTATION_TYPES"]),
                "CGC_Tumour_Types_Somatic": tumour_types,
                "CGC_Tissue_Type": tissue_types,
                "CGC_Molecular_Genetics": join_unique(subset["MOLECULAR_GENETICS"]),
                "CGC_Tissue_Keyword_Match": "Yes" if keyword_flag else "No",
            }
        )

    return pd.DataFrame(records)


def summarize_intogen_compendium(genes: list[str], intogen: pd.DataFrame) -> pd.DataFrame:
    records = []
    for gene in genes:
        subset = intogen[intogen["SYMBOL"] == gene]

        if subset.empty:
            records.append(
                {
                    "Gene": gene, "IntOGen_driver": "No", "IntOGen_Number_of_Cohorts": 0,
                    "IntOGen_Number_of_Cancer_Types": 0, "IntOGen_Cancer_Types": "", "IntOGen_Cohorts": "",
                    "IntOGen_Role": "", "IntOGen_Methods": "", "IntOGen_Min_Qvalue": np.nan,
                    "IntOGen_Max_Mutated_Samples_in_Cohort": np.nan, "IntOGen_Max_Percent_Samples_in_Cohort": np.nan,
                }
            )
            continue

        driver_subset = subset[subset["IS_DRIVER"] == True]  # noqa: E712

        records.append(
            {
                "Gene": gene,
                "IntOGen_driver": "Yes" if not driver_subset.empty else "No",
                "IntOGen_Number_of_Cohorts": driver_subset["COHORT"].nunique(),
                "IntOGen_Number_of_Cancer_Types": driver_subset["CANCER_TYPE"].nunique(),
                "IntOGen_Cancer_Types": join_unique(driver_subset["CANCER_TYPE"]),
                "IntOGen_Cohorts": join_unique(driver_subset["COHORT"]),
                "IntOGen_Role": join_unique(driver_subset["ROLE"]),
                "IntOGen_Methods": join_unique(driver_subset["METHODS"]),
                "IntOGen_Min_Qvalue": numeric_min(driver_subset["QVALUE_COMBINATION"]),
                "IntOGen_Max_Mutated_Samples_in_Cohort": numeric_max(driver_subset["SAMPLES"]),
                "IntOGen_Max_Percent_Samples_in_Cohort": numeric_max(driver_subset["%_SAMPLES_COHORT"]),
            }
        )

    return pd.DataFrame(records)


def summarize_intogen_unfiltered(genes: list[str], unfiltered: pd.DataFrame) -> pd.DataFrame:
    records = []
    for gene in genes:
        subset = unfiltered[unfiltered["SYMBOL"] == gene]

        if subset.empty:
            records.append(
                {
                    "Gene": gene, "IntOGen_Unfiltered_Any_Driver_Call": "No",
                    "IntOGen_Unfiltered_Driver_Cohorts": 0, "IntOGen_Unfiltered_Min_Tier": np.nan,
                    "IntOGen_Unfiltered_Max_Num_Papers": np.nan, "IntOGen_Unfiltered_Filters": "",
                }
            )
            continue

        driver_rows = subset[subset["DRIVER"] == True]  # noqa: E712

        records.append(
            {
                "Gene": gene,
                "IntOGen_Unfiltered_Any_Driver_Call": "Yes" if not driver_rows.empty else "No",
                "IntOGen_Unfiltered_Driver_Cohorts": driver_rows["COHORT"].nunique(),
                "IntOGen_Unfiltered_Min_Tier": numeric_min(driver_rows["TIER"]),
                "IntOGen_Unfiltered_Max_Num_Papers": numeric_max(driver_rows["NUM_PAPERS"]),
                "IntOGen_Unfiltered_Filters": join_unique(driver_rows["FILTER"]),
            }
        )

    return pd.DataFrame(records)


def classify_driver(row: pd.Series) -> str:
    cgc = row["CGC_status"] == "Yes"
    intogen_driver = row["IntOGen_driver"] == "Yes"

    if cgc and intogen_driver:
        return "Supported by CGC and IntOGen"
    if cgc:
        return "CGC only"
    if intogen_driver:
        return "IntOGen only"
    return "Not identified in CGC or IntOGen"


FINAL_COLUMNS = [
    "Gene",
    "Cancer_Driver_Identified", "Driver_Evidence_Category", "Independent_Dual_Source_Support",
    "CGC_status", "CGC_Tier", "CGC_Role_in_Cancer", "CGC_Mutation_Types", "CGC_Tumour_Types_Somatic",
    "CGC_Tissue_Type", "CGC_Tissue_Keyword_Match", "CGC_Molecular_Genetics", "CGC_COSMIC_Gene_ID",
    "IntOGen_driver", "IntOGen_Role", "IntOGen_Number_of_Cancer_Types", "IntOGen_Cancer_Types",
    "IntOGen_Number_of_Cohorts", "IntOGen_Cohorts", "IntOGen_Methods", "IntOGen_Min_Qvalue",
    "IntOGen_Max_Mutated_Samples_in_Cohort", "IntOGen_Max_Percent_Samples_in_Cohort",
    "IntOGen_Unfiltered_Any_Driver_Call", "IntOGen_Unfiltered_Driver_Cohorts",
    "IntOGen_Unfiltered_Min_Tier", "IntOGen_Unfiltered_Max_Num_Papers", "IntOGen_Unfiltered_Filters",
    "Literature_Evidence_Notes",
]

CATEGORY_PRIORITY = {
    "Supported by CGC and IntOGen": 1,
    "CGC only": 2,
    "IntOGen only": 3,
    "Not identified in CGC or IntOGen": 4,
}


def build_annotation_table(
    genes: pd.DataFrame, cosmic: pd.DataFrame, intogen_compendium: pd.DataFrame, intogen_unfiltered: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_list = genes["Gene"].tolist()

    master = (
        genes.merge(summarize_cosmic_cgc(gene_list, cosmic), on="Gene", how="left")
        .merge(summarize_intogen_compendium(gene_list, intogen_compendium), on="Gene", how="left")
        .merge(summarize_intogen_unfiltered(gene_list, intogen_unfiltered), on="Gene", how="left")
    )

    master["Driver_Evidence_Category"] = master.apply(classify_driver, axis=1)
    master["Cancer_Driver_Identified"] = np.where(
        (master["CGC_status"] == "Yes") | (master["IntOGen_driver"] == "Yes"), "Yes", "No"
    )
    master["Independent_Dual_Source_Support"] = np.where(
        (master["CGC_status"] == "Yes") & (master["IntOGen_driver"] == "Yes"), "Yes", "No"
    )

    # Placeholder for evidence gathered separately (e.g. a targeted
    # literature review); intentionally not inferred from CGC/IntOGen.
    master["Literature_Evidence_Notes"] = "To be evaluated"

    master["_priority"] = master["Driver_Evidence_Category"].map(CATEGORY_PRIORITY)
    master = master.sort_values(by=["_priority", "Gene"]).drop(columns="_priority").reset_index(drop=True)
    master = master[FINAL_COLUMNS]

    summary = pd.DataFrame(
        {
            "Metric": [
                "Total genes analyzed",
                "Genes present in COSMIC Cancer Gene Census",
                "Genes identified as drivers in IntOGen",
                "Genes supported by both CGC and IntOGen",
                "Genes supported by at least one driver resource",
                "CGC genes matching a configured tissue keyword",
            ],
            "Value": [
                len(master),
                (master["CGC_status"] == "Yes").sum(),
                (master["IntOGen_driver"] == "Yes").sum(),
                (master["Independent_Dual_Source_Support"] == "Yes").sum(),
                (master["Cancer_Driver_Identified"] == "Yes").sum(),
                (master["CGC_Tissue_Keyword_Match"] == "Yes").sum(),
            ],
        }
    )

    return master, summary


def format_excel_workbook(excel_path: Path) -> None:
    workbook = load_workbook(excel_path)
    for worksheet in workbook.worksheets:
        for cell in worksheet[1]:
            cell.font = Font(bold=True)
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        for row in worksheet.iter_rows(min_row=2):
            for cell in row:
                cell.alignment = Alignment(vertical="top", wrap_text=True)
        for column_cells in worksheet.columns:
            column_letter = get_column_letter(column_cells[0].column)
            max_length = max((len(str(c.value)) for c in column_cells if c.value is not None), default=0)
            worksheet.column_dimensions[column_letter].width = min(max(max_length + 2, 12), 45)
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions
    workbook.save(excel_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genes", type=Path, required=True, help="CSV with a single 'Gene' column.")
    parser.add_argument("--cosmic-cgc", type=Path, required=True)
    parser.add_argument("--intogen-compendium", type=Path, required=True)
    parser.add_argument("--intogen-unfiltered", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results/driver_gene_annotation"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(args.genes)
    genes["Gene"] = genes["Gene"].apply(clean_gene_symbol)
    genes = genes.drop_duplicates(subset=["Gene"]).reset_index(drop=True)

    cosmic = pd.read_csv(args.cosmic_cgc, sep="\t", low_memory=False)
    cosmic["GENE_SYMBOL"] = cosmic["GENE_SYMBOL"].apply(clean_gene_symbol)

    intogen_compendium = pd.read_csv(args.intogen_compendium, sep="\t", low_memory=False)
    intogen_compendium["SYMBOL"] = intogen_compendium["SYMBOL"].apply(clean_gene_symbol)

    intogen_unfiltered = pd.read_csv(args.intogen_unfiltered, sep="\t", low_memory=False)
    intogen_unfiltered["SYMBOL"] = intogen_unfiltered["SYMBOL"].apply(clean_gene_symbol)

    master, summary = build_annotation_table(genes, cosmic, intogen_compendium, intogen_unfiltered)

    csv_path = args.output_dir / "driver_gene_annotation.csv"
    excel_path = args.output_dir / "driver_gene_annotation.xlsx"

    master.to_csv(csv_path, index=False)

    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        master.to_excel(writer, sheet_name="Driver_Annotation", index=False)
        summary.to_excel(writer, sheet_name="Summary", index=False)

    format_excel_workbook(excel_path)

    print("=" * 70)
    print("DRIVER GENE ANNOTATION COMPLETE")
    print("=" * 70)
    print(summary.to_string(index=False))
    print("\nDriver evidence categories:")
    print(master["Driver_Evidence_Category"].value_counts().to_string())
    print(f"\nFiles saved:\n{csv_path}\n{excel_path}")


if __name__ == "__main__":
    main()
