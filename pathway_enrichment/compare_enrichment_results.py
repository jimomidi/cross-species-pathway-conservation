"""
Compare DAVID functional-annotation results between two species and
identify conserved (shared), species-A-specific, and species-B-specific
pathways.

For each of GO Biological Process, KEGG, and Reactome, this script:

  1. Reads the DAVID export for each species (one CSV per database, per
     species).
  2. Normalizes pathway/term names so that species-specific identifier
     prefixes (e.g. "hsa04151", "cfa04151", "GO:0006355~") do not prevent
     matching equivalent terms across species.
  3. Identifies shared, species-A-specific, and species-B-specific
     pathways, and computes the Jaccard index of overlap.
  4. For shared pathways, computes Pearson and Spearman correlation of
     fold enrichment between species and saves a Venn diagram and a
     fold-enrichment scatter plot.
  5. Writes per-database CSV tables and one combined, formatted Excel
     workbook.

Note on filtering: this script analyzes every pathway retained in the
input DAVID files; it does not apply an additional P-value/Benjamini/FDR
or gene-count filter. If your DAVID export already includes non-
significant terms and you want to restrict the comparison to significant
ones, filter the input files (or extend ``prepare_david_data``) before
running this script.

Usage
-----
    python compare_enrichment_results.py \\
        --species-a-dir data/enrichment/canine \\
        --species-b-dir data/enrichment/human \\
        --output-dir results/pathway_conservation
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib_venn import venn2
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy.stats import pearsonr, spearmanr

import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import FONT_FAMILY, SPECIES_A, SPECIES_B  # noqa: E402

DATABASES = {
    "GO_BP": {"filename": "go_bp.csv", "display_name": "GO Biological Process"},
    "KEGG": {"filename": "kegg.csv", "display_name": "KEGG"},
    "REACTOME": {"filename": "reactome.csv", "display_name": "Reactome"},
}

FIGURE_DPI = 600
SAVE_PDF = True

SPECIES_A_COLOR = "#87CEEB"
SPECIES_B_COLOR = "#FFB7B2"
SHARED_COLOR = "#A78BFA"

REQUIRED_COLUMNS = {"Term", "Count", "P-Value", "Fold Enrichment"}

NUMERIC_COLUMNS = [
    "Count", "List Total", "Pop Hits", "Pop Total",
    "P-Value", "Benjamini", "Fold Enrichment", "Bonferroni", "FDR", "Fisher Exact",
]

SHARED_SOURCE_COLUMNS = [
    "Matching Key", "Normalized Term", "Category", "Count", "List Total",
    "Pop Hits", "Pop Total", "P-Value", "Benjamini", "FDR",
    "Fold Enrichment", "Fold Enrichment Rank", "User Ids",
]


def read_david_file(file_path: Path) -> pd.DataFrame:
    """Read a DAVID export robustly, trying several delimiters and encodings."""
    if not file_path.exists():
        raise FileNotFoundError(f"Input file was not found:\n{file_path}")

    last_error = None
    for encoding in ["utf-8-sig", "utf-8", "latin-1"]:
        try:
            candidate = pd.read_csv(file_path, sep=None, engine="python", encoding=encoding)
            if candidate.shape[1] > 1:
                candidate.columns = candidate.columns.astype(str).str.strip()
                return candidate
        except Exception as error:  # noqa: BLE001
            last_error = error

    raise ValueError(f"Could not read file:\n{file_path}\n\nLast error: {last_error}")


def validate_columns(df: pd.DataFrame, file_path: Path) -> None:
    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        raise ValueError(
            f"The following required columns are missing from:\n{file_path}\n\n"
            f"Missing columns: {sorted(missing)}\n\nAvailable columns: {list(df.columns)}"
        )


def clean_numeric_column(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.astype(str).str.replace(",", "", regex=False).str.strip(), errors="coerce")


def choose_term_text(row: pd.Series) -> str:
    """Prefer "Full Term" when available and non-empty; otherwise use "Term"."""
    if "Full Term" in row.index:
        full_term = row["Full Term"]
        if pd.notna(full_term):
            full_term = str(full_term).strip()
            if full_term and full_term.lower() != "nan":
                return full_term
    return str(row["Term"]).strip()


def normalize_term_name(term: str, database_name: str) -> str | float:
    """Strip species-specific identifier prefixes so terms match across species.

    Examples removed: "GO:0006355~", "hsa04151:", "cfa04151:", "R-HSA-123456~".
    """
    if pd.isna(term):
        return np.nan

    term = str(term).strip()

    if "~" in term:
        parts = term.split("~", 1)
        if len(parts) == 2 and parts[1].strip():
            term = parts[1].strip()

    if database_name == "GO_BP":
        term = re.sub(r"^GO:\d+\s*[:\-–—]?\s*", "", term, flags=re.IGNORECASE)
    elif database_name == "KEGG":
        term = re.sub(r"^[A-Za-z]{2,5}\d{4,6}\s*[:\-–—]?\s*", "", term)
        term = re.sub(r"^\d{4,6}\s*[:\-–—]?\s*", "", term)
    elif database_name == "REACTOME":
        term = re.sub(r"^R-[A-Za-z]+-\d+\s*[:\-–—]?\s*", "", term, flags=re.IGNORECASE)

    term = re.sub(r"\s+", " ", term.strip()).rstrip(" .;:")
    return term


def create_matching_key(term: str) -> str | float:
    if pd.isna(term):
        return np.nan
    key = str(term).strip().lower().replace("–", "-").replace("—", "-").replace("&", "and")
    key = re.sub(r"\s+", " ", key).rstrip(" .;:")
    return key


def prepare_david_data(df: pd.DataFrame, species: str, database_name: str) -> pd.DataFrame:
    """Clean and standardize a DAVID export for cross-species matching.

    All rows with a valid pathway name and positive fold enrichment are
    retained; see the module docstring regarding significance filtering.
    """
    df = df.copy()

    for column in NUMERIC_COLUMNS:
        if column in df.columns:
            df[column] = clean_numeric_column(df[column])

    df["Original Term"] = df.apply(choose_term_text, axis=1)
    df["Normalized Term"] = df["Original Term"].apply(lambda v: normalize_term_name(v, database_name))
    df["Matching Key"] = df["Normalized Term"].apply(create_matching_key)
    df["Species"] = species
    df["Database"] = database_name

    df = df.dropna(subset=["Matching Key", "Fold Enrichment"]).copy()
    df = df[np.isfinite(df["Fold Enrichment"]) & (df["Fold Enrichment"] > 0)].copy()

    # If duplicate normalized pathway names occur, keep the row with the
    # lowest P-value, then the highest fold enrichment.
    df = df.sort_values(by=["P-Value", "Fold Enrichment"], ascending=[True, False], na_position="last")
    df = df.drop_duplicates(subset="Matching Key", keep="first").copy()

    df["Fold Enrichment Rank"] = df["Fold Enrichment"].rank(ascending=False, method="min").astype(int)

    return df.reset_index(drop=True)


def safe_correlation(x: pd.Series, y: pd.Series, method: str = "pearson") -> tuple[float, float]:
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    valid = x.notna() & y.notna()
    x, y = x[valid], y[valid]

    if len(x) < 3 or x.nunique() < 2 or y.nunique() < 2:
        return np.nan, np.nan

    try:
        result = pearsonr(x, y) if method == "pearson" else spearmanr(x, y)
        return float(result.statistic), float(result.pvalue)
    except Exception:  # noqa: BLE001
        return np.nan, np.nan


def format_statistic(value: float, digits: int = 3) -> str:
    return "NA" if pd.isna(value) else f"{value:.{digits}f}"


def format_p_value(p_value: float) -> str:
    if pd.isna(p_value):
        return "NA"
    return f"{p_value:.2e}" if p_value < 0.001 else f"{p_value:.3f}"


def save_venn_plot(a_terms: set, b_terms: set, database_display_name: str, figure_dir: Path, output_prefix: str) -> None:
    fig, ax = plt.subplots(figsize=(6, 5.5))

    venn = venn2(
        subsets=(set(a_terms), set(b_terms)),
        set_labels=(f"{SPECIES_A} enriched pathways", f"{SPECIES_B} enriched pathways"),
        set_colors=(SPECIES_A_COLOR, SPECIES_B_COLOR),
        alpha=0.75,
        ax=ax,
    )

    for label in (venn.set_labels or []):
        if label is not None:
            label.set_fontname(FONT_FAMILY)
            label.set_fontsize(14)
            label.set_fontweight("bold")
    for label in (venn.subset_labels or []):
        if label is not None:
            label.set_fontname(FONT_FAMILY)
            label.set_fontsize(16)
            label.set_fontweight("bold")

    ax.set_title(f"{database_display_name} pathway overlap", fontsize=17, fontweight="bold", pad=12)
    plt.tight_layout()

    plt.savefig(figure_dir / f"{output_prefix}_Venn.png", dpi=FIGURE_DPI, bbox_inches="tight")
    if SAVE_PDF:
        plt.savefig(figure_dir / f"{output_prefix}_Venn.pdf", bbox_inches="tight")
    plt.close(fig)


def save_scatter_plot(
    shared_df: pd.DataFrame,
    database_display_name: str,
    figure_dir: Path,
    output_prefix: str,
    pearson_r: float,
    pearson_p: float,
    spearman_rho: float,
    spearman_p: float,
) -> None:
    if shared_df.empty:
        print(f"Warning: no shared pathways for {database_display_name}; scatter plot was not created.")
        return

    x = shared_df[f"{SPECIES_A} Fold Enrichment"]
    y = shared_df[f"{SPECIES_B} Fold Enrichment"]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=90, color=SHARED_COLOR, alpha=0.80, edgecolors="black", linewidth=0.8)

    minimum_value, maximum_value = min(x.min(), y.min()), max(x.max(), y.max())
    total_range = maximum_value - minimum_value
    padding = (total_range * 0.08) if total_range else max(maximum_value * 0.10, 1)
    lower_limit, upper_limit = max(0, minimum_value - padding), maximum_value + padding

    ax.plot([lower_limit, upper_limit], [lower_limit, upper_limit], color="black", linewidth=0.9, dashes=(10, 6), label="y = x")
    ax.set_xlim(lower_limit, upper_limit)
    ax.set_ylim(lower_limit, upper_limit)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel(f"{SPECIES_A} Fold Enrichment", fontsize=16, fontweight="bold", labelpad=10)
    ax.set_ylabel(f"{SPECIES_B} Fold Enrichment", fontsize=16, fontweight="bold", labelpad=10)
    ax.set_title(f"{database_display_name} enrichment conservation", fontsize=17, fontweight="bold", pad=12)

    correlation_text = (
        f"Shared pathways: n = {len(shared_df)}\n"
        f"Pearson r = {format_statistic(pearson_r)}, P = {format_p_value(pearson_p)}\n"
        f"Spearman ρ = {format_statistic(spearman_rho)}, P = {format_p_value(spearman_p)}"
    )
    ax.text(
        0.04, 0.96, correlation_text, transform=ax.transAxes, ha="left", va="top", fontsize=11.5,
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "edgecolor": "black", "alpha": 0.85},
    )

    ax.tick_params(axis="both", which="major", labelsize=13, length=5, width=1.1)
    ax.grid(alpha=0.18, linewidth=0.7)
    ax.legend(frameon=False, fontsize=11, loc="lower right")
    plt.tight_layout()

    plt.savefig(figure_dir / f"{output_prefix}_Fold_Enrichment_Scatter.png", dpi=FIGURE_DPI, bbox_inches="tight")
    if SAVE_PDF:
        plt.savefig(figure_dir / f"{output_prefix}_Fold_Enrichment_Scatter.pdf", bbox_inches="tight")
    plt.close(fig)


def analyze_database(
    a_df: pd.DataFrame, b_df: pd.DataFrame, database_name: str, database_display_name: str, figure_dir: Path
) -> dict:
    """Run the full pathway-conservation comparison for one database."""
    print("\n" + "=" * 70)
    print(f"Analyzing: {database_display_name}")
    print("=" * 70)

    a_keys, b_keys = set(a_df["Matching Key"]), set(b_df["Matching Key"])
    shared_keys = a_keys & b_keys
    union_keys = a_keys | b_keys

    jaccard_index = len(shared_keys) / len(union_keys) if union_keys else np.nan

    rename_a = {
        "Normalized Term": f"{SPECIES_A} Term", "Category": f"{SPECIES_A} Category", "Count": f"{SPECIES_A} Count",
        "List Total": f"{SPECIES_A} List Total", "Pop Hits": f"{SPECIES_A} Pop Hits", "Pop Total": f"{SPECIES_A} Pop Total",
        "P-Value": f"{SPECIES_A} P-Value", "Benjamini": f"{SPECIES_A} Benjamini", "FDR": f"{SPECIES_A} FDR",
        "Fold Enrichment": f"{SPECIES_A} Fold Enrichment", "Fold Enrichment Rank": f"{SPECIES_A} Rank",
        "User Ids": f"{SPECIES_A} User Ids",
    }
    rename_b = {
        "Normalized Term": f"{SPECIES_B} Term", "Category": f"{SPECIES_B} Category", "Count": f"{SPECIES_B} Count",
        "List Total": f"{SPECIES_B} List Total", "Pop Hits": f"{SPECIES_B} Pop Hits", "Pop Total": f"{SPECIES_B} Pop Total",
        "P-Value": f"{SPECIES_B} P-Value", "Benjamini": f"{SPECIES_B} Benjamini", "FDR": f"{SPECIES_B} FDR",
        "Fold Enrichment": f"{SPECIES_B} Fold Enrichment", "Fold Enrichment Rank": f"{SPECIES_B} Rank",
        "User Ids": f"{SPECIES_B} User Ids",
    }

    a_shared = a_df[a_df["Matching Key"].isin(shared_keys)][
        [c for c in SHARED_SOURCE_COLUMNS if c in a_df.columns]
    ].rename(columns=rename_a)
    b_shared = b_df[b_df["Matching Key"].isin(shared_keys)][
        [c for c in SHARED_SOURCE_COLUMNS if c in b_df.columns]
    ].rename(columns=rename_b)

    shared_df = pd.merge(a_shared, b_shared, on="Matching Key", how="inner")

    if not shared_df.empty:
        shared_df["Pathway"] = shared_df[f"{SPECIES_A} Term"].fillna(shared_df[f"{SPECIES_B} Term"])
        shared_df["Rank Difference"] = (shared_df[f"{SPECIES_A} Rank"] - shared_df[f"{SPECIES_B} Rank"]).abs()
        shared_df["Mean Fold Enrichment"] = (
            shared_df[f"{SPECIES_A} Fold Enrichment"] + shared_df[f"{SPECIES_B} Fold Enrichment"]
        ) / 2
        shared_df[f"Fold Enrichment Ratio {SPECIES_A}/{SPECIES_B}"] = (
            shared_df[f"{SPECIES_A} Fold Enrichment"] / shared_df[f"{SPECIES_B} Fold Enrichment"]
        )
        shared_df = shared_df.sort_values(by=["Mean Fold Enrichment", "Rank Difference"], ascending=[False, True])

        desired_order = [
            "Pathway", "Matching Key", f"{SPECIES_A} Fold Enrichment", f"{SPECIES_B} Fold Enrichment",
            "Mean Fold Enrichment", f"Fold Enrichment Ratio {SPECIES_A}/{SPECIES_B}",
            f"{SPECIES_A} Rank", f"{SPECIES_B} Rank", "Rank Difference",
            f"{SPECIES_A} Count", f"{SPECIES_B} Count", f"{SPECIES_A} List Total", f"{SPECIES_B} List Total",
            f"{SPECIES_A} Pop Hits", f"{SPECIES_B} Pop Hits", f"{SPECIES_A} Pop Total", f"{SPECIES_B} Pop Total",
            f"{SPECIES_A} P-Value", f"{SPECIES_B} P-Value", f"{SPECIES_A} Benjamini", f"{SPECIES_B} Benjamini",
            f"{SPECIES_A} FDR", f"{SPECIES_B} FDR", f"{SPECIES_A} User Ids", f"{SPECIES_B} User Ids",
        ]
        shared_df = shared_df[[c for c in desired_order if c in shared_df.columns]]

    a_specific_df = a_df[a_df["Matching Key"].isin(a_keys - b_keys)].copy()
    b_specific_df = b_df[b_df["Matching Key"].isin(b_keys - a_keys)].copy()
    a_specific_df["Conservation Class"] = f"{SPECIES_A}-specific"
    b_specific_df["Conservation Class"] = f"{SPECIES_B}-specific"

    if not shared_df.empty:
        pearson_r, pearson_p = safe_correlation(shared_df[f"{SPECIES_A} Fold Enrichment"], shared_df[f"{SPECIES_B} Fold Enrichment"], "pearson")
        spearman_rho, spearman_p = safe_correlation(shared_df[f"{SPECIES_A} Fold Enrichment"], shared_df[f"{SPECIES_B} Fold Enrichment"], "spearman")
    else:
        pearson_r = pearson_p = spearman_rho = spearman_p = np.nan

    summary_df = pd.DataFrame(
        {
            "Database": [database_display_name],
            f"{SPECIES_A} Retained Pathways": [len(a_keys)],
            f"{SPECIES_B} Retained Pathways": [len(b_keys)],
            "Shared Conserved Pathways": [len(shared_keys)],
            f"{SPECIES_A}-Specific Pathways": [len(a_keys - b_keys)],
            f"{SPECIES_B}-Specific Pathways": [len(b_keys - a_keys)],
            "Union Pathways": [len(union_keys)],
            "Jaccard Index": [jaccard_index],
            f"Shared Percentage of {SPECIES_A} Pathways": [100 * len(shared_keys) / len(a_keys) if a_keys else np.nan],
            f"Shared Percentage of {SPECIES_B} Pathways": [100 * len(shared_keys) / len(b_keys) if b_keys else np.nan],
            "Pearson r": [pearson_r],
            "Pearson P-Value": [pearson_p],
            "Spearman rho": [spearman_rho],
            "Spearman P-Value": [spearman_p],
        }
    )

    save_venn_plot(a_keys, b_keys, database_display_name, figure_dir, database_name)
    save_scatter_plot(shared_df, database_display_name, figure_dir, database_name, pearson_r, pearson_p, spearman_rho, spearman_p)

    print(f"{SPECIES_A} retained pathways: {len(a_keys)}")
    print(f"{SPECIES_B} retained pathways: {len(b_keys)}")
    print(f"Shared pathways: {len(shared_keys)}")
    print(f"Jaccard index: {jaccard_index:.4f}" if pd.notna(jaccard_index) else "Jaccard index: not calculable")

    return {
        "summary": summary_df,
        f"{SPECIES_A.lower()}_retained": a_df,
        f"{SPECIES_B.lower()}_retained": b_df,
        "shared": shared_df,
        "a_specific": a_specific_df,
        "b_specific": b_specific_df,
    }


def format_excel_workbook(excel_path: Path) -> None:
    workbook = load_workbook(excel_path)
    header_fill = PatternFill(fill_type="solid", fgColor="D9CCFF")
    header_font = Font(name=FONT_FAMILY, bold=True, size=11)
    body_font = Font(name=FONT_FAMILY, size=10)

    for worksheet in workbook.worksheets:
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions

        for cell in worksheet[1]:
            cell.fill = header_fill
            cell.font = header_font
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

        for row in worksheet.iter_rows(min_row=2):
            for cell in row:
                cell.font = body_font
                cell.alignment = Alignment(vertical="top", wrap_text=False)

        for column_cells in worksheet.columns:
            column_letter = get_column_letter(column_cells[0].column)
            max_length = max((len(str(c.value)) for c in column_cells if c.value is not None), default=0)
            worksheet.column_dimensions[column_letter].width = max(min(max_length + 2, 55), 12)

    workbook.save(excel_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--species-a-dir", type=Path, default=Path(f"data/enrichment/{SPECIES_A.lower()}"))
    parser.add_argument("--species-b-dir", type=Path, default=Path(f"data/enrichment/{SPECIES_B.lower()}"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/pathway_conservation"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    figure_dir = args.output_dir / "figures"
    table_dir = args.output_dir / "tables"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams["font.family"] = FONT_FAMILY
    plt.rcParams["axes.linewidth"] = 1.2
    plt.rcParams["xtick.major.width"] = 1.1
    plt.rcParams["ytick.major.width"] = 1.1

    all_results = {}
    all_summaries = []

    for database_name, info in DATABASES.items():
        a_raw = read_david_file(args.species_a_dir / info["filename"])
        b_raw = read_david_file(args.species_b_dir / info["filename"])
        validate_columns(a_raw, args.species_a_dir / info["filename"])
        validate_columns(b_raw, args.species_b_dir / info["filename"])

        a_clean = prepare_david_data(a_raw, SPECIES_A, database_name)
        b_clean = prepare_david_data(b_raw, SPECIES_B, database_name)

        results = analyze_database(a_clean, b_clean, database_name, info["display_name"], figure_dir)
        all_results[database_name] = results
        all_summaries.append(results["summary"])

        for key, suffix in [
            (f"{SPECIES_A.lower()}_retained", f"{SPECIES_A}_All_Retained"),
            (f"{SPECIES_B.lower()}_retained", f"{SPECIES_B}_All_Retained"),
            ("shared", "Conserved_Shared_Pathways"),
            ("a_specific", f"{SPECIES_A}_Specific_Pathways"),
            ("b_specific", f"{SPECIES_B}_Specific_Pathways"),
            ("summary", "Summary"),
        ]:
            results[key].to_csv(table_dir / f"{database_name}_{suffix}.csv", index=False)

    combined_summary = pd.concat(all_summaries, ignore_index=True)
    combined_summary.to_csv(table_dir / "Combined_Pathway_Conservation_Summary.csv", index=False)

    excel_path = args.output_dir / "Pathway_Conservation_Analysis.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        combined_summary.to_excel(writer, sheet_name="Combined_Summary", index=False)
        for database_name, results in all_results.items():
            results["summary"].to_excel(writer, sheet_name=f"{database_name}_Summary", index=False)
            results["shared"].to_excel(writer, sheet_name=f"{database_name}_Shared", index=False)
            results["a_specific"].to_excel(writer, sheet_name=f"{database_name}_{SPECIES_A[:6]}Only", index=False)
            results["b_specific"].to_excel(writer, sheet_name=f"{database_name}_{SPECIES_B[:6]}Only", index=False)
            results[f"{SPECIES_A.lower()}_retained"].to_excel(writer, sheet_name=f"{database_name}_{SPECIES_A[:6]}All", index=False)
            results[f"{SPECIES_B.lower()}_retained"].to_excel(writer, sheet_name=f"{database_name}_{SPECIES_B[:6]}All", index=False)

    format_excel_workbook(excel_path)

    print("\n" + "=" * 70)
    print("PATHWAY CONSERVATION ANALYSIS COMPLETED")
    print("=" * 70)
    print(f"\nExcel workbook:\n{excel_path}")
    print(f"Figures folder:\n{figure_dir}")
    print(f"Tables folder:\n{table_dir}")
    print("\nCombined summary:")
    print(combined_summary.to_string(index=False))


if __name__ == "__main__":
    main()
