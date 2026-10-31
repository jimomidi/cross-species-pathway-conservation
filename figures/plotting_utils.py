"""
Shared helpers for the pathway-conservation figure scripts in this folder.

All figure scripts consume the same input: the three per-database
"conserved shared pathways" tables produced by
``pathway_enrichment/compare_enrichment_results.py``
(``GO_BP_Conserved_Shared_Pathways.csv``, ``KEGG_Conserved_Shared_Pathways.csv``,
``REACTOME_Conserved_Shared_Pathways.csv``), each containing a ``Pathway``
column and ``<species> Fold Enrichment`` columns for the two species
configured in ``config.py``.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_DATABASE_FILES = {
    "GO BP": "GO_BP_Conserved_Shared_Pathways.csv",
    "KEGG": "KEGG_Conserved_Shared_Pathways.csv",
    "Reactome": "REACTOME_Conserved_Shared_Pathways.csv",
}


def read_csv_robustly(file_path: Path) -> pd.DataFrame:
    """Read a CSV/TSV file, automatically detecting the delimiter."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file was not found:\n{file_path}")

    last_error = None
    for settings in [{"sep": None, "engine": "python"}, {"sep": ","}, {"sep": "\t"}]:
        try:
            df = pd.read_csv(file_path, **settings)
            if df.shape[1] > 1:
                return df
        except Exception as error:  # noqa: BLE001
            last_error = error

    raise ValueError(f"The file could not be read correctly:\n{file_path}\nLast error: {last_error}")


def load_conserved_pathways(
    input_dir: Path,
    species_a: str,
    species_b: str,
    database_files: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Load and combine the per-database conserved-pathway tables.

    Numeric fold-enrichment columns are coerced, rows with missing/invalid
    values are dropped, and exact duplicate rows are removed.
    """
    database_files = database_files or DEFAULT_DATABASE_FILES
    a_col, b_col = f"{species_a} Fold Enrichment", f"{species_b} Fold Enrichment"
    required_columns = ["Pathway", a_col, b_col]

    tables = []
    for database_name, filename in database_files.items():
        df = read_csv_robustly(Path(input_dir) / filename)
        df.columns = df.columns.astype(str).str.strip()

        missing = [c for c in required_columns if c not in df.columns]
        if missing:
            raise KeyError(f"Missing columns {missing} in {filename}. Available: {df.columns.tolist()}")

        df = df.copy()
        df["Database"] = database_name
        tables.append(df)

    combined = pd.concat(tables, ignore_index=True)

    for column in (a_col, b_col):
        combined[column] = pd.to_numeric(combined[column], errors="coerce")

    combined = combined.dropna(subset=["Pathway", a_col, b_col]).copy()
    combined["Pathway"] = combined["Pathway"].astype(str).str.strip()
    combined = combined[(combined[a_col] > 0) & (combined[b_col] > 0)].copy()
    combined = combined.drop_duplicates(subset=["Database", "Pathway", a_col, b_col]).reset_index(drop=True)

    if "Mean Fold Enrichment" not in combined.columns:
        combined["Mean Fold Enrichment"] = combined[[a_col, b_col]].mean(axis=1)

    return combined


def parse_gene_ids(value) -> set[str]:
    """Parse a delimited "User Ids" field into a set of uppercase gene identifiers."""
    import re

    if pd.isna(value):
        return set()

    text = str(value).strip()
    if not text:
        return set()

    text = re.sub(r"[;|/\n\r\t]+", ",", text)

    elements: list[str] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        subitems = item.split()
        elements.extend(subitems) if len(subitems) > 1 else elements.append(item)

    return {element.strip().upper() for element in elements if element.strip()}


def scale_values(values, new_min: float, new_max: float) -> np.ndarray:
    """Linearly rescale numeric values into [new_min, new_max]."""
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]

    if len(finite) == 0:
        return np.full(len(values), (new_min + new_max) / 2)

    old_min, old_max = np.nanmin(finite), np.nanmax(finite)
    if np.isclose(old_min, old_max):
        return np.full(len(values), (new_min + new_max) / 2)

    return (values - old_min) / (old_max - old_min) * (new_max - new_min) + new_min


def save_figure(fig, base_path: Path, formats: tuple[str, ...] = ("png", "pdf", "svg"), dpi: int = 600) -> list[Path]:
    """Save a figure in each of the requested formats, returning the paths written."""
    base_path = Path(base_path)
    base_path.parent.mkdir(parents=True, exist_ok=True)
    written = []

    for fmt in formats:
        output_path = base_path.with_suffix(f".{fmt}")
        save_kwargs = {"bbox_inches": "tight"}
        if fmt in ("png", "tiff"):
            save_kwargs["dpi"] = dpi
        if fmt == "tiff":
            save_kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(output_path, **save_kwargs)
        written.append(output_path)

    return written


def apply_publication_style(font_family: str = "Times New Roman") -> None:
    """Apply the shared matplotlib style used across the figure scripts."""
    import matplotlib.pyplot as plt

    plt.rcParams["font.family"] = font_family
    plt.rcParams["font.size"] = 12
    plt.rcParams["axes.labelsize"] = 15
    plt.rcParams["axes.titlesize"] = 16
    plt.rcParams["xtick.labelsize"] = 12
    plt.rcParams["ytick.labelsize"] = 12
    plt.rcParams["legend.fontsize"] = 11
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
