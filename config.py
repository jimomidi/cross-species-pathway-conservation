"""
Shared configuration for the cross-species pathway conservation toolkit.

This module centralizes the directory layout and the small set of
labels/colors that are reused across multiple scripts. Every script in
this repository can be run independently; the defaults below simply
avoid repeating the same relative paths in every file.

Directory layout (all paths are relative to the repository root and are
created automatically by the scripts that need them):

    data/                       Place raw/input files here (not tracked by git).
        biomart/                Ortholog export from Ensembl BioMart.
        variants/               Per-sample VCF files (tumor and liquid-biopsy).
        annotation/             ANNOVAR exonic-variant-function files.
        reference_mutations/    Public reference mutation catalog (e.g. COSMIC export).
        enrichment/<species>/   DAVID functional-annotation exports, one folder per species.
        driver_catalogs/        COSMIC Cancer Gene Census / IntOGen driver tables.

    results/                    All generated tables and figures (not tracked by git).

Edit SPECIES_A / SPECIES_B if you are adapting this toolkit to a
different pair of organisms. The two labels are used throughout the
toolkit for axis labels, legends, and output file names.
"""

from __future__ import annotations

from pathlib import Path

# ---------------------------------------------------------------------------
# Directory layout
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent

DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

# ---------------------------------------------------------------------------
# Species labels
#
# This toolkit was designed for a two-species comparative-genomics study
# (a companion-animal species compared with human data). Replace these
# labels to reuse the toolkit for a different species pair; they are used
# only for display purposes (axis labels, legends, and file names).
# ---------------------------------------------------------------------------
SPECIES_A = "Canine"
SPECIES_B = "Human"

# KEGG organism codes and NCBI taxonomy IDs, used by the pathway-universe
# scripts. Update these if SPECIES_A/SPECIES_B are changed.
KEGG_ORGANISM_CODE = {
    SPECIES_A: "cfa",  # Canis lupus familiaris
    SPECIES_B: "hsa",  # Homo sapiens
}

NCBI_TAXONOMY_ID = {
    SPECIES_A: "9615",
    SPECIES_B: "9606",
}

# ---------------------------------------------------------------------------
# Shared figure style
# ---------------------------------------------------------------------------
FONT_FAMILY = "Times New Roman"

DATABASE_COLORS = {
    "GO BP": "#76B7E5",
    "KEGG": "#79C99E",
    "Reactome": "#C7A0D9",
}

DATABASE_MARKERS = {
    "GO BP": "o",
    "KEGG": "s",
    "Reactome": "^",
}


def ensure_dir(path: Path) -> Path:
    """Create ``path`` (including parents) if it does not already exist."""
    path.mkdir(parents=True, exist_ok=True)
    return path
