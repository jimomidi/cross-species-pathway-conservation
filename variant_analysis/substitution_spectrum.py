"""
Compare the 12-class single-nucleotide substitution spectrum between a
reference mutation catalog (species B) and paired tumor / liquid-biopsy
VCF files (species A).

For each species-A sample, only variants shared between its tumor and
liquid-biopsy (e.g. cell-free DNA) VCF files are counted, on the
assumption that variants confirmed in both sample types are the most
reliable calls. The same mutation observed in two different samples is
counted twice, since the samples are independent.

For the species-B reference catalog (for example, an exported table from
a public somatic-mutation database such as COSMIC), the substitution type
is extracted from a CDS-mutation notation column (e.g. ``c.17G>T``), and
each genomic mutation is counted once per sample (duplicate transcript
annotations of the same mutation are removed).

Input layout
------------
VCF files are discovered automatically from ``--vcf-dir`` using the naming
convention ``<sample>_Tumor.vcf.gz`` / ``<sample>_cfDNA.vcf.gz``; any sample
missing one of the two files is skipped with a warning.

Usage
-----
    python substitution_spectrum.py \\
        --reference-mutations data/reference_mutations/reference_snvs.csv \\
        --vcf-dir data/variants \\
        --output-dir results/substitution_spectrum
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import FONT_FAMILY, SPECIES_A, SPECIES_B  # noqa: E402
from vcf_utils import filter_snvs, read_vcf_variants  # noqa: E402

MUTATION_ORDER = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G",
]

TUMOR_SUFFIX = "_Tumor.vcf.gz"
LIQUID_BIOPSY_SUFFIX = "_cfDNA.vcf.gz"


def extract_substitution(cds_mutation: str) -> str | None:
    """Extract a base-substitution class (e.g. "A>C") from CDS notation."""
    if pd.isna(cds_mutation):
        return None

    match = re.search(r"([ACGT])>([ACGT])", str(cds_mutation).upper())
    if match:
        substitution = f"{match.group(1)}>{match.group(2)}"
        if substitution in MUTATION_ORDER:
            return substitution

    return None


def summarize_reference_catalog(reference_table: pd.DataFrame) -> pd.Series:
    """Compute the 12-class substitution counts for the reference catalog.

    Multiple transcript annotations of the same genomic mutation in the
    same sample are collapsed to a single event before counting.
    """
    table = reference_table.copy()
    table["Mutation_Type"] = table["CDS Mutation"].apply(extract_substitution)
    table = table[table["Mutation_Type"].isin(MUTATION_ORDER)]

    table = table.drop_duplicates(subset=["Sample ID", "Genomic Co-ordinates"])

    counts = table["Mutation_Type"].value_counts().reindex(MUTATION_ORDER, fill_value=0)
    return counts


def discover_sample_pairs(vcf_dir: Path) -> dict[str, tuple[Path, Path]]:
    """Find matching tumor / liquid-biopsy VCF pairs in ``vcf_dir``."""
    pairs: dict[str, tuple[Path, Path]] = {}

    for tumor_path in sorted(vcf_dir.glob(f"*{TUMOR_SUFFIX}")):
        sample = tumor_path.name[: -len(TUMOR_SUFFIX)]
        liquid_biopsy_path = vcf_dir / f"{sample}{LIQUID_BIOPSY_SUFFIX}"

        if liquid_biopsy_path.exists():
            pairs[sample] = (tumor_path, liquid_biopsy_path)
        else:
            print(f"Warning: no matching liquid-biopsy VCF for sample '{sample}'; skipping.")

    return pairs


def summarize_shared_variants(sample_pairs: dict[str, tuple[Path, Path]]) -> pd.Series:
    """Compute the 12-class substitution counts for tumor-vs-liquid-biopsy shared SNVs."""
    counts: Counter[str] = Counter()

    for sample, (tumor_path, liquid_biopsy_path) in sample_pairs.items():
        tumor_snvs = filter_snvs(read_vcf_variants(tumor_path))
        liquid_biopsy_snvs = filter_snvs(read_vcf_variants(liquid_biopsy_path))
        shared = tumor_snvs & liquid_biopsy_snvs

        print(
            f"{sample}: tumor={len(tumor_snvs):,} | "
            f"liquid biopsy={len(liquid_biopsy_snvs):,} | shared={len(shared):,}"
        )

        for _, _, ref, alt in shared:
            substitution = f"{ref}>{alt}"
            if substitution in MUTATION_ORDER:
                counts[substitution] += 1

    return pd.Series(counts).reindex(MUTATION_ORDER, fill_value=0)


def build_comparison_table(reference_counts: pd.Series, species_a_counts: pd.Series) -> pd.DataFrame:
    reference_total = reference_counts.sum()
    species_a_total = species_a_counts.sum()

    return pd.DataFrame(
        {
            "Mutation": MUTATION_ORDER,
            f"{SPECIES_B}_Count": reference_counts.values,
            f"{SPECIES_B}_Percent": (reference_counts / reference_total * 100).values if reference_total else 0,
            f"{SPECIES_A}_Shared_Count": species_a_counts.values,
            f"{SPECIES_A}_Shared_Percent": (species_a_counts / species_a_total * 100).values if species_a_total else 0,
        }
    )


def plot_comparison(comparison: pd.DataFrame, output_path: Path) -> None:
    """Save a publication-style bar chart of the substitution spectrum comparison."""
    plot_df = comparison.set_index("Mutation")[[f"{SPECIES_B}_Percent", f"{SPECIES_A}_Shared_Percent"]].copy()
    plot_df.columns = [SPECIES_B, f"{SPECIES_A} (tumor ∩ liquid biopsy)"]

    plt.rcParams["font.family"] = FONT_FAMILY

    fig, ax = plt.subplots(figsize=(11, 5.5))
    plot_df.plot(kind="bar", ax=ax, width=0.75)

    ax.set_xlabel("Base Substitution")
    ax.set_ylabel("Percentage of SNVs (%)")
    ax.set_title(f"Base-Substitution Spectrum in {SPECIES_B} and {SPECIES_A} Samples")
    ax.set_xticklabels(MUTATION_ORDER, rotation=0)
    ax.legend(title="")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--reference-mutations",
        type=Path,
        default=Path("data/reference_mutations/reference_snvs.csv"),
        help="CSV export from a public mutation catalog with 'CDS Mutation', 'Sample ID', "
        "and 'Genomic Co-ordinates' columns (default: %(default)s).",
    )
    parser.add_argument(
        "--vcf-dir",
        type=Path,
        default=Path("data/variants"),
        help="Directory containing '<sample>_Tumor.vcf.gz' / '<sample>_cfDNA.vcf.gz' pairs "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/substitution_spectrum"),
        help="Directory for output tables and figures (default: %(default)s).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading reference mutation catalog for {SPECIES_B}...")
    reference_table = pd.read_csv(args.reference_mutations, low_memory=False)
    reference_counts = summarize_reference_catalog(reference_table)

    print(f"\nReading {SPECIES_A} VCF files from {args.vcf_dir}...")
    sample_pairs = discover_sample_pairs(args.vcf_dir)
    if not sample_pairs:
        raise FileNotFoundError(f"No tumor/liquid-biopsy VCF pairs were found in {args.vcf_dir}")

    species_a_counts = summarize_shared_variants(sample_pairs)

    comparison = build_comparison_table(reference_counts, species_a_counts)
    comparison_path = args.output_dir / "substitution_spectrum_comparison.csv"
    comparison.to_csv(comparison_path, index=False)

    figure_path = args.output_dir / "substitution_spectrum_comparison.png"
    plot_comparison(comparison, figure_path)

    print("\n" + "=" * 60)
    print("FINAL 12-CATEGORY COMPARISON")
    print("=" * 60)
    print(comparison.to_string(index=False))
    print(f"\n{SPECIES_B} total unique SNVs: {reference_counts.sum():,}")
    print(f"{SPECIES_A} shared (tumor ∩ liquid biopsy) total SNVs: {species_a_counts.sum():,}")
    print(f"\nSaved: {comparison_path}")
    print(f"Saved: {figure_path}")


if __name__ == "__main__":
    main()
