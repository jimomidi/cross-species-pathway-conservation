"""
Compare the functional (exonic-consequence) mutation-type spectrum between a
reference species-B annotation and species-A tumor/liquid-biopsy VCF pairs
annotated with ANNOVAR.

For each species-A sample, SNVs and indels shared between the tumor and
liquid-biopsy VCF files are matched against the corresponding ANNOVAR
``.exonic_variant_function`` files (SNVs by exact chromosome/position/
allele match; indels by a small positional window, to tolerate the
different left-alignment conventions used by VCF and ANNOVAR). Only shared
variants with an exonic annotation are counted; shared variants without one
are reported separately rather than being labeled "Other".

Input layout
------------
Annotation files are discovered from ``--annotation-dir`` using the naming
convention ``<sample>_Tumor.exonic_variant_function`` /
``<sample>_cfDNA.exonic_variant_function``, matching the VCF sample names
used by ``substitution_spectrum.py``.

Reference (species-B) counts
-----------------------------
This script does not derive species-B functional-type counts itself: bring
your own per-category counts (for example, computed from a public
annotated mutation catalog) as a two-column CSV via
``--reference-functional-types`` with columns ``Mutation_Type`` and
``Count``. See the note in the project README about this input.

Usage
-----
    python functional_mutation_types.py \\
        --annotation-dir data/annotation \\
        --vcf-dir data/variants \\
        --reference-functional-types data/reference_mutations/reference_functional_types.csv \\
        --output-dir results/functional_mutation_types
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import FONT_FAMILY, SPECIES_A, SPECIES_B  # noqa: E402
from vcf_utils import normalize_chromosome, read_vcf_variants  # noqa: E402

TUMOR_VCF_SUFFIX = "_Tumor.vcf.gz"
LIQUID_BIOPSY_VCF_SUFFIX = "_cfDNA.vcf.gz"
TUMOR_ANNOTATION_SUFFIX = "_Tumor.exonic_variant_function"
LIQUID_BIOPSY_ANNOTATION_SUFFIX = "_cfDNA.exonic_variant_function"

# Canonical functional-type categories produced by map_annovar_type().
MUTATION_TYPES = [
    "Missense substitution",
    "Synonymous substitution",
    "Nonsense substitution",
    "Frameshift insertion",
    "Frameshift deletion",
    "Inframe insertion",
    "Inframe deletion",
    "Complex mutation",
    "Other",
]

# Priority used by choose_best() to break ties when a variant matches more
# than one annotated functional type (see that function's docstring).
_TYPE_PRIORITY = {name: rank for rank, name in enumerate(MUTATION_TYPES)}


def map_annovar_type(raw_type: str) -> str:
    """Map an ANNOVAR exonic-function string to a common functional category."""
    text = str(raw_type).strip().lower()

    # "nonsynonymous" must be checked before "synonymous" (substring match).
    if "nonsynonymous snv" in text or "nonsynomymous snv" in text:
        return "Missense substitution"
    if "synonymous snv" in text:
        return "Synonymous substitution"
    if "stopgain" in text:
        return "Nonsense substitution"
    if "frameshift insertion" in text:
        return "Frameshift insertion"
    if "frameshift deletion" in text:
        return "Frameshift deletion"
    if "nonframeshift insertion" in text or "non-frameshift insertion" in text:
        return "Inframe insertion"
    if "nonframeshift deletion" in text or "non-frameshift deletion" in text:
        return "Inframe deletion"
    if "complex" in text or "delins" in text:
        return "Complex mutation"

    # "stoploss" and "unknown" are genuinely annotated types that fall
    # outside the main categories tracked here.
    return "Other"


def read_exonic_annotation(path: Path) -> pd.DataFrame:
    """Parse an ANNOVAR ``.exonic_variant_function`` file."""
    records = []

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            try:
                start, end = int(fields[4]), int(fields[5])
            except ValueError:
                continue

            records.append(
                {
                    "Chrom": normalize_chromosome(fields[3]),
                    "Start": start,
                    "End": end,
                    "Ref": str(fields[6]).upper(),
                    "Alt": str(fields[7]).upper(),
                    "Functional_Type": map_annovar_type(fields[1]),
                }
            )

    return pd.DataFrame(records)


def build_annotation_lookup(annotation: pd.DataFrame) -> tuple[dict, dict]:
    """Index an annotation table for exact (SNV) and positional (indel) matching."""
    exact_lookup: dict[tuple, list[str]] = {}
    position_lookup: dict[tuple, list[str]] = {}

    for _, row in annotation.iterrows():
        chrom, start, end = row["Chrom"], int(row["Start"]), int(row["End"])
        functional_type = row["Functional_Type"]

        exact_lookup.setdefault((chrom, start, row["Ref"], row["Alt"]), []).append(functional_type)

        # A small window tolerates the different anchor-base conventions
        # VCF and ANNOVAR use for insertions/deletions.
        for position in range(min(start, end) - 2, max(start, end) + 3):
            position_lookup.setdefault((chrom, position), []).append(functional_type)

    return exact_lookup, position_lookup


def choose_best(possible_types: list[str]) -> str:
    """Pick a single representative functional type from several candidates.

    A shared variant can match more than one annotated transcript record
    (and, for indels, more than one nearby position), each potentially
    reporting a different consequence. This resolves ties by:

    1. Taking the most frequently observed type among the candidates.
    2. Breaking remaining ties using the severity-ish order defined by
       ``MUTATION_TYPES`` (missense before synonymous, etc.), which
       matches the order the original classification was written in.

    Note: this tie-breaking function was referenced by the source project
    this script was adapted from, but was not defined anywhere in the
    provided code. The implementation here is a reasonable, explicit
    default -- confirm it matches your intended handling of multi-transcript
    annotations before relying on downstream counts.
    """
    counts = Counter(possible_types)
    max_count = max(counts.values())
    tied = [t for t, c in counts.items() if c == max_count]
    return min(tied, key=lambda t: _TYPE_PRIORITY.get(t, len(MUTATION_TYPES)))


def discover_sample_names(vcf_dir: Path) -> list[str]:
    return sorted(p.name[: -len(TUMOR_VCF_SUFFIX)] for p in vcf_dir.glob(f"*{TUMOR_VCF_SUFFIX}"))


def summarize_shared_functional_types(
    samples: list[str], vcf_dir: Path, annotation_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """For each sample, classify tumor/liquid-biopsy shared variants by functional type."""
    included_records = []
    excluded_records = []

    for sample in samples:
        tumor_vcf = vcf_dir / f"{sample}{TUMOR_VCF_SUFFIX}"
        liquid_biopsy_vcf = vcf_dir / f"{sample}{LIQUID_BIOPSY_VCF_SUFFIX}"
        tumor_annotation_path = annotation_dir / f"{sample}{TUMOR_ANNOTATION_SUFFIX}"
        liquid_biopsy_annotation_path = annotation_dir / f"{sample}{LIQUID_BIOPSY_ANNOTATION_SUFFIX}"

        missing = [
            p for p in [tumor_vcf, liquid_biopsy_vcf, tumor_annotation_path, liquid_biopsy_annotation_path] if not p.exists()
        ]
        if missing:
            print(f"Warning: sample '{sample}' is missing files, skipping: {missing}")
            continue

        shared = read_vcf_variants(tumor_vcf) & read_vcf_variants(liquid_biopsy_vcf)

        tumor_exact, tumor_pos = build_annotation_lookup(read_exonic_annotation(tumor_annotation_path))
        liquid_exact, liquid_pos = build_annotation_lookup(read_exonic_annotation(liquid_biopsy_annotation_path))

        included = excluded = 0

        for chrom, pos, ref, alt in shared:
            possible_types: list[str] = []

            if len(ref) == 1 and len(alt) == 1:
                key = (chrom, pos, ref, alt)
                possible_types.extend(tumor_exact.get(key, []))
                possible_types.extend(liquid_exact.get(key, []))
            else:
                for test_pos in range(pos - 2, pos + 3):
                    possible_types.extend(tumor_pos.get((chrom, test_pos), []))
                    possible_types.extend(liquid_pos.get((chrom, test_pos), []))

            if possible_types:
                included_records.append(
                    {
                        "Sample": sample,
                        "Chrom": chrom,
                        "Pos": pos,
                        "Ref": ref,
                        "Alt": alt,
                        "Functional_Type": choose_best(possible_types),
                    }
                )
                included += 1
            else:
                excluded_records.append(
                    {"Sample": sample, "Chrom": chrom, "Pos": pos, "Ref": ref, "Alt": alt, "Reason": "Shared but no exonic annotation"}
                )
                excluded += 1

        print(f"{sample}: shared={len(shared):,} | annotated & included={included:,} | excluded={excluded:,}")

    return pd.DataFrame(included_records), pd.DataFrame(excluded_records)


def plot_comparison(comparison: pd.DataFrame, output_path: Path) -> None:
    plot_df = comparison.set_index("Mutation_Type")[[f"{SPECIES_B}_Percent", f"{SPECIES_A}_Shared_Percent"]].copy()
    plot_df.columns = [SPECIES_B, SPECIES_A]

    plt.rcParams["font.family"] = FONT_FAMILY

    fig, ax = plt.subplots(figsize=(6, 5))
    plot_df.plot(kind="bar", ax=ax, width=0.76, color=["#20C9C3", "#7C5CFC"], edgecolor="none")

    ax.set_ylabel("Percentage of Mutations (%)", fontweight="bold", labelpad=8)
    ax.set_xticklabels(MUTATION_TYPES, rotation=45, ha="right", rotation_mode="anchor")
    ax.legend(title="", frameon=False, loc="upper center")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf-dir", type=Path, default=Path("data/variants"))
    parser.add_argument("--annotation-dir", type=Path, default=Path("data/annotation"))
    parser.add_argument(
        "--reference-functional-types",
        type=Path,
        required=True,
        help="CSV with 'Mutation_Type' and 'Count' columns for the reference species (see module docstring).",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/functional_mutation_types"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    reference = pd.read_csv(args.reference_functional_types)
    reference_counts = (
        reference.set_index("Mutation_Type")["Count"].reindex(MUTATION_TYPES, fill_value=0)
    )

    samples = discover_sample_names(args.vcf_dir)
    if not samples:
        raise FileNotFoundError(f"No '{TUMOR_VCF_SUFFIX}' files were found in {args.vcf_dir}")

    included_df, excluded_df = summarize_shared_functional_types(samples, args.vcf_dir, args.annotation_dir)

    species_a_counts = included_df["Functional_Type"].value_counts().reindex(MUTATION_TYPES, fill_value=0)

    comparison = pd.DataFrame(
        {
            "Mutation_Type": MUTATION_TYPES,
            f"{SPECIES_B}_Count": reference_counts.values,
            f"{SPECIES_B}_Percent": (reference_counts / reference_counts.sum() * 100).values if reference_counts.sum() else 0,
            f"{SPECIES_A}_Shared_Count": species_a_counts.values,
            f"{SPECIES_A}_Shared_Percent": (species_a_counts / species_a_counts.sum() * 100).values if species_a_counts.sum() else 0,
        }
    )

    comparison_path = args.output_dir / "functional_mutation_type_comparison.csv"
    comparison.to_csv(comparison_path, index=False)
    included_df.to_csv(args.output_dir / f"{SPECIES_A.lower()}_shared_exonic_functional_types.csv", index=False)
    excluded_df.to_csv(args.output_dir / f"{SPECIES_A.lower()}_shared_nonexonic_excluded.csv", index=False)

    figure_path = args.output_dir / "functional_mutation_type_comparison.png"
    plot_comparison(comparison, figure_path)

    print("\n" + "=" * 60)
    print("FINAL FUNCTIONAL TYPE COMPARISON")
    print("=" * 60)
    print(comparison.to_string(index=False))
    print(f"\nSaved: {comparison_path}")
    print(f"Saved: {figure_path}")


if __name__ == "__main__":
    main()
