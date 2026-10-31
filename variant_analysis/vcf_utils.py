"""
Shared helpers for reading gzip-compressed VCF files.

These utilities are intentionally minimal: they parse only the columns
needed by the mutation-spectrum scripts in this package (chromosome,
position, reference allele, alternate allele) and do not depend on a full
VCF-parsing library, so they have no third-party dependency beyond the
standard library.
"""

from __future__ import annotations

import gzip
from pathlib import Path

Variant = tuple[str, int, str, str]


def normalize_chromosome(chrom: str) -> str:
    """Strip a leading "chr" prefix and lowercase the chromosome name."""
    chrom = str(chrom).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom.lower()


def read_vcf_variants(path: Path) -> set[Variant]:
    """Read a gzip-compressed VCF file and return the set of called variants.

    Multiallelic records are expanded into one entry per alternate allele.
    Only positions with a well-formed REF/ALT pair are kept; records where
    REF or ALT is missing (``.``) or a spanning-deletion placeholder (``*``)
    are skipped. Chromosome names are normalized with
    :func:`normalize_chromosome` so that "chr1" and "1" are treated as the
    same chromosome.

    Returns a set of ``(chromosome, position, ref, alt)`` tuples.
    """
    variants: set[Variant] = set()

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue

            chrom = normalize_chromosome(fields[0])

            try:
                position = int(fields[1])
            except ValueError:
                continue

            ref = fields[3].upper()
            alts = fields[4].upper().split(",")

            for alt in alts:
                if ref in (".", "") or alt in (".", "", "*") or ref == alt:
                    continue
                variants.add((chrom, position, ref, alt))

    return variants


def filter_snvs(variants: set[Variant]) -> set[Variant]:
    """Restrict a variant set to single-nucleotide substitutions."""
    bases = set("ACGT")
    return {
        (chrom, pos, ref, alt)
        for chrom, pos, ref, alt in variants
        if len(ref) == 1 and len(alt) == 1 and ref in bases and alt in bases
    }
