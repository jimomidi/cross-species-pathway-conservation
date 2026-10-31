"""
Convert a raw Ensembl BioMart export into a clean CSV file.

BioMart's "Export all results to" file (typically named ``mart_export.txt``)
is usually tab-separated, but the delimiter can vary depending on how it was
downloaded. This script detects the delimiter automatically, removes exact
duplicate rows, and writes a standard UTF-8 CSV file that downstream scripts
(see ``map_orthologs.py``) can read reliably.

Usage
-----
    python convert_biomart_export.py --input data/biomart/mart_export.txt

Input
-----
A BioMart text export with a header row (comma-, tab-, or semicolon-
separated).

Output
------
A CSV file with the same base name and a ``.csv`` extension, written next to
the input file (or to ``--output`` if provided).
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


def detect_delimiter(file_path: Path, sample_size: int = 10_000) -> str:
    """Sniff the field delimiter from the first ``sample_size`` characters.

    Falls back to a tab character if the delimiter cannot be determined,
    since BioMart's default export format is tab-separated.
    """
    with open(file_path, "r", encoding="utf-8-sig", errors="replace") as handle:
        sample = handle.read(sample_size)

    try:
        return csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";"]).delimiter
    except csv.Error:
        return "\t"


def convert_biomart_export(input_path: Path, output_path: Path) -> pd.DataFrame:
    """Read a raw BioMart export and write it out as a clean CSV file."""
    delimiter = detect_delimiter(input_path)
    print(f"Detected delimiter: {delimiter!r}")

    table = pd.read_csv(
        input_path,
        sep=delimiter,
        dtype=str,
        encoding="utf-8-sig",
        engine="python",
    )

    table = table.drop_duplicates()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False, encoding="utf-8-sig")

    print(f"\nCSV file saved to:\n{output_path}")
    print(f"\nNumber of rows: {table.shape[0]}")
    print(f"Number of columns: {table.shape[1]}")

    print("\nColumn names:")
    for i, column in enumerate(table.columns, start=1):
        print(f"{i}. {column}")

    return table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/biomart/mart_export.txt"),
        help="Path to the raw BioMart export (default: %(default)s).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Path to the output CSV file (default: same name as --input, "
        "with a .csv extension).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_path = args.output or args.input.with_suffix(".csv")

    if not args.input.exists():
        raise FileNotFoundError(f"Input file was not found:\n{args.input}")

    table = convert_biomart_export(args.input, output_path)
    print(f"\nFirst 10 rows:\n{table.head(10)}")


if __name__ == "__main__":
    main()
