"""
Batch-convert white/near-white image backgrounds to transparent.

Every supported image in the input folder is converted to RGBA and any
pixel with all three RGB channels at or above ``--white-threshold`` is
made fully transparent. Results are written to a subfolder so that the
original images are never overwritten.

Usage
-----
    python make_background_transparent.py --input-dir results/figures/raw_exports
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image

SUPPORTED_EXTENSIONS = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp"}
DEFAULT_WHITE_THRESHOLD = 245
DEFAULT_DPI = (300, 300)


def make_background_transparent(image_path: Path, output_path: Path, white_threshold: int, dpi: tuple[int, int]) -> None:
    """Convert near-white pixels in ``image_path`` to transparent and save as PNG."""
    image = Image.open(image_path).convert("RGBA")
    data = np.array(image)

    rgb = data[:, :, :3]
    white_mask = np.all(rgb >= white_threshold, axis=2)
    data[white_mask, 3] = 0

    Image.fromarray(data, "RGBA").save(output_path, format="PNG", dpi=dpi, optimize=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True, help="Folder containing images to process.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Defaults to <input-dir>/transparent.")
    parser.add_argument("--white-threshold", type=int, default=DEFAULT_WHITE_THRESHOLD, help="Minimum R,G,B value (0-255) treated as white.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir or (args.input_dir / "transparent")
    output_dir.mkdir(parents=True, exist_ok=True)

    processed = 0
    for image_path in sorted(args.input_dir.iterdir()):
        if image_path.suffix.lower() not in SUPPORTED_EXTENSIONS:
            continue

        output_path = output_dir / f"{image_path.stem}_transparent.png"
        print(f"Processing: {image_path.name}")
        make_background_transparent(image_path, output_path, args.white_threshold, DEFAULT_DPI)
        print(f"Saved: {output_path.name}")
        processed += 1

    print(f"\nDone. {processed} image(s) processed.")
    print(f"Output folder: {output_dir}")


if __name__ == "__main__":
    main()
