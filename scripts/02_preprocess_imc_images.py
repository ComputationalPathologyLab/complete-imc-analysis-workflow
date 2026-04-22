#!/usr/bin/env python3
"""Run Steinbock IMC preprocessing for a workflow data directory.

This script is intended to be called from the workflow notebook. It follows the
Steinbock hands-on preprocessing command:

    steinbock preprocess imc images --hpf 50

The script uses the official Steinbock Docker image by default. Steinbock
working outputs keep the required names inside `data/`, while numbered research
copies and logs are written to `results/` and `logs/`.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Preprocess raw IMC data with Steinbock.")
    parser.add_argument(
        "--workflow-dir",
        type=Path,
        default=Path.cwd(),
        help="Workflow directory containing data/raw and data/panel.csv.",
    )
    parser.add_argument(
        "--steinbock-image",
        default="ghcr.io/bodenmillergroup/steinbock:0.16.1",
        help="Steinbock Docker image.",
    )
    parser.add_argument(
        "--hpf",
        default="50",
        help="Hot-pixel filtering threshold passed to Steinbock.",
    )
    parser.add_argument(
        "--numbered-images-copy",
        type=Path,
        default=Path("results/02_images.csv"),
        help="Numbered copy of Steinbock data/images.csv.",
    )
    parser.add_argument(
        "--numbered-inventory",
        type=Path,
        default=Path("results/03_preprocessed_tiff_inventory.csv"),
        help="Numbered inventory of generated TIFF image files.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("logs/02_preprocess_imc_images.log"),
        help="Numbered preprocessing command log.",
    )
    parser.add_argument(
        "--keep-source",
        choices=["mcd", "all"],
        default="mcd",
        help=(
            "Images to retain after Steinbock preprocessing. The template study "
            "uses MCD-derived acquisitions; ROI text files are retained in raw/ "
            "for metadata/header recovery but should not become separate images."
        ),
    )
    parser.add_argument(
        "--min-width",
        type=int,
        default=1000,
        help="Minimum retained image width in pixels. Use 0 to disable.",
    )
    parser.add_argument(
        "--min-height",
        type=int,
        default=1000,
        help="Minimum retained image height in pixels. Use 0 to disable.",
    )
    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required {description} not found: {path}")


def run_steinbock(workflow_dir: Path, image: str, hpf: str, log_file: Path) -> None:
    data_dir = workflow_dir / "data"
    command = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{data_dir.resolve()}:/data",
        "-u",
        f"{subprocess.check_output(['id', '-u'], text=True).strip()}:"
        f"{subprocess.check_output(['id', '-g'], text=True).strip()}",
        image,
        "preprocess",
        "imc",
        "images",
        "--hpf",
        hpf,
    ]

    log_file.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(command, cwd=workflow_dir, text=True, capture_output=True)

    log_file.write_text(
        "COMMAND\n"
        + " ".join(command)
        + "\n\nSTDOUT\n"
        + completed.stdout
        + "\n\nSTDERR\n"
        + completed.stderr
        + f"\n\nRETURN_CODE\n{completed.returncode}\n",
        encoding="utf-8",
    )

    if completed.returncode != 0:
        raise RuntimeError(f"Steinbock preprocessing failed. See log: {log_file}")


def copy_images_csv(workflow_dir: Path, numbered_images_copy: Path) -> None:
    source = workflow_dir / "data" / "images.csv"
    destination = workflow_dir / numbered_images_copy
    require_input(source, "Steinbock images.csv output")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def filter_preprocessed_images(
    workflow_dir: Path,
    keep_source: str,
    min_width: int,
    min_height: int,
    log_file: Path,
) -> None:
    images_csv = workflow_dir / "data" / "images.csv"
    img_dir = workflow_dir / "data" / "img"
    require_input(images_csv, "Steinbock images.csv output")
    require_input(img_dir, "Steinbock image output directory")

    with images_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames

    if fieldnames is None:
        raise ValueError(f"Could not read header from {images_csv}")

    kept_rows = []
    removed_rows = []

    for row in rows:
        remove_reasons = []
        source_file = row.get("source_file", "")
        image_name = row.get("image", "")
        width = int(float(row.get("width_px") or 0))
        height = int(float(row.get("height_px") or 0))

        if keep_source == "mcd" and not source_file.endswith(".mcd"):
            remove_reasons.append("source_not_mcd")
        if min_width and width < min_width:
            remove_reasons.append(f"width_lt_{min_width}")
        if min_height and height < min_height:
            remove_reasons.append(f"height_lt_{min_height}")

        if remove_reasons:
            removed_rows.append((row, remove_reasons))
            image_path = img_dir / image_name
            if image_path.exists():
                image_path.unlink()
        else:
            kept_rows.append(row)

    with images_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(kept_rows)

    summary_lines = [
        "",
        "POSTPROCESSING_FILTER",
        f"keep_source={keep_source}",
        f"min_width={min_width}",
        f"min_height={min_height}",
        f"images_kept={len(kept_rows)}",
        f"images_removed={len(removed_rows)}",
    ]
    for row, reasons in removed_rows:
        summary_lines.append(
            f"removed={row.get('image', '')}; source={row.get('source_file', '')}; "
            f"size={row.get('width_px', '')}x{row.get('height_px', '')}; "
            f"reason={'+'.join(reasons)}"
        )

    with log_file.open("a", encoding="utf-8") as handle:
        handle.write("\n".join(summary_lines) + "\n")


def write_tiff_inventory(workflow_dir: Path, numbered_inventory: Path) -> None:
    img_dir = workflow_dir / "data" / "img"
    require_input(img_dir, "Steinbock image output directory")

    rows = []
    for index, path in enumerate(sorted(img_dir.glob("*.tiff")), start=1):
        rows.append(
            {
                "index": index,
                "file_name": path.name,
                "relative_path": str(path.relative_to(workflow_dir)),
                "size_bytes": path.stat().st_size,
            }
        )

    destination = workflow_dir / numbered_inventory
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["index", "file_name", "relative_path", "size_bytes"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()

    require_input(workflow_dir / "data" / "raw", "raw data directory")
    require_input(workflow_dir / "data" / "panel.csv", "Steinbock panel file")

    run_steinbock(workflow_dir, args.steinbock_image, args.hpf, workflow_dir / args.log_file)
    filter_preprocessed_images(
        workflow_dir,
        args.keep_source,
        args.min_width,
        args.min_height,
        workflow_dir / args.log_file,
    )
    copy_images_csv(workflow_dir, args.numbered_images_copy)
    write_tiff_inventory(workflow_dir, args.numbered_inventory)

    print(f"Steinbock preprocessing complete: {workflow_dir / 'data' / 'img'}")
    print(f"Numbered images metadata: {workflow_dir / args.numbered_images_copy}")
    print(f"Numbered TIFF inventory: {workflow_dir / args.numbered_inventory}")
    print(f"Command log: {workflow_dir / args.log_file}")


if __name__ == "__main__":
    main()
