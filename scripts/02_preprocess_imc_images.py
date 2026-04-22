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
    copy_images_csv(workflow_dir, args.numbered_images_copy)
    write_tiff_inventory(workflow_dir, args.numbered_inventory)

    print(f"Steinbock preprocessing complete: {workflow_dir / 'data' / 'img'}")
    print(f"Numbered images metadata: {workflow_dir / args.numbered_images_copy}")
    print(f"Numbered TIFF inventory: {workflow_dir / args.numbered_inventory}")
    print(f"Command log: {workflow_dir / args.log_file}")


if __name__ == "__main__":
    main()
