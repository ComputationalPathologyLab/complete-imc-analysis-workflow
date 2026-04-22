#!/usr/bin/env python3
"""Run Steinbock Mesmer segmentation for preprocessed IMC images.

The script is intended to be called from the workflow notebook. It follows the
study-style Steinbock segmentation command:

    steinbock segment deepcell --app mesmer --minmax

Steinbock working masks are written to `data/masks/`. A numbered mask inventory
and command log are written to `results/` and `logs/`.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Segment IMC images with Steinbock Mesmer.")
    parser.add_argument(
        "--workflow-dir",
        type=Path,
        default=Path.cwd(),
        help="Workflow directory containing the data folder.",
    )
    parser.add_argument(
        "--steinbock-image",
        default="ghcr.io/bodenmillergroup/steinbock:0.16.1",
        help="Steinbock Docker image.",
    )
    parser.add_argument(
        "--token-env",
        default="DEEPCELL_ACCESS_TOKEN",
        help="Environment variable name containing the DeepCell access token.",
    )
    parser.add_argument(
        "--numbered-inventory",
        type=Path,
        default=Path("results/04_mesmer_mask_inventory.csv"),
        help="Numbered inventory of generated Mesmer mask files.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("logs/03_segment_mesmer.log"),
        help="Numbered Mesmer segmentation command log.",
    )
    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required {description} not found: {path}")


def run_steinbock_mesmer(workflow_dir: Path, image: str, token_env: str, log_file: Path) -> None:
    data_dir = workflow_dir / "data"
    user_group = (
        f"{subprocess.check_output(['id', '-u'], text=True).strip()}:"
        f"{subprocess.check_output(['id', '-g'], text=True).strip()}"
    )

    command = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{data_dir.resolve()}:/data",
        "-u",
        user_group,
    ]

    if os.environ.get(token_env):
        command.extend(["-e", token_env])

    command.extend(
        [
            image,
            "segment",
            "deepcell",
            "--app",
            "mesmer",
            "--minmax",
        ]
    )

    log_file.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(command, cwd=workflow_dir, text=True, capture_output=True)

    token_status = "present" if os.environ.get(token_env) else "not_present"
    log_file.write_text(
        "COMMAND\n"
        + " ".join(command)
        + f"\n\nTOKEN_ENVIRONMENT_VARIABLE\n{token_env}={token_status}\n"
        + "\nSTDOUT\n"
        + completed.stdout
        + "\n\nSTDERR\n"
        + completed.stderr
        + f"\n\nRETURN_CODE\n{completed.returncode}\n",
        encoding="utf-8",
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "Steinbock Mesmer segmentation failed. "
            f"This commonly occurs when DeepCell model access is unavailable. See log: {log_file}"
        )


def write_mask_inventory(workflow_dir: Path, numbered_inventory: Path) -> None:
    masks_dir = workflow_dir / "data" / "masks"
    require_input(masks_dir, "Steinbock masks output directory")

    rows = []
    for index, path in enumerate(sorted(masks_dir.glob("*.tiff")), start=1):
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

    require_input(workflow_dir / "data" / "img", "preprocessed image directory")
    require_input(workflow_dir / "data" / "images.csv", "Steinbock images.csv file")
    require_input(workflow_dir / "data" / "panel.csv", "Steinbock panel file")

    run_steinbock_mesmer(workflow_dir, args.steinbock_image, args.token_env, workflow_dir / args.log_file)
    write_mask_inventory(workflow_dir, args.numbered_inventory)

    print(f"Mesmer segmentation complete: {workflow_dir / 'data' / 'masks'}")
    print(f"Numbered mask inventory: {workflow_dir / args.numbered_inventory}")
    print(f"Command log: {workflow_dir / args.log_file}")


if __name__ == "__main__":
    main()
