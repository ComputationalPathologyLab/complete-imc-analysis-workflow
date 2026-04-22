#!/usr/bin/env python3
"""Validate expected outputs from the Steinbock reproduction workflow.

This script is intended to be called from the workflow notebook after all
Steinbock preprocessing, segmentation, measurement, neighbor, and export steps
have completed.

The script does not create scientific results. It creates a numbered validation
table documenting whether expected workflow outputs are present.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate expected workflow outputs.")
    parser.add_argument(
        "--workflow-dir",
        type=Path,
        default=Path.cwd(),
        help="Workflow directory containing data, logs, and results folders.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/11_workflow_output_validation.csv"),
        help="Numbered validation table path.",
    )
    return parser.parse_args()


def count_files(path: Path, pattern: str) -> int:
    if not path.exists():
        return 0
    return len(list(path.glob(pattern)))


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()

    checks = [
        ("raw_mcd_files", "data/raw", "*.mcd", None),
        ("raw_roi_text_files", "data/raw", "*.txt", None),
        ("working_panel_csv", "data/panel.csv", None, 1),
        ("numbered_panel_csv", "results/01_panel.csv", None, 1),
        ("preprocessed_tiff_images", "data/img", "*.tiff", None),
        ("working_images_csv", "data/images.csv", None, 1),
        ("numbered_images_csv", "results/02_images.csv", None, 1),
        ("preprocessed_tiff_inventory", "results/03_preprocessed_tiff_inventory.csv", None, 1),
        ("mesmer_masks", "data/masks", "*.tiff", None),
        ("mesmer_mask_inventory", "results/04_mesmer_mask_inventory.csv", None, 1),
        ("intensity_tables", "data/intensities", "*.csv", None),
        ("intensity_inventory", "results/05_intensity_table_inventory.csv", None, 1),
        ("regionprops_tables", "data/regionprops", "*.csv", None),
        ("regionprops_inventory", "results/06_regionprops_table_inventory.csv", None, 1),
        ("neighbor_tables", "data/neighbors", "*.csv", None),
        ("neighbor_inventory", "results/07_neighbors_table_inventory.csv", None, 1),
        ("working_cells_csv", "data/cells.csv", None, 1),
        ("numbered_cells_csv", "results/08_cells.csv", None, 1),
        ("working_cells_h5ad", "data/cells.h5ad", None, 1),
        ("numbered_cells_h5ad", "results/09_cells.h5ad", None, 1),
        ("working_graphml_files", "data/graphs", "*.graphml", None),
        ("numbered_graphml_files", "results/10_graphml", "*.graphml", None),
        ("numbered_graphml_inventory", "results/10_graphml_inventory.csv", None, 1),
    ]

    rows = []
    for check_name, relative_path, pattern, expected_count in checks:
        path = workflow_dir / relative_path
        if pattern is None:
            observed_count = 1 if path.exists() else 0
        else:
            observed_count = count_files(path, pattern)

        if expected_count is None:
            passed = observed_count > 0
            expected = ">0"
        else:
            passed = observed_count == expected_count
            expected = str(expected_count)

        rows.append(
            {
                "check_name": check_name,
                "relative_path": relative_path,
                "pattern": pattern or "",
                "expected_count": expected,
                "observed_count": observed_count,
                "passed": passed,
            }
        )

    output = workflow_dir / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "check_name",
                "relative_path",
                "pattern",
                "expected_count",
                "observed_count",
                "passed",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    passed_count = sum(row["passed"] for row in rows)
    print(f"Validation checks passed: {passed_count}/{len(rows)}")
    print(f"Validation table: {output}")


if __name__ == "__main__":
    main()
