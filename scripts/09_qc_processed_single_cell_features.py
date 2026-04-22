#!/usr/bin/env python3
"""Create QC summaries and figures for processed single-cell features."""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import Counter
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="QC processed single-cell features before downstream analysis."
    )
    parser.add_argument(
        "--workflow-dir",
        type=Path,
        default=Path.cwd(),
        help="Workflow directory containing results and figures folders.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("results/11_processed_single_cell_features.csv"),
        help="Processed single-cell feature table, relative to workflow-dir.",
    )
    parser.add_argument(
        "--panel",
        type=Path,
        default=Path("data/panel.csv"),
        help="Panel CSV used to identify marker columns, relative to workflow-dir.",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=Path("results/12_processed_single_cell_qc_summary.csv"),
        help="Overall QC summary CSV, relative to workflow-dir.",
    )
    parser.add_argument(
        "--marker-output",
        type=Path,
        default=Path("results/12_processed_marker_qc_summary.csv"),
        help="Marker-level QC summary CSV, relative to workflow-dir.",
    )
    parser.add_argument(
        "--cell-count-figure",
        type=Path,
        default=Path("figures/01_cells_per_image.png"),
        help="Cell count figure path, relative to workflow-dir.",
    )
    parser.add_argument(
        "--area-figure",
        type=Path,
        default=Path("figures/02_area_distribution_after_filtering.png"),
        help="Area distribution figure path, relative to workflow-dir.",
    )
    parser.add_argument(
        "--marker-figure",
        type=Path,
        default=Path("figures/03_marker_summary_after_transformation.png"),
        help="Marker summary figure path, relative to workflow-dir.",
    )
    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required {description} not found: {path}")


def parse_float(value: str, column: str) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(f"Could not parse numeric value in column {column}: {value}") from exc


def percentile(values: list[float], percentile_value: float) -> float:
    if not values:
        return math.nan
    sorted_values = sorted(values)
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = (percentile_value / 100.0) * (len(sorted_values) - 1)
    lower_index = math.floor(position)
    upper_index = math.ceil(position)
    if lower_index == upper_index:
        return sorted_values[int(position)]
    lower_value = sorted_values[lower_index]
    upper_value = sorted_values[upper_index]
    return lower_value + (upper_value - lower_value) * (position - lower_index)


def read_marker_columns(panel_path: Path) -> list[str]:
    with panel_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        markers = [row["name"].strip() for row in reader if row.get("keep", "1") == "1"]
    return [marker for marker in markers if marker]


def read_processed_table(
    input_path: Path,
    marker_columns: list[str],
) -> tuple[list[str], Counter[str], list[float], dict[str, list[float]], dict[str, int]]:
    image_counts: Counter[str] = Counter()
    areas: list[float] = []
    marker_values = {marker: [] for marker in marker_columns}
    missing_or_invalid = {marker: 0 for marker in marker_columns}

    with input_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        missing_markers = [marker for marker in marker_columns if marker not in fieldnames]
        if missing_markers:
            raise ValueError("Missing marker columns: " + ", ".join(missing_markers))
        if "Image" not in fieldnames or "area" not in fieldnames:
            raise ValueError("Processed table must contain Image and area columns.")

        for row in reader:
            image_counts[row["Image"]] += 1
            areas.append(parse_float(row["area"], "area"))
            for marker in marker_columns:
                try:
                    value = float(row[marker])
                except ValueError:
                    missing_or_invalid[marker] += 1
                    continue
                if math.isnan(value) or math.isinf(value):
                    missing_or_invalid[marker] += 1
                    continue
                marker_values[marker].append(value)

    return fieldnames, image_counts, areas, marker_values, missing_or_invalid


def write_overall_summary(
    path: Path,
    fieldnames: list[str],
    image_counts: Counter[str],
    areas: list[float],
    marker_columns: list[str],
) -> None:
    rows = [
        {"metric": "total_cells", "value": sum(image_counts.values())},
        {"metric": "image_count", "value": len(image_counts)},
        {"metric": "table_columns", "value": len(fieldnames)},
        {"metric": "marker_columns", "value": len(marker_columns)},
        {"metric": "minimum_area", "value": min(areas) if areas else math.nan},
        {"metric": "median_area", "value": percentile(areas, 50)},
        {"metric": "area_p99", "value": percentile(areas, 99)},
        {"metric": "maximum_area", "value": max(areas) if areas else math.nan},
    ]
    for image, count in sorted(image_counts.items()):
        rows.append({"metric": f"cells_in_{image}", "value": count})

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
        writer.writeheader()
        writer.writerows(rows)


def write_marker_summary(
    path: Path,
    marker_values: dict[str, list[float]],
    missing_or_invalid: dict[str, int],
) -> list[dict[str, float | str | int]]:
    rows = []
    for marker, values in marker_values.items():
        rows.append(
            {
                "marker": marker,
                "valid_values": len(values),
                "missing_or_invalid_values": missing_or_invalid[marker],
                "minimum": min(values) if values else math.nan,
                "median": percentile(values, 50),
                "p99": percentile(values, 99),
                "maximum": max(values) if values else math.nan,
            }
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "marker",
                "valid_values",
                "missing_or_invalid_values",
                "minimum",
                "median",
                "p99",
                "maximum",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows


def save_cell_count_figure(path: Path, image_counts: Counter[str]) -> None:
    labels = list(sorted(image_counts))
    values = [image_counts[label] for label in labels]
    path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 4.8))
    plt.bar(labels, values, color="#3B6EA8")
    plt.ylabel("Cells")
    plt.xlabel("Image")
    plt.title("Cells Per Image")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def save_area_figure(path: Path, areas: list[float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(7.5, 4.8))
    plt.hist(areas, bins=60, color="#5D8A66", edgecolor="white")
    plt.xlabel("Area (pixels)")
    plt.ylabel("Cells")
    plt.title("Area Distribution After Filtering")
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def save_marker_figure(path: Path, marker_rows: list[dict[str, float | str | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    markers = [str(row["marker"]) for row in marker_rows]
    medians = [float(row["median"]) for row in marker_rows]
    p99_values = [float(row["p99"]) for row in marker_rows]

    plt.figure(figsize=(12, 5.2))
    x_positions = range(len(markers))
    plt.plot(x_positions, medians, marker="o", linewidth=1.5, label="Median", color="#3B6EA8")
    plt.plot(x_positions, p99_values, marker="o", linewidth=1.5, label="P99", color="#B85C38")
    plt.ylabel("Transformed Marker Value")
    plt.xlabel("Marker")
    plt.title("Marker Summary After Censoring and Arcsinh Transformation")
    plt.xticks(list(x_positions), markers, rotation=70, ha="right")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()
    input_path = workflow_dir / args.input
    panel_path = workflow_dir / args.panel

    require_input(input_path, "processed single-cell feature table")
    require_input(panel_path, "panel file")

    marker_columns = read_marker_columns(panel_path)
    fieldnames, image_counts, areas, marker_values, missing_or_invalid = read_processed_table(
        input_path,
        marker_columns,
    )
    marker_rows = write_marker_summary(workflow_dir / args.marker_output, marker_values, missing_or_invalid)
    write_overall_summary(workflow_dir / args.summary_output, fieldnames, image_counts, areas, marker_columns)
    save_cell_count_figure(workflow_dir / args.cell_count_figure, image_counts)
    save_area_figure(workflow_dir / args.area_figure, areas)
    save_marker_figure(workflow_dir / args.marker_figure, marker_rows)

    print(f"Total cells: {sum(image_counts.values())}")
    print(f"Images: {len(image_counts)}")
    print(f"Marker columns checked: {len(marker_columns)}")
    print(f"Overall QC summary: {workflow_dir / args.summary_output}")
    print(f"Marker QC summary: {workflow_dir / args.marker_output}")
    print(f"Cell count figure: {workflow_dir / args.cell_count_figure}")
    print(f"Area figure: {workflow_dir / args.area_figure}")
    print(f"Marker figure: {workflow_dir / args.marker_figure}")


if __name__ == "__main__":
    main()
