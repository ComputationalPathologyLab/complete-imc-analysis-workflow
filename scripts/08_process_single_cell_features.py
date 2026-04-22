#!/usr/bin/env python3
"""Process exported single-cell features according to the template paper.

Processing steps:

1. Filter out segmented objects with area smaller than the selected threshold.
2. Censor marker intensity columns to their selected upper percentile.
3. Apply arcsinh transformation to marker intensity columns with the selected
   cofactor.

The script keeps the raw exported `data/cells.csv` unchanged and writes a new
processed table for downstream analysis.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter, censor, and arcsinh-transform exported single-cell features."
    )
    parser.add_argument(
        "--workflow-dir",
        type=Path,
        default=Path.cwd(),
        help="Workflow directory containing data/cells.csv and data/panel.csv.",
    )
    parser.add_argument(
        "--input-cells",
        type=Path,
        default=Path("data/cells.csv"),
        help="Input exported Steinbock single-cell CSV, relative to workflow-dir.",
    )
    parser.add_argument(
        "--panel",
        type=Path,
        default=Path("data/panel.csv"),
        help="Panel CSV used to identify marker columns, relative to workflow-dir.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/11_processed_single_cell_features.csv"),
        help="Output processed single-cell CSV, relative to workflow-dir.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("results/11_processed_single_cell_features_summary.csv"),
        help="Processing summary CSV, relative to workflow-dir.",
    )
    parser.add_argument(
        "--area-column",
        default="area",
        help="Column containing segmented object area in pixels.",
    )
    parser.add_argument(
        "--min-area",
        type=float,
        default=4.0,
        help="Minimum retained object area in pixels. Objects below this value are removed.",
    )
    parser.add_argument(
        "--percentile",
        type=float,
        default=99.0,
        help="Upper percentile used for marker censoring.",
    )
    parser.add_argument(
        "--cofactor",
        type=float,
        default=1.0,
        help="Cofactor used for arcsinh transformation.",
    )
    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required {description} not found: {path}")


def read_panel_markers(panel_path: Path) -> list[str]:
    with panel_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        markers = [row["name"].strip() for row in reader if row.get("keep", "1") == "1"]

    markers = [marker for marker in markers if marker]
    if not markers:
        raise ValueError(f"No marker names were found in panel file: {panel_path}")
    return markers


def parse_float(value: str, column: str) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(f"Could not parse numeric value in column {column}: {value}") from exc


def calculate_percentile(values: list[float], percentile_value: float) -> float:
    if not values:
        raise ValueError("Cannot calculate percentile from an empty value list.")

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
    fraction = position - lower_index
    return lower_value + (upper_value - lower_value) * fraction


def validate_inputs(
    fieldnames: list[str],
    marker_columns: list[str],
    area_column: str,
    percentile: float,
    cofactor: float,
) -> None:
    missing_markers = [column for column in marker_columns if column not in fieldnames]
    if missing_markers:
        raise ValueError(
            "The following panel marker columns are missing from the exported cells table: "
            + ", ".join(missing_markers)
        )

    if area_column not in fieldnames:
        raise ValueError(f"Area column not found in exported cells table: {area_column}")

    if not 0 < percentile <= 100:
        raise ValueError("Percentile must be greater than 0 and less than or equal to 100.")

    if cofactor <= 0:
        raise ValueError("Cofactor must be greater than 0.")


def read_cells(cells_path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with cells_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        rows = list(reader)

    if fieldnames is None:
        raise ValueError(f"Could not read header from cells table: {cells_path}")

    return fieldnames, rows


def process_single_cell_features(
    cells: list[dict[str, str]],
    marker_columns: list[str],
    area_column: str,
    min_area: float,
    percentile: float,
    cofactor: float,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    input_cell_count = len(cells)

    filtered = [
        row.copy()
        for row in cells
        if parse_float(row[area_column], area_column) >= min_area
    ]
    filtered_cell_count = len(filtered)

    if not filtered:
        raise ValueError("No cells remain after area filtering.")

    marker_percentiles: dict[str, float] = {}
    for marker in marker_columns:
        values = [parse_float(row[marker], marker) for row in filtered]
        threshold = calculate_percentile(values, percentile)
        marker_percentiles[marker] = threshold
        for row in filtered:
            value = parse_float(row[marker], marker)
            censored = min(value, threshold)
            row[marker] = str(math.asinh(censored / cofactor))

    summary_rows = [
        {"metric": "input_cells", "value": input_cell_count},
        {"metric": "retained_cells_after_area_filter", "value": filtered_cell_count},
        {"metric": "removed_cells_area_lt_min_area", "value": input_cell_count - filtered_cell_count},
        {"metric": "min_area_pixels", "value": min_area},
        {"metric": "censoring_percentile", "value": percentile},
        {"metric": "arcsinh_cofactor", "value": cofactor},
        {"metric": "marker_columns_processed", "value": len(marker_columns)},
    ]
    for marker, threshold in marker_percentiles.items():
        summary_rows.append({"metric": f"{marker}_p{percentile:g}_censor_threshold", "value": threshold})

    return filtered, summary_rows


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()
    input_cells = workflow_dir / args.input_cells
    panel = workflow_dir / args.panel
    output = workflow_dir / args.output
    summary = workflow_dir / args.summary

    require_input(input_cells, "exported single-cell table")
    require_input(panel, "panel file")

    marker_columns = read_panel_markers(panel)
    fieldnames, cells = read_cells(input_cells)
    validate_inputs(fieldnames, marker_columns, args.area_column, args.percentile, args.cofactor)

    processed, processing_summary = process_single_cell_features(
        cells,
        marker_columns,
        args.area_column,
        args.min_area,
        args.percentile,
        args.cofactor,
    )

    write_rows(output, fieldnames, processed)
    write_rows(summary, ["metric", "value"], processing_summary)

    print(f"Input cells: {len(cells)}")
    print(f"Retained cells after area filtering: {len(processed)}")
    print(f"Removed cells with area < {args.min_area:g}: {len(cells) - len(processed)}")
    print(f"Marker columns processed: {len(marker_columns)}")
    print(f"Processed table: {output}")
    print(f"Processing summary: {summary}")


if __name__ == "__main__":
    main()
