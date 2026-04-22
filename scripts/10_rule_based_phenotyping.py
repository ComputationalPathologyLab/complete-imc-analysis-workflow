#!/usr/bin/env python3
"""Assign simplified rule-based phenotypes to processed single-cell features."""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

Path(".matplotlib").mkdir(exist_ok=True)
Path(".cache").mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib").resolve()))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(".cache").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


LINEAGE_MARKERS = [
    "CD138",
    "CD38",
    "IRF4",
    "CD3",
    "CD4",
    "CD8",
    "CD45",
    "CD68",
    "CD11b",
    "CD11c",
    "HLA-DR",
    "MPO",
    "CD34",
    "Vimentin",
    "CollagenTypeI",
    "Perilipin",
    "CathepsinK",
    "RUNX2",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assign transparent rule-based phenotypes.")
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
        help="Processed single-cell table, relative to workflow-dir.",
    )
    parser.add_argument(
        "--output-cells",
        type=Path,
        default=Path("results/13_cells_with_phenotypes.csv"),
        help="Output table with phenotype labels, relative to workflow-dir.",
    )
    parser.add_argument(
        "--composition-image",
        type=Path,
        default=Path("results/14_phenotype_composition_by_image.csv"),
        help="Phenotype composition by image, relative to workflow-dir.",
    )
    parser.add_argument(
        "--composition-category",
        type=Path,
        default=Path("results/15_phenotype_composition_by_category.csv"),
        help="Phenotype composition by category, relative to workflow-dir.",
    )
    parser.add_argument(
        "--thresholds-output",
        type=Path,
        default=Path("results/16_rule_based_phenotyping_thresholds.csv"),
        help="Marker positivity thresholds, relative to workflow-dir.",
    )
    parser.add_argument(
        "--counts-figure",
        type=Path,
        default=Path("figures/04_phenotype_counts_by_image.png"),
        help="Phenotype count figure path, relative to workflow-dir.",
    )
    parser.add_argument(
        "--category-figure",
        type=Path,
        default=Path("figures/05_phenotype_composition_by_category.png"),
        help="Phenotype composition figure path, relative to workflow-dir.",
    )
    parser.add_argument(
        "--positive-percentile",
        type=float,
        default=80.0,
        help="Percentile used to define marker positivity from processed marker values.",
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
    sorted_values = sorted(values)
    if not sorted_values:
        return math.nan
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = (percentile_value / 100.0) * (len(sorted_values) - 1)
    lower_index = math.floor(position)
    upper_index = math.ceil(position)
    if lower_index == upper_index:
        return sorted_values[int(position)]
    lower = sorted_values[lower_index]
    upper = sorted_values[upper_index]
    return lower + (upper - lower) * (position - lower_index)


def read_cells(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        rows = list(reader)
    if fieldnames is None:
        raise ValueError(f"Could not read header from {path}")
    return fieldnames, rows


def calculate_thresholds(
    rows: list[dict[str, str]],
    markers: list[str],
    positive_percentile: float,
) -> dict[str, float]:
    thresholds = {}
    for marker in markers:
        values = [parse_float(row[marker], marker) for row in rows]
        thresholds[marker] = percentile(values, positive_percentile)
    return thresholds


def is_positive(row: dict[str, str], marker: str, thresholds: dict[str, float]) -> bool:
    return parse_float(row[marker], marker) >= thresholds[marker]


def category_from_image(image_name: str) -> str:
    if "_MGUS" in image_name:
        return "MGUS"
    if "_SMM" in image_name:
        return "SMM"
    if "_UB" in image_name:
        return "UB"
    if "_B" in image_name:
        return "B"
    return "Unknown"


def assign_phenotype(row: dict[str, str], thresholds: dict[str, float]) -> str:
    positive = lambda marker: is_positive(row, marker, thresholds)

    if positive("CD138") and (positive("CD38") or positive("IRF4")):
        return "Plasma_cell_like"
    if positive("CD3") and positive("CD4") and not positive("CD8"):
        return "CD4_T_cell_like"
    if positive("CD3") and positive("CD8") and not positive("CD4"):
        return "CD8_T_cell_like"
    if positive("CD3"):
        return "T_cell_like"
    if positive("MPO") and positive("CD11b"):
        return "Neutrophil_like"
    if positive("CD68") and (positive("CD11b") or positive("CD45")):
        return "Macrophage_monocyte_like"
    if positive("HLA-DR") and positive("CD11c"):
        return "Dendritic_myeloid_like"
    if positive("CD45") and (positive("CD68") or positive("CD11b") or positive("HLA-DR")):
        return "Myeloid_cell_like"
    if positive("CD34"):
        return "CD34_positive_cell_like"
    if positive("Perilipin"):
        return "Adipocyte_like"
    if positive("CathepsinK"):
        return "Osteoclast_like"
    if positive("RUNX2") and positive("CollagenTypeI"):
        return "Osteoblast_stromal_like"
    if positive("Vimentin") and positive("CollagenTypeI"):
        return "Stromal_like"
    return "Unknown"


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_composition_rows(
    rows: list[dict[str, str]],
    group_column: str,
    phenotype_column: str = "phenotype",
) -> list[dict[str, object]]:
    group_totals = Counter(row[group_column] for row in rows)
    grouped_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        grouped_counts[row[group_column]][row[phenotype_column]] += 1

    output_rows = []
    for group in sorted(grouped_counts):
        for phenotype, count in sorted(grouped_counts[group].items()):
            total = group_totals[group]
            output_rows.append(
                {
                    group_column: group,
                    "phenotype": phenotype,
                    "cell_count": count,
                    "fraction": count / total if total else 0,
                }
            )
    return output_rows


def write_thresholds(path: Path, thresholds: dict[str, float], positive_percentile: float) -> None:
    rows = [
        {
            "marker": marker,
            "positive_percentile": positive_percentile,
            "positive_threshold": threshold,
        }
        for marker, threshold in sorted(thresholds.items())
    ]
    write_csv(path, ["marker", "positive_percentile", "positive_threshold"], rows)


def plot_counts_by_image(path: Path, rows: list[dict[str, str]]) -> None:
    image_counts: dict[str, Counter[str]] = defaultdict(Counter)
    phenotypes = sorted({row["phenotype"] for row in rows})
    images = sorted({row["Image"] for row in rows})
    for row in rows:
        image_counts[row["Image"]][row["phenotype"]] += 1

    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(11, 5.5))
    bottoms = [0] * len(images)
    for phenotype in phenotypes:
        values = [image_counts[image][phenotype] for image in images]
        plt.bar(images, values, bottom=bottoms, label=phenotype)
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
    plt.ylabel("Cells")
    plt.xlabel("Image")
    plt.title("Rule-Based Phenotype Counts By Image")
    plt.xticks(rotation=45, ha="right")
    plt.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def plot_composition_by_category(path: Path, composition_rows: list[dict[str, object]]) -> None:
    categories = sorted({str(row["category"]) for row in composition_rows})
    phenotypes = sorted({str(row["phenotype"]) for row in composition_rows})
    fractions = {
        (str(row["category"]), str(row["phenotype"])): float(row["fraction"])
        for row in composition_rows
    }

    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8.5, 5.2))
    bottoms = [0.0] * len(categories)
    for phenotype in phenotypes:
        values = [fractions.get((category, phenotype), 0.0) for category in categories]
        plt.bar(categories, values, bottom=bottoms, label=phenotype)
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
    plt.ylabel("Fraction of Cells")
    plt.xlabel("Category")
    plt.title("Rule-Based Phenotype Composition By Category")
    plt.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()
    input_path = workflow_dir / args.input
    require_input(input_path, "processed single-cell feature table")

    fieldnames, rows = read_cells(input_path)
    missing_markers = [marker for marker in LINEAGE_MARKERS if marker not in fieldnames]
    if missing_markers:
        raise ValueError("Missing required lineage markers: " + ", ".join(missing_markers))

    thresholds = calculate_thresholds(rows, LINEAGE_MARKERS, args.positive_percentile)
    output_rows: list[dict[str, str]] = []
    for row in rows:
        output_row = row.copy()
        output_row["category"] = category_from_image(row["Image"])
        output_row["phenotype"] = assign_phenotype(row, thresholds)
        output_rows.append(output_row)

    output_fieldnames = fieldnames + ["category", "phenotype"]
    write_csv(workflow_dir / args.output_cells, output_fieldnames, output_rows)

    composition_by_image = build_composition_rows(output_rows, "Image")
    write_csv(
        workflow_dir / args.composition_image,
        ["Image", "phenotype", "cell_count", "fraction"],
        composition_by_image,
    )

    composition_by_category = build_composition_rows(output_rows, "category")
    write_csv(
        workflow_dir / args.composition_category,
        ["category", "phenotype", "cell_count", "fraction"],
        composition_by_category,
    )

    write_thresholds(workflow_dir / args.thresholds_output, thresholds, args.positive_percentile)
    plot_counts_by_image(workflow_dir / args.counts_figure, output_rows)
    plot_composition_by_category(workflow_dir / args.category_figure, composition_by_category)

    phenotype_counts = Counter(row["phenotype"] for row in output_rows)
    print(f"Cells phenotyped: {len(output_rows)}")
    print(f"Phenotype classes assigned: {len(phenotype_counts)}")
    for phenotype, count in sorted(phenotype_counts.items()):
        print(f"{phenotype}: {count}")
    print(f"Phenotyped table: {workflow_dir / args.output_cells}")
    print(f"Composition by image: {workflow_dir / args.composition_image}")
    print(f"Composition by category: {workflow_dir / args.composition_category}")
    print(f"Threshold table: {workflow_dir / args.thresholds_output}")


if __name__ == "__main__":
    main()
