#!/usr/bin/env python3
"""Compute phenotype-neighbor enrichment relative to phenotype abundance."""

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute spatial phenotype enrichment.")
    parser.add_argument("--workflow-dir", type=Path, default=Path.cwd())
    parser.add_argument(
        "--cells",
        type=Path,
        default=Path("results/13_cells_with_phenotypes.csv"),
    )
    parser.add_argument(
        "--interactions-image",
        type=Path,
        default=Path("results/17_spatial_phenotype_interactions_by_image.csv"),
    )
    parser.add_argument(
        "--interactions-category",
        type=Path,
        default=Path("results/18_spatial_phenotype_interactions_by_category.csv"),
    )
    parser.add_argument(
        "--output-image",
        type=Path,
        default=Path("results/20_spatial_phenotype_enrichment_by_image.csv"),
    )
    parser.add_argument(
        "--output-category",
        type=Path,
        default=Path("results/21_spatial_phenotype_enrichment_by_category.csv"),
    )
    parser.add_argument(
        "--heatmap",
        type=Path,
        default=Path("figures/07_spatial_phenotype_enrichment_heatmap.png"),
    )
    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required {description} not found: {path}")


def image_stem(image_name: str) -> str:
    return image_name[:-5] if image_name.endswith(".tiff") else image_name


def read_abundance(cells_path: Path) -> tuple[dict[str, Counter[str]], dict[str, Counter[str]], set[str]]:
    by_image: dict[str, Counter[str]] = defaultdict(Counter)
    by_category: dict[str, Counter[str]] = defaultdict(Counter)
    phenotypes: set[str] = set()

    with cells_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"Image", "category", "phenotype"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing columns from phenotyped table: " + ", ".join(sorted(missing)))

        for row in reader:
            image = image_stem(row["Image"])
            category = row["category"]
            phenotype = row["phenotype"]
            by_image[image][phenotype] += 1
            by_category[category][phenotype] += 1
            phenotypes.add(phenotype)

    return by_image, by_category, phenotypes


def expected_pair_fraction(counts: Counter[str], phenotype_a: str, phenotype_b: str) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    p_a = counts[phenotype_a] / total
    p_b = counts[phenotype_b] / total
    if phenotype_a == phenotype_b:
        return p_a * p_b
    return 2 * p_a * p_b


def enrich_rows(
    interactions_path: Path,
    group_column: str,
    abundances: dict[str, Counter[str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with interactions_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {group_column, "phenotype_a", "phenotype_b", "edge_count", "fraction_of_group_edges"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing columns from interaction table: " + ", ".join(sorted(missing)))

        for row in reader:
            group = row[group_column]
            phenotype_a = row["phenotype_a"]
            phenotype_b = row["phenotype_b"]
            observed_fraction = float(row["fraction_of_group_edges"])
            expected_fraction = expected_pair_fraction(abundances[group], phenotype_a, phenotype_b)
            observed_expected_ratio = (
                observed_fraction / expected_fraction if expected_fraction > 0 else math.nan
            )
            log2_enrichment = (
                math.log2(observed_expected_ratio)
                if observed_expected_ratio > 0 and not math.isnan(observed_expected_ratio)
                else math.nan
            )

            rows.append(
                {
                    group_column: group,
                    "phenotype_a": phenotype_a,
                    "phenotype_b": phenotype_b,
                    "edge_count": int(row["edge_count"]),
                    "observed_edge_fraction": observed_fraction,
                    "expected_edge_fraction_from_abundance": expected_fraction,
                    "observed_expected_ratio": observed_expected_ratio,
                    "log2_enrichment": log2_enrichment,
                    "group_cell_count": sum(abundances[group].values()),
                    "phenotype_a_cell_count": abundances[group][phenotype_a],
                    "phenotype_b_cell_count": abundances[group][phenotype_b],
                }
            )
    return rows


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_category_heatmap(path: Path, rows: list[dict[str, object]], phenotypes: set[str]) -> None:
    # Aggregate category-level enrichment by phenotype pair using edge-count weighting.
    weighted: dict[tuple[str, str], list[tuple[float, int]]] = defaultdict(list)
    for row in rows:
        log2_value = float(row["log2_enrichment"])
        edge_count = int(row["edge_count"])
        if math.isnan(log2_value):
            continue
        key = (str(row["phenotype_a"]), str(row["phenotype_b"]))
        weighted[key].append((log2_value, edge_count))
        if key[0] != key[1]:
            weighted[(key[1], key[0])].append((log2_value, edge_count))

    labels = sorted(phenotypes)
    matrix: list[list[float]] = []
    for row_label in labels:
        matrix_row = []
        for col_label in labels:
            values = weighted.get((row_label, col_label), [])
            if not values:
                matrix_row.append(0.0)
                continue
            total_edges = sum(edge_count for _value, edge_count in values)
            if total_edges == 0:
                matrix_row.append(0.0)
            else:
                matrix_row.append(sum(value * edge_count for value, edge_count in values) / total_edges)
        matrix.append(matrix_row)

    max_abs = max((abs(value) for row in matrix for value in row), default=1.0)
    max_abs = max(max_abs, 1.0)

    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 8.5))
    plt.imshow(matrix, cmap="coolwarm", vmin=-max_abs, vmax=max_abs)
    plt.colorbar(label="log2 Observed / Expected")
    plt.xticks(range(len(labels)), labels, rotation=70, ha="right", fontsize=7)
    plt.yticks(range(len(labels)), labels, fontsize=7)
    plt.title("Spatial Phenotype Enrichment")
    plt.tight_layout()
    plt.savefig(path, dpi=170)
    plt.close()


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()
    cells_path = workflow_dir / args.cells
    interactions_image = workflow_dir / args.interactions_image
    interactions_category = workflow_dir / args.interactions_category

    require_input(cells_path, "phenotyped cell table")
    require_input(interactions_image, "image-level interaction table")
    require_input(interactions_category, "category-level interaction table")

    abundance_by_image, abundance_by_category, phenotypes = read_abundance(cells_path)
    image_rows = enrich_rows(interactions_image, "Image", abundance_by_image)
    category_rows = enrich_rows(interactions_category, "category", abundance_by_category)

    fieldnames_image = [
        "Image",
        "phenotype_a",
        "phenotype_b",
        "edge_count",
        "observed_edge_fraction",
        "expected_edge_fraction_from_abundance",
        "observed_expected_ratio",
        "log2_enrichment",
        "group_cell_count",
        "phenotype_a_cell_count",
        "phenotype_b_cell_count",
    ]
    fieldnames_category = fieldnames_image.copy()
    fieldnames_category[0] = "category"

    write_csv(workflow_dir / args.output_image, fieldnames_image, image_rows)
    write_csv(workflow_dir / args.output_category, fieldnames_category, category_rows)
    plot_category_heatmap(workflow_dir / args.heatmap, category_rows, phenotypes)

    print(f"Image-level enrichment rows: {len(image_rows)}")
    print(f"Category-level enrichment rows: {len(category_rows)}")
    print(f"Phenotypes represented: {len(phenotypes)}")
    print(f"Image-level enrichment table: {workflow_dir / args.output_image}")
    print(f"Category-level enrichment table: {workflow_dir / args.output_category}")
    print(f"Enrichment heatmap: {workflow_dir / args.heatmap}")


if __name__ == "__main__":
    main()
