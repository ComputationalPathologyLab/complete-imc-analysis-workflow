#!/usr/bin/env python3
"""Summarize spatial neighbor interactions between rule-based phenotypes."""

from __future__ import annotations

import argparse
import csv
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
    parser = argparse.ArgumentParser(
        description="Summarize phenotype-phenotype neighbor interactions."
    )
    parser.add_argument("--workflow-dir", type=Path, default=Path.cwd())
    parser.add_argument(
        "--cells",
        type=Path,
        default=Path("results/13_cells_with_phenotypes.csv"),
        help="Phenotyped single-cell table, relative to workflow-dir.",
    )
    parser.add_argument(
        "--neighbors-dir",
        type=Path,
        default=Path("data/neighbors"),
        help="Directory containing Steinbock neighbor CSV files, relative to workflow-dir.",
    )
    parser.add_argument(
        "--output-image",
        type=Path,
        default=Path("results/17_spatial_phenotype_interactions_by_image.csv"),
    )
    parser.add_argument(
        "--output-category",
        type=Path,
        default=Path("results/18_spatial_phenotype_interactions_by_category.csv"),
    )
    parser.add_argument(
        "--output-summary",
        type=Path,
        default=Path("results/19_spatial_phenotype_interaction_summary.csv"),
    )
    parser.add_argument(
        "--heatmap",
        type=Path,
        default=Path("figures/06_spatial_phenotype_interaction_heatmap.png"),
    )
    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required {description} not found: {path}")


def image_stem(image_name: str) -> str:
    return image_name[:-5] if image_name.endswith(".tiff") else image_name


def canonical_pair(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))


def read_cell_lookup(cells_path: Path) -> tuple[dict[tuple[str, str], dict[str, str]], set[str]]:
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    phenotypes: set[str] = set()
    with cells_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"Image", "Object", "category", "phenotype"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing required columns from phenotyped table: " + ", ".join(sorted(missing)))

        for row in reader:
            stem = image_stem(row["Image"])
            lookup[(stem, row["Object"])] = {
                "category": row["category"],
                "phenotype": row["phenotype"],
            }
            phenotypes.add(row["phenotype"])
    return lookup, phenotypes


def summarize_interactions(
    lookup: dict[tuple[str, str], dict[str, str]],
    neighbors_dir: Path,
) -> tuple[Counter[tuple[str, str, str]], Counter[tuple[str, str, str]], dict[str, int]]:
    by_image: Counter[tuple[str, str, str]] = Counter()
    by_category: Counter[tuple[str, str, str]] = Counter()
    diagnostics = {
        "neighbor_files": 0,
        "raw_directed_edges": 0,
        "unique_undirected_edges": 0,
        "edges_skipped_missing_cell": 0,
    }

    for neighbor_file in sorted(neighbors_dir.glob("*.csv")):
        diagnostics["neighbor_files"] += 1
        image = neighbor_file.stem
        seen_edges: set[tuple[str, str]] = set()
        with neighbor_file.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                diagnostics["raw_directed_edges"] += 1
                obj_a = row["Object"]
                obj_b = row["Neighbor"]
                if obj_a == obj_b:
                    continue
                edge = tuple(sorted((obj_a, obj_b)))
                if edge in seen_edges:
                    continue
                seen_edges.add(edge)

                cell_a = lookup.get((image, obj_a))
                cell_b = lookup.get((image, obj_b))
                if cell_a is None or cell_b is None:
                    diagnostics["edges_skipped_missing_cell"] += 1
                    continue

                phenotype_a, phenotype_b = canonical_pair(cell_a["phenotype"], cell_b["phenotype"])
                category = cell_a["category"]
                by_image[(image, phenotype_a, phenotype_b)] += 1
                by_category[(category, phenotype_a, phenotype_b)] += 1

        diagnostics["unique_undirected_edges"] += len(seen_edges)

    return by_image, by_category, diagnostics


def rows_from_counter(
    counter: Counter[tuple[str, str, str]],
    group_name: str,
) -> list[dict[str, object]]:
    totals: Counter[str] = Counter()
    for key, count in counter.items():
        totals[key[0]] += count

    rows = []
    for (group, phenotype_a, phenotype_b), count in sorted(counter.items()):
        total = totals[group]
        rows.append(
            {
                group_name: group,
                "phenotype_a": phenotype_a,
                "phenotype_b": phenotype_b,
                "edge_count": count,
                "fraction_of_group_edges": count / total if total else 0,
            }
        )
    return rows


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, diagnostics: dict[str, int], by_image_rows: list[dict[str, object]]) -> None:
    rows = [{"metric": key, "value": value} for key, value in diagnostics.items()]
    rows.append({"metric": "image_level_interaction_rows", "value": len(by_image_rows)})
    write_csv(path, ["metric", "value"], rows)


def plot_heatmap(path: Path, by_category: Counter[tuple[str, str, str]], phenotypes: set[str]) -> None:
    # Aggregate across categories for a compact overview figure.
    matrix_counts: Counter[tuple[str, str]] = Counter()
    for (_category, phenotype_a, phenotype_b), count in by_category.items():
        matrix_counts[(phenotype_a, phenotype_b)] += count
        if phenotype_a != phenotype_b:
            matrix_counts[(phenotype_b, phenotype_a)] += count

    labels = sorted(phenotypes)
    if not labels:
        return

    matrix = [[matrix_counts[(row_label, col_label)] for col_label in labels] for row_label in labels]

    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 8.5))
    plt.imshow(matrix, cmap="viridis")
    plt.colorbar(label="Neighbor Edge Count")
    plt.xticks(range(len(labels)), labels, rotation=70, ha="right", fontsize=7)
    plt.yticks(range(len(labels)), labels, fontsize=7)
    plt.title("Spatial Phenotype Neighbor Interactions")
    plt.tight_layout()
    plt.savefig(path, dpi=170)
    plt.close()


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()
    cells_path = workflow_dir / args.cells
    neighbors_dir = workflow_dir / args.neighbors_dir
    require_input(cells_path, "phenotyped cell table")
    require_input(neighbors_dir, "neighbor directory")

    lookup, phenotypes = read_cell_lookup(cells_path)
    by_image, by_category, diagnostics = summarize_interactions(lookup, neighbors_dir)

    by_image_rows = rows_from_counter(by_image, "Image")
    by_category_rows = rows_from_counter(by_category, "category")

    write_csv(
        workflow_dir / args.output_image,
        ["Image", "phenotype_a", "phenotype_b", "edge_count", "fraction_of_group_edges"],
        by_image_rows,
    )
    write_csv(
        workflow_dir / args.output_category,
        ["category", "phenotype_a", "phenotype_b", "edge_count", "fraction_of_group_edges"],
        by_category_rows,
    )
    write_summary(workflow_dir / args.output_summary, diagnostics, by_image_rows)
    plot_heatmap(workflow_dir / args.heatmap, by_category, phenotypes)

    print(f"Neighbor files processed: {diagnostics['neighbor_files']}")
    print(f"Raw directed edges read: {diagnostics['raw_directed_edges']}")
    print(f"Unique undirected edges considered: {diagnostics['unique_undirected_edges']}")
    print(f"Edges skipped because phenotype lookup was missing: {diagnostics['edges_skipped_missing_cell']}")
    print(f"Image-level interaction rows: {len(by_image_rows)}")
    print(f"Category-level interaction rows: {len(by_category_rows)}")
    print(f"Interaction heatmap: {workflow_dir / args.heatmap}")


if __name__ == "__main__":
    main()
