#!/usr/bin/env python3
"""Summarize spatial niches around Plasma_cell_like cells."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path

Path(".matplotlib").mkdir(exist_ok=True)
Path(".cache").mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib").resolve()))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(".cache").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


FOCUS = "Plasma_cell_like"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze Plasma_cell_like spatial neighbor niches.")
    parser.add_argument("--workflow-dir", type=Path, default=Path.cwd())
    parser.add_argument("--composition-category", type=Path, default=Path("results/15_phenotype_composition_by_category.csv"))
    parser.add_argument("--enrichment-image", type=Path, default=Path("results/20_spatial_phenotype_enrichment_by_image.csv"))
    parser.add_argument("--enrichment-category", type=Path, default=Path("results/21_spatial_phenotype_enrichment_by_category.csv"))
    parser.add_argument("--min-summary-edges", type=int, default=25)
    parser.add_argument("--output-image", type=Path, default=Path("results/22_plasma_cell_niche_by_image.csv"))
    parser.add_argument("--output-category", type=Path, default=Path("results/23_plasma_cell_niche_by_category.csv"))
    parser.add_argument("--output-summary", type=Path, default=Path("results/24_plasma_cell_niche_summary.csv"))
    parser.add_argument("--figure", type=Path, default=Path("figures/08_plasma_cell_neighbor_enrichment.png"))
    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required {description} not found: {path}")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plasma_partner(row: dict[str, str]) -> str | None:
    a = row["phenotype_a"]
    b = row["phenotype_b"]
    if a == FOCUS and b == FOCUS:
        return FOCUS
    if a == FOCUS:
        return b
    if b == FOCUS:
        return a
    return None


def filter_plasma_rows(rows: list[dict[str, str]], group_column: str) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for row in rows:
        partner = plasma_partner(row)
        if partner is None:
            continue
        output.append(
            {
                group_column: row[group_column],
                "neighbor_phenotype": partner,
                "edge_count": int(row["edge_count"]),
                "observed_edge_fraction": float(row["observed_edge_fraction"]),
                "expected_edge_fraction_from_abundance": float(row["expected_edge_fraction_from_abundance"]),
                "observed_expected_ratio": float(row["observed_expected_ratio"]),
                "log2_enrichment": float(row["log2_enrichment"]),
                "group_cell_count": int(row["group_cell_count"]),
                "plasma_cell_count": int(row["phenotype_a_cell_count"] if row["phenotype_a"] == FOCUS else row["phenotype_b_cell_count"]),
                "neighbor_cell_count": int(row["phenotype_b_cell_count"] if row["phenotype_a"] == FOCUS else row["phenotype_a_cell_count"]),
            }
        )
    return output


def summarize_category_rows(rows: list[dict[str, object]], min_summary_edges: int) -> list[dict[str, object]]:
    summary: list[dict[str, object]] = []
    categories = sorted({str(row["category"]) for row in rows})
    for category in categories:
        category_rows = [row for row in rows if row["category"] == category]
        reliable_rows = [row for row in category_rows if int(row["edge_count"]) >= min_summary_edges]
        ranking_rows = reliable_rows if reliable_rows else category_rows
        enriched = [row for row in category_rows if float(row["log2_enrichment"]) > 0]
        depleted = [row for row in category_rows if float(row["log2_enrichment"]) < 0]
        top_enriched = max(ranking_rows, key=lambda row: float(row["log2_enrichment"])) if ranking_rows else None
        top_depleted = min(ranking_rows, key=lambda row: float(row["log2_enrichment"])) if ranking_rows else None
        summary.append(
            {
                "category": category,
                "plasma_neighbor_pairs": len(category_rows),
                "minimum_edges_for_top_neighbor_summary": min_summary_edges,
                "plasma_neighbor_pairs_passing_minimum_edges": len(reliable_rows),
                "enriched_pairs_log2_gt_0": len(enriched),
                "depleted_pairs_log2_lt_0": len(depleted),
                "top_enriched_neighbor": top_enriched["neighbor_phenotype"] if top_enriched else "",
                "top_enriched_log2": top_enriched["log2_enrichment"] if top_enriched else math.nan,
                "top_enriched_edges": top_enriched["edge_count"] if top_enriched else 0,
                "top_depleted_neighbor": top_depleted["neighbor_phenotype"] if top_depleted else "",
                "top_depleted_log2": top_depleted["log2_enrichment"] if top_depleted else math.nan,
                "top_depleted_edges": top_depleted["edge_count"] if top_depleted else 0,
            }
        )
    return summary


def add_plasma_abundance(summary: list[dict[str, object]], composition_rows: list[dict[str, str]]) -> None:
    abundance = {}
    for row in composition_rows:
        if row["phenotype"] == FOCUS:
            abundance[row["category"]] = {
                "plasma_cell_count": int(row["cell_count"]),
                "plasma_cell_fraction": float(row["fraction"]),
            }
    for row in summary:
        row.update(abundance.get(str(row["category"]), {"plasma_cell_count": 0, "plasma_cell_fraction": 0.0}))


def plot_plasma_enrichment(path: Path, rows: list[dict[str, object]]) -> None:
    categories = sorted({str(row["category"]) for row in rows})
    partners = sorted({str(row["neighbor_phenotype"]) for row in rows})
    values = {
        (str(row["category"]), str(row["neighbor_phenotype"])): float(row["log2_enrichment"])
        for row in rows
    }
    matrix = [[values.get((category, partner), 0.0) for partner in partners] for category in categories]
    max_abs = max((abs(value) for row in matrix for value in row), default=1.0)
    max_abs = max(max_abs, 1.0)

    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(12, 4.8))
    plt.imshow(matrix, cmap="coolwarm", vmin=-max_abs, vmax=max_abs, aspect="auto")
    plt.colorbar(label="log2 Observed / Expected")
    plt.xticks(range(len(partners)), partners, rotation=70, ha="right", fontsize=8)
    plt.yticks(range(len(categories)), categories, fontsize=9)
    plt.title("Plasma-cell-like Neighbor Enrichment by Category")
    plt.tight_layout()
    plt.savefig(path, dpi=170)
    plt.close()


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()
    composition_category = workflow_dir / args.composition_category
    enrichment_image = workflow_dir / args.enrichment_image
    enrichment_category = workflow_dir / args.enrichment_category

    require_input(composition_category, "phenotype composition by category")
    require_input(enrichment_image, "image-level enrichment table")
    require_input(enrichment_category, "category-level enrichment table")

    image_rows = filter_plasma_rows(read_rows(enrichment_image), "Image")
    category_rows = filter_plasma_rows(read_rows(enrichment_category), "category")
    summary_rows = summarize_category_rows(category_rows, args.min_summary_edges)
    add_plasma_abundance(summary_rows, read_rows(composition_category))

    image_fields = [
        "Image",
        "neighbor_phenotype",
        "edge_count",
        "observed_edge_fraction",
        "expected_edge_fraction_from_abundance",
        "observed_expected_ratio",
        "log2_enrichment",
        "group_cell_count",
        "plasma_cell_count",
        "neighbor_cell_count",
    ]
    category_fields = image_fields.copy()
    category_fields[0] = "category"
    summary_fields = [
        "category",
        "plasma_cell_count",
        "plasma_cell_fraction",
        "plasma_neighbor_pairs",
        "minimum_edges_for_top_neighbor_summary",
        "plasma_neighbor_pairs_passing_minimum_edges",
        "enriched_pairs_log2_gt_0",
        "depleted_pairs_log2_lt_0",
        "top_enriched_neighbor",
        "top_enriched_log2",
        "top_enriched_edges",
        "top_depleted_neighbor",
        "top_depleted_log2",
        "top_depleted_edges",
    ]

    write_csv(workflow_dir / args.output_image, image_fields, image_rows)
    write_csv(workflow_dir / args.output_category, category_fields, category_rows)
    write_csv(workflow_dir / args.output_summary, summary_fields, summary_rows)
    plot_plasma_enrichment(workflow_dir / args.figure, category_rows)

    print(f"Plasma niche image rows: {len(image_rows)}")
    print(f"Plasma niche category rows: {len(category_rows)}")
    print(f"Plasma niche summary rows: {len(summary_rows)}")
    print(f"Plasma niche figure: {workflow_dir / args.figure}")


if __name__ == "__main__":
    main()
