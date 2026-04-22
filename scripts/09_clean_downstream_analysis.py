#!/usr/bin/env python3
"""Run clean downstream analysis from reproduced Steinbock outputs.

This script performs a compact representative-subset analysis:

1. QC of exported single-cell data.
2. Rule-based phenotyping from marker expression.
3. Phenotype composition by diagnostic category.
4. Phenotype-aware spatial interaction summaries from neighbor files.

The workflow is intentionally compact and produces only final downstream
outputs.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


PHENOTYPE_RULES = [
    ("plasma_myeloma_like", ["CD138", "CD38", "IRF4"], 2),
    ("t_cell_like", ["CD3", "CD4", "CD8"], 1),
    ("myeloid_macrophage_like", ["CD68", "CD11b", "CD11c", "MPO"], 2),
    ("stromal_bone_like", ["Vimentin", "CollagenTypeI", "RUNX2", "CathepsinK"], 2),
    ("proliferating_like", ["Ki67"], 1),
]

REGIONPROPS_COLUMNS = {"area", "centroid-0", "centroid-1", "axis_major_length", "axis_minor_length", "eccentricity"}
NON_MARKER_COLUMNS = {"Image", "Object"} | REGIONPROPS_COLUMNS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run clean downstream IMC analysis.")
    parser.add_argument("--workflow-dir", type=Path, default=Path.cwd())
    parser.add_argument("--cells-csv", type=Path, default=Path("results/08_cells.csv"))
    parser.add_argument("--neighbors-dir", type=Path, default=Path("data/neighbors"))
    parser.add_argument("--threshold-quantile", type=float, default=0.75)
    parser.add_argument("--cofactor", type=float, default=1.0)
    return parser.parse_args()


def numeric(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return float("nan")


def quantile(values: list[float], q: float) -> float:
    clean = sorted(value for value in values if not math.isnan(value))
    if not clean:
        return float("nan")
    index = (len(clean) - 1) * q
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return clean[int(index)]
    return clean[lower] * (upper - index) + clean[upper] * (index - lower)


def parse_image_metadata(image_name: str) -> dict[str, str]:
    stem = image_name.replace(".tiff", "")
    parts = stem.split("_")
    category = "Unknown"
    imc_id = "Unknown"
    for part in parts:
        if part.startswith("IMC"):
            imc_id = part
        if part in {"UB", "MGUS", "SMM", "B"}:
            category = part
    category_index = parts.index(category) if category in parts else len(parts) - 1
    case_id = "_".join(parts[: category_index + 1])
    roi_id = stem.replace(case_id + "_", "") if stem.startswith(case_id + "_") else ""
    return {"case_id": case_id, "imc_id": imc_id, "category": category, "roi_id": roi_id}


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_cells(path: Path, cofactor: float) -> tuple[list[dict[str, str]], list[str], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        original_fieldnames = list(reader.fieldnames)
        marker_columns = [col for col in original_fieldnames if col not in NON_MARKER_COLUMNS]
        rows = []
        for row in reader:
            row.update(parse_image_metadata(row["Image"]))
            for marker in marker_columns:
                value = numeric(row[marker])
                row[f"asinh_{marker}"] = math.asinh(value / cofactor) if not math.isnan(value) else float("nan")
            rows.append(row)
    return rows, original_fieldnames, marker_columns


def marker_thresholds(rows: list[dict[str, str]], q: float) -> dict[str, float]:
    selected_markers = sorted({marker for _, markers, _ in PHENOTYPE_RULES for marker in markers})
    thresholds = {}
    for marker in selected_markers:
        col = f"asinh_{marker}"
        if col in rows[0]:
            thresholds[marker] = quantile([numeric(row[col]) for row in rows], q)
    return thresholds


def assign_phenotype(row: dict[str, str], thresholds: dict[str, float]) -> tuple[str, str]:
    passing = []
    for phenotype, markers, minimum in PHENOTYPE_RULES:
        available = [marker for marker in markers if marker in thresholds]
        score = sum(numeric(row[f"asinh_{marker}"]) >= thresholds[marker] for marker in available)
        if score >= minimum:
            passing.append((phenotype, score))
    if not passing:
        return "unassigned", "no_rule_passed"
    passing.sort(key=lambda item: item[1], reverse=True)
    if len(passing) > 1 and passing[0][1] == passing[1][1]:
        return passing[0][0], "multiple_rules_tied"
    return passing[0][0], "single_best_rule"


def add_phenotypes(rows: list[dict[str, str]], thresholds: dict[str, float]) -> None:
    for row in rows:
        phenotype, note = assign_phenotype(row, thresholds)
        row["phenotype"] = phenotype
        row["phenotype_note"] = note


def qc_summary(rows: list[dict[str, str]], marker_columns: list[str]) -> list[dict[str, object]]:
    category_counts = Counter(row["category"] for row in rows)
    image_counts = Counter(row["Image"] for row in rows)
    summary = [
        {"metric": "cells", "group": "all", "value": len(rows)},
        {"metric": "images", "group": "all", "value": len(image_counts)},
        {"metric": "categories", "group": "all", "value": len(category_counts)},
        {"metric": "marker_columns", "group": "all", "value": len(marker_columns)},
    ]
    for category, count in sorted(category_counts.items()):
        summary.append({"metric": "cells_by_category", "group": category, "value": count})
    return summary


def composition_by_category(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["category"]].append(row)
    phenotypes = sorted({row["phenotype"] for row in rows})
    out = []
    for category, group_rows in sorted(grouped.items()):
        counts = Counter(row["phenotype"] for row in group_rows)
        total = len(group_rows)
        for phenotype in phenotypes:
            count = counts[phenotype]
            out.append(
                {
                    "category": category,
                    "phenotype": phenotype,
                    "cell_count": count,
                    "fraction": count / total if total else 0,
                }
            )
    return out


def normalize_object_id(value: str) -> str:
    try:
        return str(int(float(value)))
    except ValueError:
        return value


def neighbor_columns(fieldnames: list[str]) -> tuple[str, str]:
    lower = {name.casefold(): name for name in fieldnames}
    if "object" in lower and "neighbor" in lower:
        return lower["object"], lower["neighbor"]
    if len(fieldnames) >= 2:
        return fieldnames[0], fieldnames[1]
    raise ValueError("Neighbor table has fewer than two columns.")


def spatial_interactions(rows: list[dict[str, str]], neighbors_dir: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    phenotype_by_cell = {
        (row["Image"], normalize_object_id(row["Object"])): row["phenotype"]
        for row in rows
    }
    phenotype_counts = Counter(row["phenotype"] for row in rows)
    pair_counts: Counter[tuple[str, str]] = Counter()

    for path in sorted(neighbors_dir.glob("*.csv")):
        image = f"{path.stem}.tiff"
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                continue
            left_col, right_col = neighbor_columns(reader.fieldnames)
            for edge in reader:
                left = normalize_object_id(edge[left_col])
                right = normalize_object_id(edge[right_col])
                left_pheno = phenotype_by_cell.get((image, left))
                right_pheno = phenotype_by_cell.get((image, right))
                if left_pheno is None or right_pheno is None:
                    continue
                pair_counts[tuple(sorted([left_pheno, right_pheno]))] += 1

    total_edges = sum(pair_counts.values())
    total_cells = sum(phenotype_counts.values())
    interaction_rows = []
    for (phenotype_a, phenotype_b), count in sorted(pair_counts.items()):
        interaction_rows.append(
            {
                "phenotype_a": phenotype_a,
                "phenotype_b": phenotype_b,
                "edge_count": count,
                "edge_fraction": count / total_edges if total_edges else 0,
            }
        )

    enrichment_rows = []
    phenotypes = sorted(phenotype_counts)
    for i, phenotype_a in enumerate(phenotypes):
        for phenotype_b in phenotypes[i:]:
            observed = pair_counts.get(tuple(sorted([phenotype_a, phenotype_b])), 0)
            p_a = phenotype_counts[phenotype_a] / total_cells if total_cells else 0
            p_b = phenotype_counts[phenotype_b] / total_cells if total_cells else 0
            expected_fraction = p_a * p_b if phenotype_a == phenotype_b else 2 * p_a * p_b
            expected_edges = expected_fraction * total_edges
            enrichment_rows.append(
                {
                    "phenotype_a": phenotype_a,
                    "phenotype_b": phenotype_b,
                    "observed_edges": observed,
                    "expected_edges": expected_edges,
                    "observed_expected_ratio": observed / expected_edges if expected_edges else 0,
                }
            )
    return interaction_rows, enrichment_rows


def plot_cell_counts(path: Path, summary_rows: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    rows = [row for row in summary_rows if row["metric"] == "cells_by_category"]
    labels = [str(row["group"]) for row in rows]
    values = [int(row["value"]) for row in rows]
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(labels, values, color="#3d5a80")
    ax.set_xlabel("Category")
    ax.set_ylabel("Cell count")
    ax.set_title("Cell Counts by Category")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_composition(path: Path, rows: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    categories = sorted({str(row["category"]) for row in rows})
    phenotypes = sorted({str(row["phenotype"]) for row in rows})
    fractions = {(row["category"], row["phenotype"]): float(row["fraction"]) for row in rows}
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(max(7, len(categories) * 1.4), 5))
    bottom = [0.0] * len(categories)
    for phenotype in phenotypes:
        values = [fractions.get((category, phenotype), 0.0) for category in categories]
        ax.bar(categories, values, bottom=bottom, label=phenotype)
        bottom = [left + value for left, value in zip(bottom, values)]
    ax.set_xlabel("Category")
    ax.set_ylabel("Fraction of cells")
    ax.set_title("Phenotype Composition by Category")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_enrichment(path: Path, rows: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    phenotypes = sorted({str(row["phenotype_a"]) for row in rows} | {str(row["phenotype_b"]) for row in rows})
    index = {phenotype: idx for idx, phenotype in enumerate(phenotypes)}
    matrix = [[0.0 for _ in phenotypes] for _ in phenotypes]
    for row in rows:
        a = str(row["phenotype_a"])
        b = str(row["phenotype_b"])
        value = float(row["observed_expected_ratio"])
        matrix[index[a]][index[b]] = value
        matrix[index[b]][index[a]] = value
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(max(6, len(phenotypes) * 1.0), max(5, len(phenotypes) * 0.8)))
    im = ax.imshow(matrix, cmap="coolwarm", vmin=0)
    ax.set_xticks(range(len(phenotypes)))
    ax.set_xticklabels(phenotypes, rotation=45, ha="right")
    ax.set_yticks(range(len(phenotypes)))
    ax.set_yticklabels(phenotypes)
    ax.set_title("Phenotype Spatial Enrichment")
    fig.colorbar(im, ax=ax, label="Observed / expected")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()
    cells_csv = workflow_dir / args.cells_csv
    neighbors_dir = workflow_dir / args.neighbors_dir

    if not cells_csv.exists():
        raise FileNotFoundError(f"Required cells table not found: {cells_csv}")
    if not neighbors_dir.exists():
        raise FileNotFoundError(f"Required neighbors directory not found: {neighbors_dir}")

    rows, original_fieldnames, marker_columns = load_cells(cells_csv, args.cofactor)
    thresholds = marker_thresholds(rows, args.threshold_quantile)
    add_phenotypes(rows, thresholds)

    qc_rows = qc_summary(rows, marker_columns)
    composition_rows = composition_by_category(rows)
    interaction_rows, enrichment_rows = spatial_interactions(rows, neighbors_dir)

    output_fieldnames = original_fieldnames + ["case_id", "imc_id", "category", "roi_id"] + [
        f"asinh_{marker}" for marker in marker_columns
    ] + ["phenotype", "phenotype_note"]

    write_csv(workflow_dir / "results/12_cell_qc_summary.csv", qc_rows, ["metric", "group", "value"])
    write_csv(workflow_dir / "results/13_cells_with_phenotypes.csv", rows, output_fieldnames)
    write_csv(workflow_dir / "results/14_phenotype_composition_by_category.csv", composition_rows, ["category", "phenotype", "cell_count", "fraction"])
    write_csv(workflow_dir / "results/15_phenotype_spatial_interactions.csv", interaction_rows, ["phenotype_a", "phenotype_b", "edge_count", "edge_fraction"])
    write_csv(workflow_dir / "results/16_phenotype_spatial_enrichment.csv", enrichment_rows, ["phenotype_a", "phenotype_b", "observed_edges", "expected_edges", "observed_expected_ratio"])

    plot_cell_counts(workflow_dir / "figures/01_cell_counts_by_category.png", qc_rows)
    plot_composition(workflow_dir / "figures/02_phenotype_composition_by_category.png", composition_rows)
    plot_enrichment(workflow_dir / "figures/03_spatial_enrichment_heatmap.png", enrichment_rows)

    print("Clean downstream analysis complete.")
    print("results/12_cell_qc_summary.csv")
    print("results/13_cells_with_phenotypes.csv")
    print("results/14_phenotype_composition_by_category.csv")
    print("results/15_phenotype_spatial_interactions.csv")
    print("results/16_phenotype_spatial_enrichment.csv")
    print("figures/01_cell_counts_by_category.png")
    print("figures/02_phenotype_composition_by_category.png")
    print("figures/03_spatial_enrichment_heatmap.png")


if __name__ == "__main__":
    main()
