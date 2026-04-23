#!/usr/bin/env python3
"""Run the complete reusable IMC workflow in numbered script order."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


STEPS = [
    ("01", "Create marker panel", "01_create_panel_from_raw.py"),
    ("02", "Preprocess IMC images", "02_preprocess_imc_images.py"),
    ("03", "Segment images with Mesmer", "03_segment_mesmer.py"),
    ("04", "Measure marker intensities", "04_measure_intensities.py"),
    ("05", "Measure region properties", "05_measure_regionprops.py"),
    ("06", "Measure spatial neighbors", "06_measure_neighbors.py"),
    ("07", "Export single-cell data", "07_export_data.py"),
    ("08", "Process single-cell features", "08_process_single_cell_features.py"),
    ("09", "Run processed feature QC", "09_qc_processed_single_cell_features.py"),
    ("10", "Assign rule-based phenotypes", "10_rule_based_phenotyping.py"),
    ("11", "Summarize spatial phenotype interactions", "11_spatial_phenotype_interactions.py"),
    ("12", "Calculate spatial phenotype enrichment", "12_spatial_phenotype_enrichment.py"),
    ("13", "Summarize plasma-cell-like spatial niche", "13_plasma_cell_niche_analysis.py"),
    ("14", "Generate final interpretation report", "14_generate_final_interpretation_report.py"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the complete IMC reproduction workflow.")
    parser.add_argument("--workflow-dir", type=Path, default=Path.cwd())
    parser.add_argument("--start-at", choices=[step_id for step_id, _, _ in STEPS], default="01")
    parser.add_argument("--stop-after", choices=[step_id for step_id, _, _ in STEPS], default="14")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def selected_steps(start_at: str, stop_after: str) -> list[tuple[str, str, str]]:
    start_index = next(index for index, step in enumerate(STEPS) if step[0] == start_at)
    stop_index = next(index for index, step in enumerate(STEPS) if step[0] == stop_after)
    if start_index > stop_index:
        raise ValueError("--start-at must be earlier than or equal to --stop-after")
    return STEPS[start_index : stop_index + 1]


def require_raw_inputs(workflow_dir: Path) -> None:
    raw_dir = workflow_dir / "data" / "raw"
    if not raw_dir.exists():
        raise FileNotFoundError(f"Raw data directory not found: {raw_dir}")
    raw_files = [path for path in raw_dir.iterdir() if path.is_file()]
    if not raw_files:
        raise FileNotFoundError(f"No raw input files found in: {raw_dir}")


def run_step(workflow_dir: Path, script_name: str, dry_run: bool) -> None:
    script_path = workflow_dir / "scripts" / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Workflow script not found: {script_path}")
    cmd = [sys.executable, str(script_path), "--workflow-dir", str(workflow_dir)]
    print(" ".join(cmd))
    if dry_run:
        return
    subprocess.run(cmd, cwd=workflow_dir, check=True)


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()
    require_raw_inputs(workflow_dir)

    steps = selected_steps(args.start_at, args.stop_after)
    print(f"Workflow directory: {workflow_dir}")
    print(f"Selected steps: {steps[0][0]} through {steps[-1][0]}")
    if args.dry_run:
        print("Dry run enabled. Commands will be displayed without execution.")

    for step_id, description, script_name in steps:
        print(f"\n[{step_id}] {description}")
        run_step(workflow_dir, script_name, args.dry_run)

    print("\nWorkflow run complete.")
    print(f"Final report path: {workflow_dir / 'results' / '25_final_interpretation_summary.md'}")


if __name__ == "__main__":
    main()
