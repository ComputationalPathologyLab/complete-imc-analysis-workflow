# Step-by-Step IMC Reproduction

![Workflow](https://img.shields.io/badge/workflow-Steinbock%20IMC-2f6f73)
![Status](https://img.shields.io/badge/status-active%20reproduction-3b7a57)
![Data](https://img.shields.io/badge/data-representative%20template%20subset-5b6f95)
![Segmentation](https://img.shields.io/badge/segmentation-DeepCell%20Mesmer-7a5c99)
![Features](https://img.shields.io/badge/features-arcsinh%20cofactor%201-8a6f3d)
![Phenotyping](https://img.shields.io/badge/phenotyping-rule--based%20approximation-9a5c5c)

## Overview

This directory contains a notebook-first reproduction workspace for a template
imaging mass cytometry (IMC) study of the multiple myeloma bone marrow
microenvironment. The workflow follows the structure of the Steinbock hands-on
pipeline, then adds documented single-cell feature processing, quality control,
and transparent first-pass phenotyping.

The immediate aim is to reproduce the study methodology on a representative
template subset. The longer-term aim is to keep the workflow sufficiently clear,
modular, and documented so that the same structure can be adapted to independent
IMC datasets.

## Study-Inspired Workflow

The workflow is organized as a linear sequence from raw IMC files to
analysis-ready single-cell tables.

```text
raw MCD / ROI text files
-> panel generation
-> Steinbock preprocessing
-> Mesmer segmentation
-> intensity, morphology, and neighbor measurement
-> CSV / AnnData / GraphML export
-> single-cell feature processing
-> processed feature QC
-> rule-based first-pass phenotyping
-> spatial phenotype interactions
-> abundance-normalized spatial enrichment
-> plasma-cell-like niche summary
-> final interpretation report
```

## Reproduction Principles

- All generated files remain inside `step_by_step_reproduction/`.
- Raw input data are read from `data/raw/`.
- Processed outputs from the original study are not used as inputs.
- Original study outputs may be used only for comparison and validation.
- Workflow execution is notebook-first and stepwise.
- Reusable scripts are stored in `scripts/` and called from notebook cells.
- Research outputs are numbered chronologically in `results/` and `figures/`.
- Interpretation boundaries are documented next to the relevant workflow step.

## Reuse on New Data

To apply the same workflow to another compatible IMC dataset, place the raw
`.mcd` files and ROI `.txt` files in `data/raw/`, then run the complete workflow
runner from the project directory:

```bash
python3 scripts/15_run_complete_workflow.py
```

For a non-executing command preview:

```bash
python3 scripts/15_run_complete_workflow.py --dry-run
```

For downstream regeneration after `data/cells.csv` and neighbor tables already
exist:

```bash
python3 scripts/15_run_complete_workflow.py --start-at 08
```

The runner executes the numbered scripts in order and regenerates the final
Markdown interpretation report at:

```text
results/25_final_interpretation_summary.md
```

## Notebooks

| Notebook | Role |
|---|---|
| `notebooks/01_Complete IMC Data Analysis Workflow.ipynb` | Raw processing, segmentation, measurement, neighbor construction, and export |
| `notebooks/02_single_cell_feature_processing.ipynb` | Area filtering, 99th percentile censoring, and arcsinh transformation |
| `notebooks/03_processed_single_cell_qc.ipynb` | Technical QC of processed single-cell features |
| `notebooks/04_rule_based_phenotyping.ipynb` | Transparent first-pass rule-based phenotype assignment |
| `notebooks/05_spatial_phenotype_interactions.ipynb` | Raw phenotype-neighbor interaction summaries |
| `notebooks/06_spatial_phenotype_enrichment.ipynb` | Observed-versus-expected spatial enrichment analysis |
| `notebooks/07_plasma_cell_spatial_niche_summary.ipynb` | Plasma-cell-like focused spatial niche summary |
| `notebooks/08_final_analysis_report.ipynb` | Data-driven final interpretation report generation |

Each notebook contains professional Markdown documentation before executable
code cells. The notebooks are intended to serve as the readable workflow record;
the scripts provide reusable implementations of the documented steps.

## Representative Dataset

The current representative subset contains four IMC cases, each with one raw
`.mcd` file and two ROI text files.

| Case | Category | Raw `.mcd` | ROI text files |
|---|---|---:|---:|
| `TS-373_IMC01_UB` | UB | 1 | 2 |
| `TS-373_IMC02_MGUS` | MGUS | 1 | 2 |
| `TS-373_IMC05_MGUS` | MGUS | 1 | 2 |
| `TS-373_IMC09_B` | B | 1 | 2 |

Current input total:

- 4 raw `.mcd` files
- 8 ROI `.txt` files
- 8 retained 1000 x 1000 MCD-derived ROI images after preprocessing cleanup

## Directory Structure

| Path | Purpose |
|---|---|
| `data/raw/` | Representative raw `.mcd` files and ROI text files |
| `data/panel.csv` | Steinbock-compatible marker/channel panel |
| `data/img/` | Preprocessed multi-channel TIFF images |
| `data/masks/` | Mesmer segmentation masks |
| `data/intensities/` | Per-cell marker intensity tables |
| `data/regionprops/` | Per-cell morphology and coordinate tables |
| `data/neighbors/` | Spatial cell-neighbor edge lists |
| `data/cells.csv` | Steinbock-exported single-cell table |
| `data/cells.h5ad` | Steinbock-exported AnnData object |
| `data/graphs/` | Spatial GraphML exports |
| `notebooks/` | Documented notebook workflow |
| `scripts/` | Numbered reusable helper scripts |
| `logs/` | Command logs and troubleshooting records |
| `results/` | Numbered tables and analysis outputs |
| `figures/` | Numbered QC and summary figures |

## Segmentation Panel

The panel contains 42 IMC channels. Mesmer segmentation markers are assigned in
`data/panel.csv` using the `deepcell` column.

| Marker role | Channel | Marker | `deepcell` value |
|---|---|---|---:|
| Nuclear | `Yb171` | `HistoneH3` | 1 |
| Nuclear | `Ir191` | `191Ir` | 1 |
| Nuclear | `Ir193` | `193Ir` | 1 |
| Membrane | `Sm152` | `CD45` | 2 |
| Membrane | `Er170` | `CD3` | 2 |
| Membrane | `Yb173` | `CD98` | 2 |
| Membrane | `Yb176` | `CD138` | 2 |

## Workflow Modules

### 1. Panel Generation

Script:

```text
scripts/01_create_panel_from_raw.py
```

Outputs:

```text
data/panel.csv
results/01_panel.csv
```

Purpose:

- Parse ROI text-file headers.
- Extract channel and marker names.
- Create a Steinbock-compatible panel.
- Assign nuclear and membrane markers for Mesmer segmentation.

### 2. IMC Preprocessing

Script:

```text
scripts/02_preprocess_imc_images.py
```

Command implemented:

```bash
steinbock preprocess imc images --hpf 50
```

Outputs:

```text
data/img/*.tiff
data/images.csv
results/02_images.csv
results/03_preprocessed_tiff_inventory.csv
```

Purpose:

- Convert raw IMC acquisitions to multi-channel TIFF images.
- Apply hot-pixel filtering with `--hpf 50`.
- Retain the study-style image set by keeping MCD-derived 1000 x 1000 ROIs.
- Remove tiny test acquisitions and duplicate images generated directly from ROI
  text files.

### 3. Mesmer Segmentation

Script:

```text
scripts/03_segment_mesmer.py
```

Command implemented:

```bash
steinbock segment deepcell --app mesmer --minmax
```

Outputs:

```text
data/masks/*.tiff
results/04_mesmer_mask_inventory.csv
```

Purpose:

- Generate whole-cell segmentation masks.
- Use combined nuclear and membrane marker inputs from `data/panel.csv`.
- Save one segmentation mask per retained ROI image.

### 4. Single-Cell Measurement

Scripts:

```text
scripts/04_measure_intensities.py
scripts/05_measure_regionprops.py
scripts/06_measure_neighbors.py
```

Commands implemented:

```bash
steinbock measure intensities
steinbock measure regionprops
steinbock measure neighbors --type expansion --dmax 4
```

Outputs:

```text
data/intensities/*.csv
data/regionprops/*.csv
data/neighbors/*.csv
results/05_intensity_table_inventory.csv
results/06_regionprops_table_inventory.csv
results/07_neighbors_table_inventory.csv
```

Purpose:

- Measure per-cell marker intensities from pixels inside each segmentation mask.
- Measure cell morphology and coordinates, including area and centroid.
- Construct spatial neighbor edge lists using boundary expansion up to 4 pixels.

### 5. Data Export

Script:

```text
scripts/07_export_data.py
```

Commands implemented:

```bash
steinbock export csv intensities regionprops -o cells.csv
steinbock export anndata --intensities intensities --data regionprops --neighbors neighbors -o cells.h5ad
steinbock export graphs --format graphml --data intensities --data regionprops
```

Outputs:

```text
data/cells.csv
data/cells.h5ad
data/graphs/*.graphml
results/08_cells.csv
results/09_cells.h5ad
results/10_graphml/
results/10_graphml_inventory.csv
```

Purpose:

- Combine marker intensities and region properties into one cell-level table.
- Export an AnnData object for downstream single-cell workflows.
- Export attributed spatial graphs for graph-based analysis.

### 6. Single-Cell Feature Processing

Script:

```text
scripts/08_process_single_cell_features.py
```

Outputs:

```text
results/11_processed_single_cell_features.csv
results/11_processed_single_cell_features_summary.csv
```

Purpose:

- Filter out cells with `area < 4` pixels.
- Censor each marker column to its own 99th percentile.
- Apply arcsinh transformation to marker columns with cofactor 1.
- Preserve the raw exported `data/cells.csv` while saving a processed
  downstream table.

Template-paper alignment:

```text
area filtering
-> marker censoring at the 99th percentile
-> arcsinh transformation with cofactor 1
```

### 7. Processed Feature QC

Script:

```text
scripts/09_qc_processed_single_cell_features.py
```

Outputs:

```text
results/12_processed_single_cell_qc_summary.csv
results/12_processed_marker_qc_summary.csv
figures/01_cells_per_image.png
figures/02_area_distribution_after_filtering.png
figures/03_marker_summary_after_transformation.png
```

Purpose:

- Confirm total processed cell counts.
- Summarize cells per ROI image.
- Verify that the area filter was applied.
- Inspect transformed marker ranges.
- Generate technical QC figures before phenotyping.

### 8. Rule-Based Phenotyping

Script:

```text
scripts/10_rule_based_phenotyping.py
```

Outputs:

```text
results/13_cells_with_phenotypes.csv
results/14_phenotype_composition_by_image.csv
results/15_phenotype_composition_by_category.csv
results/16_rule_based_phenotyping_thresholds.csv
figures/04_phenotype_counts_by_image.png
figures/05_phenotype_composition_by_category.png
```

Purpose:

- Assign transparent first-pass phenotype labels from processed marker values.
- Calculate marker positivity thresholds from selected lineage markers.
- Apply a documented rule hierarchy for approximate cell classes.
- Summarize phenotype composition by ROI image and disease/sample category.

Boundary:

This is a simplified rule-based approximation. It does not reproduce the
template paper's full expert annotation, optimized XGBoost classification, and
FlowSOM review workflow. Labels use the suffix `_like` to indicate approximate
marker-pattern assignments.

### 9. Spatial Phenotype Interactions

Script:

```text
scripts/11_spatial_phenotype_interactions.py
```

Outputs:

```text
results/17_spatial_phenotype_interactions_by_image.csv
results/18_spatial_phenotype_interactions_by_category.csv
results/19_spatial_phenotype_interaction_summary.csv
figures/06_spatial_phenotype_interaction_heatmap.png
```

Purpose:

- Join rule-based phenotype labels to cell-neighbor edge lists.
- Collapse reciprocal directed neighbor edges into unique undirected edges.
- Summarize which approximate phenotypes occur next to each other by ROI image
  and category.
- Provide the raw spatial interaction layer used by enrichment analysis.

### 10. Spatial Phenotype Enrichment

Script:

```text
scripts/12_spatial_phenotype_enrichment.py
```

Outputs:

```text
results/20_spatial_phenotype_enrichment_by_image.csv
results/21_spatial_phenotype_enrichment_by_category.csv
figures/07_spatial_phenotype_enrichment_heatmap.png
```

Purpose:

- Compare observed phenotype-neighbor fractions with abundance-based expected
  fractions.
- Calculate observed-to-expected ratios and log2 enrichment values.
- Identify spatial phenotype pairs that occur more or less often than expected
  from phenotype abundance.

Boundary:

This enrichment is an abundance-normalized exploratory statistic. It is not a
permutation-based spatial null model.

### 11. Plasma-Cell-Like Niche Summary

Script:

```text
scripts/13_plasma_cell_niche_analysis.py
```

Outputs:

```text
results/22_plasma_cell_niche_by_image.csv
results/23_plasma_cell_niche_by_category.csv
results/24_plasma_cell_niche_summary.csv
figures/08_plasma_cell_neighbor_enrichment.png
```

Purpose:

- Extract spatial enrichment rows involving `Plasma_cell_like` cells.
- Summarize enriched and depleted plasma-cell-like neighbor phenotypes.
- Apply a minimum-edge threshold for top-neighbor interpretation so rare
  one-off edges do not dominate the summary.

### 12. Final Interpretation Report

Script:

```text
scripts/14_generate_final_interpretation_report.py
```

Output:

```text
results/25_final_interpretation_summary.md
```

Purpose:

- Combine processed feature summaries, QC outputs, phenotype composition,
  spatial enrichment, and plasma-cell-like niche summaries.
- Generate a reusable Markdown report that updates when the pipeline is rerun on
  another compatible IMC dataset.
- Keep the final conclusion cautious by documenting annotation and statistical
  limitations.

### 13. Complete Workflow Runner

Script:

```text
scripts/15_run_complete_workflow.py
```

Purpose:

- Execute the numbered workflow scripts in order.
- Support complete reruns from raw data or downstream-only reruns from processed
  single-cell exports.
- Provide a dry-run mode that prints the commands without executing them.

## Numbered Results

| Output | Description |
|---|---|
| `results/01_panel.csv` | Numbered panel copy |
| `results/02_images.csv` | Numbered image metadata copy |
| `results/03_preprocessed_tiff_inventory.csv` | Preprocessed TIFF inventory |
| `results/04_mesmer_mask_inventory.csv` | Segmentation mask inventory |
| `results/05_intensity_table_inventory.csv` | Marker intensity table inventory |
| `results/06_regionprops_table_inventory.csv` | Region property table inventory |
| `results/07_neighbors_table_inventory.csv` | Spatial neighbor table inventory |
| `results/08_cells.csv` | Exported single-cell CSV |
| `results/09_cells.h5ad` | Exported AnnData object |
| `results/10_graphml_inventory.csv` | GraphML export inventory |
| `results/11_processed_single_cell_features.csv` | Processed single-cell feature table |
| `results/11_processed_single_cell_features_summary.csv` | Processing summary |
| `results/12_processed_single_cell_qc_summary.csv` | Processed table QC summary |
| `results/12_processed_marker_qc_summary.csv` | Marker-level QC summary |
| `results/13_cells_with_phenotypes.csv` | Processed cells with rule-based phenotype labels |
| `results/14_phenotype_composition_by_image.csv` | Phenotype composition by ROI image |
| `results/15_phenotype_composition_by_category.csv` | Phenotype composition by category |
| `results/16_rule_based_phenotyping_thresholds.csv` | Marker positivity thresholds |
| `results/17_spatial_phenotype_interactions_by_image.csv` | Phenotype-neighbor interactions by ROI image |
| `results/18_spatial_phenotype_interactions_by_category.csv` | Phenotype-neighbor interactions by category |
| `results/19_spatial_phenotype_interaction_summary.csv` | Spatial interaction processing summary |
| `results/20_spatial_phenotype_enrichment_by_image.csv` | Observed-versus-expected enrichment by ROI image |
| `results/21_spatial_phenotype_enrichment_by_category.csv` | Observed-versus-expected enrichment by category |
| `results/22_plasma_cell_niche_by_image.csv` | Plasma-cell-like niche enrichment by ROI image |
| `results/23_plasma_cell_niche_by_category.csv` | Plasma-cell-like niche enrichment by category |
| `results/24_plasma_cell_niche_summary.csv` | Plasma-cell-like niche interpretation summary |
| `results/25_final_interpretation_summary.md` | Final data-driven interpretation report |

## Figures

| Figure | Description |
|---|---|
| `figures/01_cells_per_image.png` | Processed cell counts per image |
| `figures/02_area_distribution_after_filtering.png` | Area distribution after filtering |
| `figures/03_marker_summary_after_transformation.png` | Marker medians and high-end transformed values |
| `figures/04_phenotype_counts_by_image.png` | Rule-based phenotype counts by image |
| `figures/05_phenotype_composition_by_category.png` | Rule-based phenotype fractions by category |
| `figures/06_spatial_phenotype_interaction_heatmap.png` | Spatial phenotype interaction heatmap |
| `figures/07_spatial_phenotype_enrichment_heatmap.png` | Spatial phenotype enrichment heatmap |
| `figures/08_plasma_cell_neighbor_enrichment.png` | Plasma-cell-like neighbor enrichment heatmap |

## Current Status

| Stage | Status |
|---|---|
| Workspace initialization | Complete |
| Raw representative input setup | Complete |
| Panel generation | Complete |
| Steinbock preprocessing | Complete |
| Mesmer segmentation | Complete |
| Intensity measurement | Complete |
| Region property measurement | Complete |
| Neighbor graph construction | Complete |
| CSV / AnnData / GraphML export | Complete |
| Single-cell feature processing | Complete |
| Processed feature QC | Complete |
| Rule-based phenotyping | Complete |
| Spatial phenotype interaction analysis | Complete |
| Spatial phenotype enrichment analysis | Complete |
| Plasma-cell-like niche summary | Complete |
| Final interpretation report | Complete |

## Current Summary

| Metric | Value |
|---|---:|
| Retained ROI images | 8 |
| Exported cells before feature filtering | 50,877 |
| Processed cells after area filtering | 50,849 |
| Cells removed with `area < 4` | 28 |
| Marker columns processed | 42 |
| Rule-based phenotype classes assigned | 14 |
| Unique undirected spatial edges analyzed | 127,090 |
| Final interpretation report | `results/25_final_interpretation_summary.md` |

## Interpretation Boundaries

- The raw-to-image, segmentation, measurement, neighbor, and export stages follow
  a Steinbock-style IMC workflow inspired by the template paper.
- The single-cell feature-processing step implements the template-paper
  description for area filtering, 99th percentile censoring, and arcsinh
  transformation.
- Rule-based phenotyping is a transparent first-pass approximation for the
  representative subset.
- Rule-based labels should not be described as exact reproduction of the
  template paper's final cell annotations.
- Composition tables and phenotype figures are exploratory until validated
  against expert annotation, reference labels, or a supervised annotation model.
- Spatial enrichment is abundance-normalized but not yet permutation-tested.
- The final report is data-driven and reusable, but the biological strength of
  its conclusions depends on phenotype validation and appropriate statistical
  testing for the target dataset.

## Suggested Next Step

The current workflow now reaches an exploratory final interpretation report. The
next methodological improvement should be either reference-guided phenotype
validation or permutation-based spatial enrichment testing. Bone-mask and
distance-to-bone analysis can be added when reproducible bone annotations are
available for the target dataset.
