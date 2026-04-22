# Workflow Status Log

This file records the stepwise status of the representative IMC study
reproduction. Each entry documents the input, action, output, and validation
status for a completed step.

## Step 1: Workspace Initialization

Status: Complete

Action:

- Created a dedicated reproduction workspace at `step_by_step_reproduction/`.
- Created subdirectories for raw data, notebooks, scripts, logs, results, and
  figures.

Output:

```text
step_by_step_reproduction/
  data/raw/
  notebooks/
  scripts/
  logs/
  results/
  figures/
```

Interpretation:

- The reproduction has an isolated working directory.
- All generated files should remain inside this directory.

## Step 2: Raw Input Data Check

Status: Complete

Input:

```text
data/raw/
```

Observed files:

| Case | Category | Raw `.mcd` | ROI text files |
|---|---|---:|---:|
| `TS-373_IMC01_UB` | UB | 1 | 2 |
| `TS-373_IMC02_MGUS` | MGUS | 1 | 2 |
| `TS-373_IMC05_MGUS` | MGUS | 1 | 2 |
| `TS-373_IMC09_B` | B | 1 | 2 |

Validation:

- Each selected case has one `.mcd` file.
- Each selected case has two ROI text files.
- The folder is suitable for Steinbock preprocessing.

## Step 3: Panel File Generation

Status: Complete

Input:

```text
data/raw/TS-373_IMC01_UB_ROI_1.txt
```

Action:

- Parsed the ROI text-file header to extract IMC marker/channel names.
- Created a Steinbock-compatible panel file with a reusable script.

Script:

```text
scripts/01_create_panel_from_raw.py
```

Notebook:

```text
notebooks/01_template_study_reproduction_workflow.ipynb
```

Output:

```text
data/panel.csv
results/01_panel.csv
```

Validation:

| Check | Result |
|---|---:|
| ROI text files checked | 8 |
| Raw channels detected | 42 |
| Panel rows created | 42 |
| Missing channels from panel | 0 |
| Extra channels in panel | 0 |

Segmentation marker assignment:

| Marker role | Channel | Marker | `deepcell` value |
|---|---|---|---:|
| Nuclear marker | `Yb171` | `HistoneH3` | 1 |
| Membrane marker | `Yb173` | `CD98` | 2 |

Interpretation:

- The panel file matches the raw data channel structure.
- The panel file is ready for Steinbock preprocessing and Mesmer segmentation.

## Next Planned Step

Step 4: Preprocess raw IMC data with Steinbock.

Command:

```bash
steinbock preprocess imc images --hpf 50
```

Expected outputs:

```text
data/img/*.tiff
data/images.csv
logs/02_preprocess_imc_images.log
results/02_images.csv
results/03_preprocessed_tiff_inventory.csv
```

Script:

```text
scripts/02_preprocess_imc_images.py
```

Notebook:

```text
notebooks/01_template_study_reproduction_workflow.ipynb
```

## Planned Step 5: Mesmer Cell Segmentation

Status: Ready after preprocessing

Command:

```bash
steinbock segment deepcell --app mesmer --minmax
```

Expected outputs:

```text
data/masks/*.tiff
logs/03_segment_mesmer.log
results/04_mesmer_mask_inventory.csv
```

Script:

```text
scripts/03_segment_mesmer.py
```

Notebook:

```text
notebooks/01_template_study_reproduction_workflow.ipynb
```

Interpretation:

- This step implements the intended Mesmer segmentation command.
- If DeepCell model access is unavailable, the command log should be retained as
  documentation of the access limitation.

## Planned Step 6: Marker Intensity Measurement

Status: Ready after segmentation

Command:

```bash
steinbock measure intensities
```

Expected outputs:

```text
data/intensities/*.csv
logs/04_measure_intensities.log
results/05_intensity_table_inventory.csv
```

Script:

```text
scripts/04_measure_intensities.py
```

Notebook:

```text
notebooks/01_template_study_reproduction_workflow.ipynb
```

Interpretation:

- This step measures mean marker intensity for each segmented cell and each IMC
  channel.

## Planned Step 7: Region Property Measurement

Status: Ready after segmentation

Command:

```bash
steinbock measure regionprops
```

Expected outputs:

```text
data/regionprops/*.csv
logs/05_measure_regionprops.log
results/06_regionprops_table_inventory.csv
```

Script:

```text
scripts/05_measure_regionprops.py
```

Notebook:

```text
notebooks/01_template_study_reproduction_workflow.ipynb
```

Interpretation:

- This step measures per-cell geometric properties, including area and centroid
  coordinates.

## Planned Step 8: Spatial Neighbor Measurement

Status: Ready after segmentation

Command:

```bash
steinbock measure neighbors --type expansion --dmax 4
```

Expected outputs:

```text
data/neighbors/*.csv
logs/06_measure_neighbors.log
results/07_neighbors_table_inventory.csv
```

Script:

```text
scripts/06_measure_neighbors.py
```

Notebook:

```text
notebooks/01_template_study_reproduction_workflow.ipynb
```

Interpretation:

- This step creates spatial cell-neighbor edge lists by expanding cell masks up
  to 4 pixels and connecting cells whose expanded regions touch.

## Planned Step 9: Data Export

Status: Ready after intensity, region property, and neighbor measurements

Commands:

```bash
steinbock export csv intensities regionprops -o cells.csv
steinbock export anndata --intensities intensities --data regionprops --neighbors neighbors -o cells.h5ad
steinbock export graphs --format graphml --data intensities --data regionprops
```

Expected outputs:

```text
data/cells.csv
data/cells.h5ad
data/graphs/*.graphml
logs/07_export_data.log
results/08_cells.csv
results/09_cells.h5ad
results/10_graphml/*.graphml
results/10_graphml_inventory.csv
```

Script:

```text
scripts/07_export_data.py
```

Notebook:

```text
notebooks/01_template_study_reproduction_workflow.ipynb
```

Interpretation:

- This step creates analysis-ready single-cell and spatial graph outputs.
- Steinbock-required working filenames remain in `data/`.
- Numbered research copies are stored in `results/`.

## Planned Step 10: Workflow Output Validation

Status: Ready after data export

Expected output:

```text
results/11_workflow_output_validation.csv
```

Script:

```text
scripts/08_validate_workflow_outputs.py
```

Notebook:

```text
notebooks/01_template_study_reproduction_workflow.ipynb
```

Interpretation:

- This step checks whether expected workflow files exist.
- This step is a file-completion validation, not a scientific quality-control
  analysis.

## Planned Step 11: Clean Downstream Analysis

Status: Ready in clean downstream notebook

Expected outputs:

```text
results/12_cell_qc_summary.csv
results/13_cells_with_phenotypes.csv
results/14_phenotype_composition_by_category.csv
results/15_phenotype_spatial_interactions.csv
results/16_phenotype_spatial_enrichment.csv
figures/01_cell_counts_by_category.png
figures/02_phenotype_composition_by_category.png
figures/03_spatial_enrichment_heatmap.png
```

Script:

```text
scripts/09_clean_downstream_analysis.py
```

Notebook:

```text
notebooks/02_clean_downstream_analysis.ipynb
```

Interpretation:

- This step performs compact QC, rule-based phenotyping, phenotype composition,
  and phenotype-aware spatial interaction analysis.
- This step is a clean representative-subset downstream analysis, not exact
  reproduction of the full model-based phenotyping/spatial workflow.
