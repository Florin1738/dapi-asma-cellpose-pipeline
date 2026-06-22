# DAPI Target-Stain Normalization

This repository will contain a reproducible local pipeline for fluorescence microscopy images where DAPI is used to count nuclei and a separate target channel is measured for signal.

The main endpoint is:

```text
target_integrated_intensity_per_DAPI_positive_nucleus =
    target_integrated_background_corrected / filtered_DAPI_positive_nucleus_count
```

This pipeline uses DAPI to count DAPI-positive nuclei. It does not normalize by DAPI fluorescence intensity. The main endpoint is target-channel integrated intensity per DAPI-positive nucleus.

## Current Status

The Cellpose-first groundwork is now runnable for nuclei counting on the local sample data. The project has:

- a project-local Python environment with Cellpose 4.2.1.1 and PyTorch
- a local `cpsam_v2` model cache in `.models/cellpose`
- data inventory and channel-preview tooling for the `ApYYM20AGGSMA_02` folder
- a count-only Cellpose CLI for the candidate DAPI channel
- per-position masks, CSV summaries, and QC montages for all 12 sample positions

The full target-channel intensity normalization endpoint is still future work. Also, the sample files are RGB pseudocolor TIFF exports, so `CH4` is treated as candidate DAPI until microscope metadata or acquisition notes confirm the channel identity.

## Why Cellpose First

Cellpose is the v1 segmentation target because it is a mature Python package, supports pretrained segmentation workflows, has a documented nuclei-oriented path, and can run on both Apple Silicon and Windows/Linux NVIDIA environments through PyTorch.

StarDist remains useful for comparison, especially for compact nuclei, but it brings TensorFlow environment constraints. CellProfiler and QuPath are valuable GUI/reference tools but are not the default runtime path for this folder-based CLI project.

## Current Count Command

The current implemented command counts nuclei from the candidate DAPI channel:

```bash
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_cellpose_counts.py \
  --input ApYYM20AGGSMA_02 \
  --output output/cellpose_counts \
  --channel CH4 \
  --model cpsam_v2 \
  --gpu
```

PowerShell equivalent for a Windows NVIDIA setup:

```powershell
$env:CELLPOSE_LOCAL_MODELS_PATH = "$PWD\.models\cellpose"
python scripts\run_cellpose_counts.py `
  --input ApYYM20AGGSMA_02 `
  --output output\cellpose_counts `
  --channel CH4 `
  --model cpsam_v2 `
  --gpu
```

Current output summary:

```text
output/cellpose_counts/summaries/nucleus_counts.csv
output/cellpose_counts/summaries/per_nucleus_locations.csv
output/cellpose_counts/masks/
output/cellpose_counts/qc/
output/cellpose_counts/qc_contact_sheet.png
output/cellpose_counts/logs/config_resolved.yaml
output/cellpose_counts/logs/run_log.txt
```

Validate the generated count artifacts:

```bash
.venv/bin/python scripts/validate_cellpose_counts.py --output output/cellpose_counts
```

Windows PowerShell:

```powershell
python scripts\validate_cellpose_counts.py --output output\cellpose_counts
```

## Current Target Normalization Command

The current provisional target-normalization command measures full-field candidate `CH2` signal and divides by the candidate `CH4` Cellpose nucleus count:

```bash
.venv/bin/python scripts/run_target_normalization.py \
  --input ApYYM20AGGSMA_02 \
  --counts output/cellpose_counts \
  --output output/target_normalization \
  --target-channel CH2 \
  --dapi-channel CH4 \
  --background-percentile 10
```

Validate the generated normalized-intensity artifacts:

```bash
.venv/bin/python scripts/validate_target_normalization.py --output output/target_normalization
```

Current target-normalization outputs:

```text
output/target_normalization/summaries/image_level_summary.csv
output/target_normalization/summaries/well_level_summary.csv
output/target_normalization/plots/normalized_intensity_by_well.png
output/target_normalization/plots/target_integrated_vs_nucleus_count.png
output/target_normalization/qc/
output/target_normalization/qc_contact_sheet.png
output/target_normalization/logs/config_resolved.yaml
```

## Planned Full Pipeline Command

```bash
python scripts/run_pipeline.py \
  --input "/path/to/images" \
  --output "/path/to/output" \
  --dapi-channel 0 \
  --target-channel 1 \
  --backend cellpose \
  --diameter auto \
  --background percentile \
  --background-percentile 10 \
  --save-qc
```

## Planned Outputs

```text
output/
  logs/
    run_log.txt
    config_resolved.yaml
  summaries/
    image_level_summary.csv
    per_nucleus_locations.csv
    roi_level_summary.csv
    validation_metrics.csv
  masks/
    <image_id>_nuclei_labels.tif
  qc/
    <image_id>_dapi_nucleus_outlines.png
    <image_id>_target_with_nucleus_outlines.png
    <image_id>_dapi_numbered_centroids.png
    <image_id>_segmentation_montage.png
  validation/
    <image_id>_validation_overlay.png
```

## Documentation Map

- [Software research](docs/SOFTWARE_RESEARCH.md): Cellpose, StarDist, CellProfiler, QuPath, and platform choices.
- [Installation](docs/INSTALLATION.md): Mac Apple Silicon and Windows NVIDIA setup plans.
- [Image inputs](docs/IMAGE_INPUTS.md): supported formats, channel and Z-stack policy.
- [Output contract](docs/OUTPUT_CONTRACT.md): CSV columns and file layout.
- [QC and visualization](docs/QC_AND_VISUALIZATION.md): overlays, montages, and how outputs will be reviewed.
- [Validation](docs/VALIDATION.md): manual centroid and mask validation protocol.
- [Data inventory](docs/DATA_INVENTORY.md): current local sample-data shape, without committing private images.
- [Data organization](docs/DATA_ORGANIZATION.md): detailed map of the current `ApYYM20AGGSMA_02` folder and channel-role caveats.
- [Data QC review](docs/DATA_QC_REVIEW.md): representative channel-preview images for visual logic checks.
- [Cellpose smoke test](docs/CELLPOSE_SMOKE_TEST.md): one-image Cellpose setup and segmentation result.
- [Cellpose nuclei count run](docs/CELLPOSE_NUCLEI_COUNT_RUN.md): all-position count-only Cellpose run, output paths, counts, and caveats.
- [Target normalization run](docs/TARGET_NORMALIZATION_RUN.md): provisional CH2-per-CH4-nucleus measurements, visualizations, and caveats.
- [Roadmap](docs/ROADMAP.md): implementation sequence.

## Near-Term Milestones

1. Confirm channel identities from acquisition metadata or raw microscope exports.
2. Add target-channel intensity measurement using the Cellpose nuclei masks.
3. Add area and border filtering to produce raw and filtered nucleus counts.
4. Add manual centroid validation against reviewer annotations.
5. Generate an HTML QC report for batch review.
6. Evaluate StarDist only after Cellpose count QC is accepted.
