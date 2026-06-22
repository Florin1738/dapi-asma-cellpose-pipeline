# DAPI Target-Stain Normalization

This repository will contain a reproducible local pipeline for fluorescence microscopy images where DAPI is used to count nuclei and a separate target channel is measured for signal.

The main endpoint is:

```text
target_integrated_intensity_per_DAPI_positive_nucleus =
    target_integrated_background_corrected / filtered_DAPI_positive_nucleus_count
```

This pipeline uses DAPI to count DAPI-positive nuclei. It does not normalize by DAPI fluorescence intensity. The main endpoint is target-channel integrated intensity per DAPI-positive nucleus.

## Current Status

This is a groundwork pass. The repository now captures the scientific contract, install strategy, output schema, validation plan, and implementation roadmap. The Cellpose runtime is not installed yet, model weights are not downloaded yet, and the analysis CLI is not implemented yet.

## Why Cellpose First

Cellpose is the v1 segmentation target because it is a mature Python package, supports pretrained segmentation workflows, has a documented nuclei-oriented path, and can run on both Apple Silicon and Windows/Linux NVIDIA environments through PyTorch.

StarDist remains useful for comparison, especially for compact nuclei, but it brings TensorFlow environment constraints. CellProfiler and QuPath are valuable GUI/reference tools but are not the default runtime path for this folder-based CLI project.

## Planned Command

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
- [Cellpose smoke test](docs/CELLPOSE_SMOKE_TEST.md): local Cellpose setup result and current model-download blocker.
- [Roadmap](docs/ROADMAP.md): implementation sequence.

## Near-Term Milestones

1. Build synthetic image generation and measurement math tests.
2. Implement TIFF/OME-TIFF discovery and channel extraction.
3. Implement Cellpose nuclei segmentation with saved masks and centroids.
4. Export image-level measurements and QC overlays.
5. Add manual centroid validation.
6. Evaluate StarDist only after Cellpose is working and QC-reviewed.
