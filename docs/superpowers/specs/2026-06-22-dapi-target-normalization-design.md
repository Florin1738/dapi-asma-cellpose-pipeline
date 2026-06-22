# DAPI Target Normalization Design

Date: 2026-06-22

## Goal

Build a reproducible local pipeline that counts DAPI-positive nuclei, measures target-channel signal, corrects target background, and reports target integrated intensity per DAPI-positive nucleus.

## Scientific Contract

The endpoint is:

```text
target_integrated_intensity_per_DAPI_positive_nucleus =
    target_integrated_background_corrected / filtered_DAPI_positive_nucleus_count
```

DAPI is a nuclear counterstain for cell/nucleus counting. The pipeline must not divide by DAPI fluorescence intensity and must not call the endpoint "DAPI-normalized intensity."

## Primary Architecture

The pipeline will be a Python CLI that processes image folders and writes auditable outputs. The first backend is Cellpose. Each run writes the resolved configuration, logs, CSV summaries, label masks, and QC overlays.

Modules will be split by responsibility:

- image discovery and loading
- axis/channel extraction
- segmentation backend adapter
- object filtering
- target measurement and background correction
- QC rendering
- validation
- CLI/config orchestration

## Platform Strategy

This repository should work on:

- Apple Silicon Mac for development, small data, synthetic tests, and pilot runs.
- Windows NVIDIA workstation for larger batch processing.

The Mac setup uses `uv` and Python 3.12. The Windows setup uses Python 3.12, a virtual environment, and the official PyTorch CUDA wheel selected for the machine.

## Inputs

v1 supports TIFF and OME-TIFF. Proprietary formats are documented as optional and should not block the first version. Ambiguous axes must stop execution with a clear message.

## Outputs

The required output contract is defined in `docs/OUTPUT_CONTRACT.md`. CSV columns must use stable names, especially `target_integrated_intensity_per_DAPI_positive_nucleus`.

## QC and Validation

Every run must save QC overlays. Precision/recall/F1 are available only when manual centroid or manual mask ground truth is supplied. Visual QC must not be represented as quantitative validation.

## Non-Goals for v1

- No web app.
- No biological interpretation.
- No required CellProfiler, QuPath, Docker, StarDist, or proprietary image-reader dependency.
- No DAPI-intensity normalization.
- No whole-cell segmentation from DAPI alone.

## Approval State

The user approved a Cellpose-first groundwork pass and asked for setup, visualization, validation, and project documentation to be prepared before heavy installs.

