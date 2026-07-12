# AGENTS.md

Project-level instructions for Codex and other coding agents working in this repository.

## Core Scientific Contract

- The primary endpoint is `target_integrated_intensity_per_DAPI_positive_nucleus`.
- DAPI is used to segment/count DAPI-positive nuclei. Do not describe the output as "DAPI-normalized intensity" because that implies division by DAPI fluorescence brightness.
- The normalization denominator is the filtered count of DAPI-positive nuclei in the image or ROI.
- Do not use DAPI intensity itself for normalization unless a future user request explicitly changes the biological question.
- Do not interpret biological significance from the measurement. This repository computes auditable image-analysis metrics.
- Current user-facing Cellpose retained-region endpoint is DAPI-anchored by default: in reports, figures, and workbooks, "Cellpose" means Cellpose objects retained only if at least one DAPI-positive nucleus centroid falls inside the object.
- Do not present anchored-vs-nonanchored Cellpose comparisons in user-facing outputs unless a future user explicitly requests a method-comparison analysis.
- User-facing citations must use this format unless explicitly changed: first three authors, last author, PMID. Do not include journal or year in the citation string unless the user asks for it.

## Implementation Priorities

- v1 backend: Cellpose, using a pretrained nuclei-capable model.
- StarDist is a secondary comparison backend, not required for the first runnable version.
- CellProfiler and QuPath are reference/GUI workflows, not default runtime dependencies.
- Use project-local environments. Do not install into global Python.
- On this Apple Silicon Mac, prefer `uv` with Python 3.12 for local development.
- On a Windows NVIDIA workstation, use a Python 3.12 virtual environment plus the official PyTorch CUDA install command for that machine.

## Data and File Handling

- Never commit real microscopy images, outputs, masks, model weights, or lab-private annotations.
- Support TIFF/OME-TIFF first. Treat `.czi`, `.nd2`, and `.lif` as optional/proprietary inputs that may require conversion.
- Do not silently guess ambiguous dimensions. If axes cannot be inferred, fail with a clear message and ask for explicit axis/channel flags.
- Record interpreted axes, channel selection, Z projection mode, segmentation parameters, and background settings in run logs.

## Quality Control and Validation

- Every analysis run must produce visual QC overlays before results are considered reviewable.
- Visual inspection is qualitative QC only. Do not report precision, recall, F1, false positives, or false negatives unless manual ground truth exists.
- Manual centroid validation must use a configurable pixel radius and must export true positives, false positives, false negatives, precision, recall, F1, and count error.
- Manual mask validation must use an IoU threshold and must export matching metrics plus mean IoU for matched objects.

## Methodology Logging Contract

- Methodology documentation lives in `docs/methodology/`.
- Every new analysis method, threshold, segmentation parameter, QC rule, or source-backed scientific assumption must be logged in `docs/methodology/run_logs/YYYY-MM-DD-<short-topic>.md`.
- Each methodology log must record: purpose, input data, command or script, output paths, exact parameters, threshold rationale, sources consulted, QC visuals, validation status, and whether the method is accepted, exploratory, or rejected.
- Do not promote exploratory diagnostics into default analysis outputs without updating `docs/methodology/CURRENT_METHODS.md` and `docs/methodology/THRESHOLDS_AND_PARAMETERS.md`.
- Current cell-area status: true whole-cell segmentation is not validated from the available DAPI+aSMA-only images. DAPI-only spatial expansion/territory masks are exploratory and must not be described as true cell masks.
- CellProfiler-style secondary-object workflows may be evaluated as reference workflows, but using CH2/aSMA itself as the cell-boundary channel risks biasing segmentation toward aSMA-high cells.

## Engineering Rules

- Keep modules small and named by responsibility: image I/O, channel extraction, segmentation, measurement, QC rendering, validation, CLI.
- Prefer reproducible CLI commands and config files over notebook-only workflows.
- Add synthetic test data early so the pipeline can be tested without private microscopy data.
- Use explicit, stable CSV column names from `docs/OUTPUT_CONTRACT.md`.
- Run tests or document why tests could not be run before claiming functionality is complete.
