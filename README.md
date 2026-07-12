# Cellpose DAPI / aSMA Per-Nucleus Quantification

Self-contained local tool for measuring alpha-smooth muscle actin (aSMA) signal relative to DAPI-positive nucleus count in fluorescence microscopy images.

The primary endpoint is:

```text
target_integrated_intensity_per_DAPI_positive_nucleus =
    target_integrated_intensity / filtered_DAPI_positive_nucleus_count
```

DAPI is used to segment and count DAPI-positive nuclei. This is not DAPI fluorescence normalization: DAPI brightness is not used as the denominator.

## Quick Start For Mac And Windows Users

The easiest way to try the tool is to download the release zip for your operating system from the GitHub Releases page.

1. Download one zip:
   - Mac: `cellpose-dapi-asma-pipeline-mac-*.zip`
   - Windows: `cellpose-dapi-asma-pipeline-windows-*.zip`
2. Unzip it.
3. Double-click the run launcher:
   - Mac: `Run Cellpose DAPI aSMA Pipeline.command`
   - Windows: `Run Cellpose DAPI aSMA Pipeline Windows.cmd`
4. On first launch, wait while the tool installs its local Python environment, Cellpose, and the Cellpose model. This needs an internet connection and can take several minutes.
5. In the app window, choose the image folder, choose where results should go, and press **Run analysis**.

You do not need to install Python, Cellpose, CellProfiler, Docker, conda, or plugins yourself. The launcher creates a project-local `.venv/` and a project-local `.models/cellpose/cpsam_v2` model cache.

### macOS Security Prompt

On first launch, macOS may block the `.command` file because it was downloaded from the internet. Use:

```text
Right-click or Control-click the launcher -> Open -> Open
```

You normally only need to do that once.

### Windows Notes

The Windows launcher uses PowerShell internally and creates `.venv\Scripts\python.exe` inside the unzipped folder. It does not need administrator rights for the project environment. A Windows NVIDIA machine can be configured for GPU acceleration, but the default setup remains a local Cellpose install either way.

## Expected Image Folder

The routine runner expects TIFF image pairs organized by acquisition and field:

```text
Parent folder/
  AcquisitionName/
    XY01/*_CH2.tif
    XY01/*_CH4.tif
    XY02/*_CH2.tif
    XY02/*_CH4.tif
```

Default channel mapping:

- `CH2`: aSMA target channel
- `CH4`: DAPI nuclei channel

Use the app's channel selectors if your dataset uses different channel names. The tool supports TIFF/OME-TIFF first. Proprietary formats such as `.czi`, `.nd2`, and `.lif` should be converted before routine use unless an explicit converter path is added.

## Main Outputs

Each run creates a timestamped output folder. Start with:

```text
final/START_HERE_RUN_SUMMARY.html
final/workbooks/cellpose_background_corrected_pi_style_summary.xlsx
final/tables/cellpose_user_friendly_endpoint_summary.csv
final/tables/cellpose_full_plate_endpoint_summary.csv
final/figures/
logs/user_cellpose_batch_run.yaml
```

Every analysis run is expected to produce visual QC overlays before results are considered reviewable. Visual QC is qualitative only; precision, recall, F1, false positives, and false negatives require manual ground truth.

## Scientific Scope

- The pipeline computes auditable image-analysis metrics.
- It does not interpret biological significance.
- Cellpose is the default backend.
- In user-facing Cellpose outputs, "Cellpose" means Cellpose objects retained only when at least one DAPI-positive nucleus centroid falls inside the object.
- Current cell-area style masks are candidate aSMA-associated regions, not validated whole-cell masks.

## Developer Setup From Git

Clone the repository, then create a project-local environment. On macOS and Linux, `uv` with Python 3.12 is preferred:

```bash
uv python install 3.12
uv venv --python 3.12 .venv
.venv/bin/python -m pip install -e '.[dev]'
```

To install the Cellpose runtime too:

```bash
.venv/bin/python -m pip install -e '.[dev,cellpose]'
```

Run tests:

```bash
.venv/bin/python -m pytest
```

Build release zips:

```bash
.venv/bin/python scripts/make_release.py --version v1.0
```

The release builder excludes real microscopy images, generated outputs, masks, model weights, manual annotations, and report artifacts.

## Command-Line Runner

The double-click app delegates to the same batch runner:

```bash
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_user_cellpose_batch.py \
  --input /path/to/image_parent_folder \
  --output /path/to/results_parent_folder \
  --target-channel CH2 \
  --dapi-channel CH4
```

For a quick smoke test on a large folder:

```bash
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_user_cellpose_batch.py \
  --input /path/to/one_acquisition_folder \
  --output /path/to/results_parent_folder \
  --target-channel CH2 \
  --dapi-channel CH4 \
  --max-images-per-acquisition 1 \
  --skip-figures
```

## Repository Data Policy

This repository is intended to contain source code, tests, configuration, setup launchers, and documentation. It must not contain:

- real microscopy images
- generated masks, overlays, tables, workbooks, or reports derived from private lab data
- Cellpose model weights
- lab-private annotations
- virtual environments

Folder-level README files are kept in `data/`, `output/`, `manual_validation/`, and `reports/` so local users understand those directories, but their contents are ignored by git.

## Documentation

- [Installation](docs/INSTALLATION.md): Mac, Windows, and developer setup.
- [Nontechnical Cellpose batch runner](docs/NONTECHNICAL_CELLPOSE_BATCH_RUNNER.md): double-click workflow and expected outputs.
- [Zero-to-run distribution](docs/ZERO_TO_RUN_DISTRIBUTION.md): release zip contents and maintainer smoke tests.
- [Release packaging](docs/RELEASE_PACKAGING.md): maintainer steps for building and uploading GitHub release zips.
- [Output contract](docs/OUTPUT_CONTRACT.md): stable output columns and file layout.
- [Image inputs](docs/IMAGE_INPUTS.md): supported inputs, channel selection, and axis policy.
- [QC and visualization](docs/QC_AND_VISUALIZATION.md): visual QC expectations.
- [Validation](docs/VALIDATION.md): manual centroid and mask validation protocols.
- [Methodology](docs/methodology/README.md): accepted, exploratory, and rejected method records.
