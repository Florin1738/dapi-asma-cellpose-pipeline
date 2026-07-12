# Installation And Runtime Setup

This project has two setup paths:

- Routine users download a release zip and double-click a launcher.
- Developers clone the repository and create a project-local Python environment.

Do not install into system Python. Do not commit `.venv/`, Cellpose model weights, real microscopy images, or generated outputs.

## Routine Mac Install

1. Download `cellpose-dapi-asma-pipeline-mac-*.zip` from the GitHub Releases page.
2. Unzip it.
3. Double-click `Run Cellpose DAPI aSMA Pipeline.command`.
4. On first launch, wait for setup to finish. It installs `uv` for the current user if needed, creates `.venv/`, installs the project and Cellpose, downloads `cpsam_v2`, verifies imports, and opens the app.
5. Pick the image folder and results folder in the app, then press **Run analysis**.

If macOS blocks the launcher, right-click or Control-click it, choose **Open**, then confirm **Open**. This is the usual Gatekeeper prompt for an unsigned downloaded script.

The optional `Setup Cellpose DAPI aSMA Pipeline.command` launcher performs the same setup without starting a run. Use it only to preinstall or repair the environment.

## Routine Windows Install

1. Download `cellpose-dapi-asma-pipeline-windows-*.zip` from the GitHub Releases page.
2. Unzip it.
3. Double-click `Run Cellpose DAPI aSMA Pipeline Windows.cmd`.
4. On first launch, wait for setup to finish. It installs `uv` for the current user if needed, creates `.venv\Scripts\python.exe`, installs the project and Cellpose, downloads `cpsam_v2`, verifies imports, and opens the app.
5. Pick the image folder and results folder in the app, then press **Run analysis**.

The optional `Setup Cellpose DAPI aSMA Pipeline Windows.cmd` launcher performs the same setup without starting a run. Use it only to preinstall or repair the environment.

For Windows NVIDIA machines, GPU acceleration depends on the installed NVIDIA driver and PyTorch build. The default release setup is still operational without CellProfiler, Docker, conda, or admin rights. A maintainer can install the machine-appropriate CUDA-enabled PyTorch wheel before routine users run large batches.

## Input Requirements

The default workflow expects TIFF image pairs organized by acquisition and field:

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

Use the app channel selectors when a dataset has different channel names. If axes or channels cannot be inferred, the pipeline should fail clearly rather than silently guessing.

## Developer Setup

Use Python 3.12 in a project-local environment. On macOS, `uv` is preferred:

```bash
uv python install 3.12
uv venv --python 3.12 .venv
.venv/bin/python -m pip install -e '.[dev]'
```

Install Cellpose runtime dependencies when you need to run segmentation:

```bash
.venv/bin/python -m pip install -e '.[dev,cellpose]'
```

Optional Mac constraints from the original pilot environment:

```bash
.venv/bin/python -m pip install -c constraints/cellpose-mac-2026-06-22.txt -e '.[dev,cellpose]'
```

Run tests:

```bash
.venv/bin/python -m pytest
```

Run the batch CLI:

```bash
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_user_cellpose_batch.py \
  --input /path/to/image_parent_folder \
  --output /path/to/results_parent_folder \
  --target-channel CH2 \
  --dapi-channel CH4
```

## Model Cache

Cellpose downloads model weights on first use. This project uses a project-local model cache:

```text
.models/cellpose/cpsam_v2
```

The expected `cpsam_v2` SHA-256 recorded in `scripts/pipeline_env.py` is:

```text
0f1cc3f7ecdd8a037a57c6c48d9d8921391be4cbce3fa9f13c3e3a2e1253c667
```

Model weights are not committed to git.

## Release Packaging

Build online installer zips:

```bash
.venv/bin/python scripts/make_release.py --version v1.0
```

Outputs:

```text
dist/cellpose-dapi-asma-pipeline-mac-v1.0.zip
dist/cellpose-dapi-asma-pipeline-windows-v1.0.zip
```

These zips include the launchers, source, scripts, docs, configs, examples, and tests. They exclude private data, generated outputs, model weights, virtual environments, masks, figures, workbooks, and report artifacts.

For an offline or faster first run, maintainers may add:

```bash
.venv/bin/python scripts/make_release.py --version v1.0 --with-model
```

Use `--with-wheelhouse` only when a compatible local `wheelhouse/` has already been prepared.

## Source Links

- Cellpose installation docs: https://cellpose.readthedocs.io/en/latest/installation.html
- Cellpose model docs: https://cellpose.readthedocs.io/en/latest/models.html
- PyTorch local install selector: https://pytorch.org/get-started/locally/
