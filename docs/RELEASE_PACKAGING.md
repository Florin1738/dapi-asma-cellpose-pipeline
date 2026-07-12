# Release Packaging

This document is for maintainers publishing a GitHub release for routine Mac
and Windows users.

## Build Release Zips

From a developer environment:

```bash
.venv/bin/python scripts/make_release.py --version v1.0
```

Expected outputs:

```text
dist/cellpose-dapi-asma-pipeline-mac-v1.0.zip
dist/cellpose-dapi-asma-pipeline-windows-v1.0.zip
```

These are online installers. First launch downloads Python packages and the
Cellpose model, so the user's machine needs internet access.

For a larger zip that includes the Cellpose model:

```bash
.venv/bin/python scripts/make_release.py --version v1.0 --with-model
```

Use `--with-wheelhouse` only when a compatible `wheelhouse/` has already been
prepared for the target operating system.

## Publish A GitHub Release

1. Commit and push the source code.
2. Create a GitHub release tag such as `v1.0`.
3. Upload both zip files from `dist/` as release assets.
4. In the release notes, tell users to download the zip for their operating
   system, unzip it, and double-click the `Run Cellpose DAPI aSMA Pipeline`
   launcher.

Do not commit the zip files to git. They are release assets.

## Smoke Tests Before Sharing

Check that the zips contain only distributable files:

```bash
.venv/bin/python -m pytest tests/test_make_release.py tests/test_setup_launchers.py tests/test_gui_launchers.py -q
```

Optional Mac dry run when local sample data are present:

```bash
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_user_cellpose_batch.py \
  --input data/aSMA_DAPI_plates \
  --output output/user_cellpose_batch/macos_dry_run \
  --target-channel CH2 \
  --dapi-channel CH4 \
  --dry-run
```

Optional one-image smoke test:

```bash
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_user_cellpose_batch.py \
  --input /path/to/one_acquisition_folder \
  --output output/user_cellpose_batch/smoke \
  --target-channel CH2 \
  --dapi-channel CH4 \
  --max-images-per-acquisition 1 \
  --skip-figures
```

## Privacy Check

Release zips and git commits must not include:

- real microscopy images
- generated masks, overlays, tables, workbooks, or reports derived from private
  lab data
- Cellpose model weights unless explicitly building an offline release asset
- manual annotation artifacts
- virtual environments
