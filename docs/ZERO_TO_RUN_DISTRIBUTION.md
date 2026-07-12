# Zero-To-Run Distribution Plan

Goal: a new user with no project checkout and no Python/Cellpose knowledge can become operational with the fewest steps.

## Recommended User Flow

1. User downloads the Mac or Windows project zip.
2. User unzips it.
3. User double-clicks setup once:
   - Mac: `Setup Cellpose DAPI aSMA Pipeline.command`
   - Windows: `Setup Cellpose DAPI aSMA Pipeline Windows.cmd`
4. User double-clicks the run launcher:
   - Mac: `Run Cellpose DAPI aSMA Pipeline.command`
   - Windows: `Run Cellpose DAPI aSMA Pipeline Windows.cmd`
5. User selects the image-data folder and output parent folder.
6. User starts with `final/START_HERE_RUN_SUMMARY.html`.

## What Setup Does

The setup launchers create everything locally inside the project folder:

- installs or locates `uv` for the current user account
- installs Python 3.12 into a project `.venv`
- installs this project plus Cellpose into `.venv`
- downloads `cpsam_v2` into `.models/cellpose/cpsam_v2`
- verifies imports
- runs a discovery dry run when bundled plate data are present

No CellProfiler, Docker, conda, Homebrew, admin rights, or global Python install is required.

## Release Folder Contents

Minimum online installer zip:

```text
Setup Cellpose DAPI aSMA Pipeline.command
Setup Cellpose DAPI aSMA Pipeline Windows.cmd
Run Cellpose DAPI aSMA Pipeline.command
Run Cellpose DAPI aSMA Pipeline Windows.cmd
pyproject.toml
constraints/
scripts/
src/
docs/
README.md
```

Optional faster/offline zip additions:

```text
.models/cellpose/cpsam_v2
wheelhouse/
```

If `wheelhouse/` contains compatible wheels, setup installs from it instead of downloading Python packages. The model file is reused when present and checksum-valid.

## Platform Notes

Mac and Windows must create their own `.venv`. Do not copy a Mac `.venv` to Windows or a Windows `.venv` to Mac.

For Windows NVIDIA machines, GPU acceleration depends on the installed PyTorch build and NVIDIA driver. The default setup is still operational without CellProfiler or Docker; a maintainer can preinstall the official CUDA-enabled PyTorch build for that workstation before routine users run data.

## Maintainer Smoke Tests

Mac:

```bash
.venv/bin/python scripts/run_user_cellpose_batch.py \
  --input data/aSMA_DAPI_plates \
  --output output/user_cellpose_batch/macos_dry_run_placeholder \
  --dry-run
```

Windows:

```powershell
.\.venv\Scripts\python.exe scripts\run_user_cellpose_batch.py `
  --input data\aSMA_DAPI_plates `
  --output output\user_cellpose_batch\windows_dry_run_placeholder `
  --dry-run
```

Expected current discovery:

```text
4 acquisition folders
227 image pairs
```

Then run a one-image smoke test with `--max-images-per-acquisition 1 --skip-figures` before handing the zip to routine users.
