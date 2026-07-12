# Nontechnical Cellpose Batch Runner

This is the lowest-friction production path for a user who should not have to manage Cellpose commands, CellProfiler plugins, Docker, or postprocessing scripts.

## New User From Zero

Give the user one project zip/folder for their operating system. They unzip it and double-click a single file. There is no separate setup step to remember: the first launch installs everything automatically, then opens the app. Later launches open the app immediately.

The folder must contain:

- Mac: `Run Cellpose DAPI aSMA Pipeline.command`  (the one to double-click)
- Windows: `Run Cellpose DAPI aSMA Pipeline Windows.cmd`  (the one to double-click)
- Manual-repair fallback: `Setup Cellpose DAPI aSMA Pipeline.command` / `... Windows.cmd`
- `scripts/`, `src/`, and `docs/`

Expected user burden from zero:

1. Download the project zip for Mac or Windows.
2. Unzip it.
3. Double-click the run launcher.
4. On the first launch only, wait a few minutes while it installs itself (needs internet).
5. In the window that opens, pick the image folder and the results folder, then press **Run analysis**.

The first launch creates `.venv/` and downloads `.models/cellpose/cpsam_v2`. A Mac `.venv/` is not portable to Windows; each OS builds its own. The install needs internet unless the maintainer bundles a local `wheelhouse/` and `.models/cellpose/cpsam_v2` in the zip. It does not require CellProfiler, Docker, conda, Homebrew, or admin rights.

## What The User Does

1. Double-click the run launcher in the project folder. On the first launch it installs itself (a few minutes), then the app window opens automatically.
2. In the app: click **Browse…** next to *Image folder* and pick the plate or acquisition data.
3. Click **Browse…** next to *Results go to* and pick where results should be saved.
4. Press **Run analysis** and watch the progress log fill in.
5. When it finishes, the results folder opens automatically — start with `START_HERE_RUN_SUMMARY.html` and the Excel workbook in `final/workbooks/`.

The app scans the image folder as soon as it is chosen and shows how many acquisition folders and image pairs it found, so the user can confirm the right folder before running. Scientific parameters are hidden behind a collapsible **▸ Advanced options** panel and default to the validated values; hovering any option shows a plain-language description. Most users never open it.

The input folder can be a full multi-plate folder or one acquisition folder. The runner discovers acquisition folders by finding the folder one level above `XY##` image folders, so plate names, drug-folder names, and acquisition names can vary. The stable requirement is that each acquisition contains `XY##` folders with matching channel TIFFs.

For the July 2026 download shape, the user can either pick the top extracted folder once or run each acquisition folder separately. The four acquisition folders are the folders named like:

```text
.../Drug 1-3/AmiYM0236Agg523CFBSMA/
.../Drug 1-3/AmiYM0236Agg523CFBSMA_01/
.../Drug 4-6/APYPIKfaldstatCFBSMA/
.../Drug 4-6/APYPIKfaldstatCFBSMA_01/
```

Each of those folders directly contains many `XY##` folders. Running one acquisition folder creates one Excel workbook for that acquisition. Running the top folder creates one mirrored output folder and one Excel workbook per discovered acquisition.

The general expected shape is:

```text
Any parent folders/
  AcquisitionName/
    XY01/*_CH2.tif
    XY01/*_CH4.tif
    XY02/*_CH2.tif
    XY02/*_CH4.tif
```

The app has channel selectors for **aSMA target channel** and **DAPI nuclei channel**. Defaults are `CH2` for aSMA target and `CH4` for DAPI, but these must be changed when a dataset uses a different channel order. For the July 2026 download, visual channel review found `CH4 = DAPI`, `CH2 = aSMA target`, and `CH1 = green ECM not analyzed`.

## First Run On Mac

Double-click `Run Cellpose DAPI aSMA Pipeline.command`. On the first launch a Terminal window opens and installs the environment (this needs internet and takes a few minutes); a dialog explains this is normal. When it finishes, the app window opens automatically. On later launches the app opens straight away.

**If macOS blocks the file** ("cannot be opened because it is from an unidentified developer"), right-click the launcher, choose **Open**, then confirm **Open** in the dialog. You only need to do this once.

In the app, click **Browse…** for the *Image folder* and pick either the folder that contains `Plate 1` and `Plate 2`, or a single acquisition folder that contains `XY##` folders. Click **Browse…** for *Results go to* and pick a destination such as Desktop or a Results folder. Press **Run analysis**. The runner creates a new timestamped folder inside your chosen results folder, named like `Cellpose_DAPI_aSMA_YYYYMMDD_HHMMSS`.

If a run fails, the progress log in the app window contains the error details; send that text plus the run output folder if one was created.

Full plates can take a while because Cellpose loads the model and segments every field. A one-image smoke test takes roughly a minute on this Mac; a full 227-field plate set is expected to take much longer. To try a quick test first, open **Advanced options** and set *Max images per acquisition* to `2`.

## First Run On Windows

Double-click `Run Cellpose DAPI aSMA Pipeline Windows.cmd`. On the first launch a Command Prompt window opens and installs the environment (needs internet, a few minutes); a dialog explains this is normal. When it finishes, the app window opens automatically. On later launches the app opens straight away.

The first launch installs `uv` for the current user if needed, creates `.venv\Scripts\python.exe`, installs Cellpose, downloads and checksum-verifies `cpsam_v2`, and runs a discovery check. After it, the Windows folder contains:

- `.venv\Scripts\python.exe`
- Cellpose installed in that Windows `.venv`
- `.models\cellpose\cpsam_v2`

The app itself is identical to the Mac app: pick the image folder and results folder, optionally open **Advanced options**, and press **Run analysis**. The Windows launcher does not require CellProfiler, Docker, or manual plugin setup. If a run fails, the app's progress log contains the error details; send that text plus the run output folder if one was created.

Before giving the Windows folder to a routine user, the maintainer should test:

```powershell
.\.venv\Scripts\python.exe scripts\run_user_cellpose_batch.py --input data\aSMA_DAPI_plates --output output\user_cellpose_batch\windows_dry_run_placeholder --target-channel CH2 --dapi-channel CH4 --dry-run
```

Expected discovery depends on the selected folder. For the July 2026 download, the usable top folder has four acquisition folders and 367 image pairs. Then run one small real smoke test from the double-click launcher or with:

```powershell
.\.venv\Scripts\python.exe scripts\run_user_cellpose_batch.py --input "data\aSMA_DAPI_plates\plate 1\ApYYM20AGGSMA_01" --output output\user_cellpose_batch\windows_smoke --target-channel CH2 --dapi-channel CH4 --max-images-per-acquisition 1 --skip-figures
```

## Main Outputs

Each run creates a timestamped folder named like `Cellpose_DAPI_aSMA_YYYYMMDD_HHMMSS`.

The user-facing outputs are in `final/`:

- `START_HERE_RUN_SUMMARY.html`: plain-language run summary.
- `tables/cellpose_user_friendly_endpoint_summary.csv`: compact table with the primary endpoint.
- `tables/cellpose_full_plate_endpoint_summary.csv`: full audit table.
- `workbooks/cellpose_background_corrected_pi_style_summary.xlsx`: simple Excel workbook for the run. The filename is retained for compatibility; current values are raw retained target-channel/aSMA integrated intensity because `background_value_per_px = 0.0`. The workbook includes a `Channel Mapping` sheet that records which channel was used for DAPI and which was used for the target.
- `figures/`: visual QC overlays and endpoint plots, when figure rendering is enabled.

Intermediate masks, per-image QC panels, and logs are also retained:

- `cellpose_counts/`: DAPI nuclei masks, count tables, and nuclei QC.
- `cellpose_cell_regions/`: Cellpose target+DAPI region masks, metrics, and region QC.
- `logs/user_cellpose_batch_run.yaml`: machine-readable run record.

## Scientific Endpoint

The user-friendly CSV exports:

```text
target_integrated_intensity_per_DAPI_positive_nucleus
```

This is computed from the DAPI-anchored Cellpose region target-channel/aSMA integrated intensity divided by the DAPI-positive nucleus count. It is not DAPI fluorescence normalization. DAPI brightness is not used as the denominator.

More explicitly, the denominator is `filtered_DAPI_positive_nucleus_count`; current nucleus filtering is none, so this equals the Cellpose CH4/DAPI label count.

The current Cellpose region policy is DAPI-anchored: Cellpose objects are retained only if at least one DAPI-positive nucleus centroid falls inside the object.

## Technical Maintainer Notes

The setup launchers automate the commands below. They are retained here for auditability and for maintainers preparing a release zip.

For a new Mac, the equivalent maintainer commands are:

```bash
uv venv --python 3.12 .venv
.venv/bin/python -m pip install -e '.[cellpose]'
mkdir -p .models/cellpose
curl -L -o .models/cellpose/cpsam_v2.download https://huggingface.co/mouseland/cellpose-sam/resolve/main/cpsam_v2
mv .models/cellpose/cpsam_v2.download .models/cellpose/cpsam_v2
```

For a Windows workstation, the equivalent maintainer commands are:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -e ".[cellpose]"
New-Item -ItemType Directory -Force .models\cellpose | Out-Null
curl.exe -L -o .models\cellpose\cpsam_v2.download https://huggingface.co/mouseland/cellpose-sam/resolve/main/cpsam_v2
Move-Item .models\cellpose\cpsam_v2.download .models\cellpose\cpsam_v2
```

For a Windows NVIDIA workstation, a maintainer can install the machine-appropriate CUDA-enabled PyTorch build before installing this project. The default setup remains operational without CellProfiler or Docker, but GPU acceleration depends on the local PyTorch/driver combination.

## Why Not CellProfiler For This Production Runner

CellProfiler is a real application and can batch images. It is useful as a reference or expert QC tool, but it is not the recommended low-burden production runner here.

| Tool | What it offers | Why it is not the routine path here |
| --- | --- | --- |
| CellProfiler | GUI pipeline editing, batch processing, masks/overlays, object tables, `.cppipe` reproducibility | Reproducing this endpoint requires staged inputs, GUI module settings, plugin/Docker or CellProfiler-Python setup for RunCellpose, and separate project postprocessing |
| Dedicated launcher | One folder picker workflow, project-local Cellpose, direct CSV/workbook/QC export | Less flexible than CellProfiler's GUI, but much lower burden for routine runs |

The dedicated launcher uses the existing project Cellpose pipeline directly. That keeps the user action to folder picking while preserving the repository's output contract, logs, masks, QC overlays, and primary endpoint.

## Command-Line Equivalent

The launcher runs the same command shape as:

```bash
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_user_cellpose_batch.py \
  --input data/aSMA_DAPI_plates \
  --output output/user_cellpose_batch/example_run \
  --target-channel CH2 \
  --dapi-channel CH4
```

For a quick smoke test:

```bash
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_user_cellpose_batch.py \
  --input "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01" \
  --output output/user_cellpose_batch/smoke \
  --target-channel CH2 \
  --dapi-channel CH4 \
  --max-images-per-acquisition 1 \
  --skip-figures
```
