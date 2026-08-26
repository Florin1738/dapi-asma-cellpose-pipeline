# Parallel per-acquisition Cellpose runner

## Purpose and scope

`scripts/run_user_cellpose_parallel_batch.py` schedules the existing validated
Cellpose workflow separately for every discovered acquisition folder. It does
not change segmentation, measurement, normalization, workbook contents, or QC.
It only allows a bounded number of independent acquisition processes to run at
the same time.

This is the supported route when each acquisition must retain its own workbook
and QC bundle. Do not point the ordinary `run_user_cellpose_batch.py` at the
parent folder when a single workbook per acquisition is required: one ordinary
CLI invocation produces one combined run-level workbook.

## Output contract

Choose a new output root. For each acquisition, the runner creates:

```text
<output-root>/
  START_HERE_PARALLEL_RUN.md
  _orchestration/
    parallel_run_manifest.json
    logs/
      001__Plate_1__<acquisition>.log
  Plate_1/
    <acquisition>/
      final/
        workbooks/cellpose_background_corrected_pi_style_summary.xlsx
        tables/
        figures/
      cellpose_counts/
      cellpose_cell_regions/
      logs/
```

The workbook, QC figures, masks, method parameters, and child log remain
isolated under that acquisition's output directory. The orchestration manifest
records the exact child command, source folder, channel mapping, status,
duration, expected workbook path, and log location for every acquisition.

The runner refuses a non-empty output root so a benchmark or production run
cannot silently mix with previous results.

## Preflight and controlled benchmark

Always validate discovery and paths before writing any data:

```powershell
& .\.venv\Scripts\python.exe .\scripts\run_user_cellpose_parallel_batch.py `
  --input "<parent folder containing acquisitions>" `
  --output "<new output root>" `
  --target-channel CH2 `
  --dapi-channel CH4 `
  --gpu `
  --workers 2 `
  --dry-run
```

For a controlled concurrency test on a single GPU, use a fresh benchmark output
root and a small per-acquisition image cap. This tests two independently loaded
Cellpose processes without conflating the result with the production outputs:

```powershell
& .\.venv\Scripts\python.exe .\scripts\run_user_cellpose_parallel_batch.py `
  --input "<parent folder containing two representative acquisitions>" `
  --output "<new benchmark output root>" `
  --target-channel CH2 `
  --dapi-channel CH4 `
  --gpu `
  --workers 2 `
  --max-images-per-acquisition 12
```

Review both child logs, the two visual-QC bundles, and the manifest. Compare
wall-clock time against a one-worker benchmark over the same image count. Use
two workers for production only if both runs complete on CUDA, visual QC is
acceptable, and throughput improves. Keep `--workers 1` if GPU-memory pressure,
errors, or negligible speed-up are observed.

Add `--acquisition "<exact folder name>"` once for each folder in a small,
restartable batch. The option selects only exact acquisition-folder or relative
display names and rejects an ambiguous or unknown name. For example:

```powershell
  --acquisition "003-7b TW aSMA" `
  --acquisition "0037a a-SMA TW"
```

## Production command

After the benchmark is accepted, use a different, empty production output root
and omit `--max-images-per-acquisition`:

```powershell
& .\.venv\Scripts\python.exe .\scripts\run_user_cellpose_parallel_batch.py `
  --input "<parent folder containing acquisitions>" `
  --output "<new production output root>" `
  --target-channel CH2 `
  --dapi-channel CH4 `
  --gpu `
  --workers 2
```

The default is fail-safe: a failed child stops scheduling folders that have not
started, while already-running isolated folders are allowed to finish. Add
`--continue-on-error` only when the operator explicitly wants independent later
acquisitions to proceed after a failure.
