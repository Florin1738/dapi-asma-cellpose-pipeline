# Project Contents

This file is a quick map for the GitHub repository and release zips. It separates distributable source files from local data and generated analysis artifacts.

## Start Here

- `README.md`: project overview, quick-start instructions, endpoint definition, and developer setup.
- `READ ME FIRST.txt`: plain-language instructions bundled into release zips.
- `docs/INSTALLATION.md`: Mac, Windows, and developer setup details.
- `docs/NONTECHNICAL_CELLPOSE_BATCH_RUNNER.md`: routine double-click app workflow and expected outputs.
- `docs/RELEASE_PACKAGING.md`: maintainer steps for building and publishing release zips.
- `AGENTS.md`: scientific and engineering rules for future coding agents.

## User Launchers

- `Run Cellpose DAPI aSMA Pipeline.command`: Mac double-click launcher.
- `Run Cellpose DAPI aSMA Pipeline Windows.cmd`: Windows double-click launcher.
- `Setup Cellpose DAPI aSMA Pipeline.command`: optional Mac setup/repair launcher.
- `Setup Cellpose DAPI aSMA Pipeline Windows.cmd`: optional Windows setup/repair launcher.

The run launchers self-heal: if the local environment is missing, they run setup before opening the app.

## Active Code

- `src/dapi_norm/`: Python package for image I/O, Cellpose wrappers, measurement, QC rendering, validation, and report-table generation.
- `scripts/`: reproducible command-line entry points plus Mac/Windows setup helpers.
- `configs/`: example configuration files.
- `examples/`: small command examples.
- `tests/`: synthetic and artifact-level tests that should not require private microscopy data.

## Documentation

- `docs/OUTPUT_CONTRACT.md`: stable output columns and artifact expectations.
- `docs/IMAGE_INPUTS.md`: supported inputs, channel policy, and axis-handling rules.
- `docs/QC_AND_VISUALIZATION.md`: visual QC expectations.
- `docs/VALIDATION.md`: manual centroid and mask validation protocols.
- `docs/methodology/`: method status, thresholds, and run-log templates.

## Local-Only Folders

These folders may exist in a working checkout but their contents are intentionally ignored by git:

- `data/`: private microscopy inputs.
- `output/`: generated analysis outputs.
- `manual_validation/`: manual annotation and validation packages.
- `reports/`: generated report packages, workbooks, figures, and renders.
- `.models/`: Cellpose model cache.
- `.venv/`: project-local Python environment.
- `dist/`: release zip outputs.

Only folder-level README files are tracked for these local-only directories.

## Git Tracking Policy

Track source code, tests, configs, small README/index files, setup launchers, and public-facing documentation.

Do not track real microscopy images, masks, generated figures, generated tables, Excel workbooks, report PDFs/DOCX files, Cellpose model weights, virtual environments, or lab-private annotations.
