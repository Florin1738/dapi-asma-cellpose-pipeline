# Scripts Folder

This folder contains command-line entry points for reproducible analysis tasks.

## Common Entry Points

- `run_pi_simple_summary.py`: generate the PI-requested workbook and side plots.
- `run_cellpose_counts.py`: count DAPI-positive nuclei with Cellpose.
- `run_cellpose_cell_regions.py`: run the selected/provisional Cellpose CH2+CH4 candidate aSMA-positive region workflow.
- `render_method_review_panels.py`: render qualitative method-comparison overlays.
- `prepare_manual_validation_package.py`: create manual-validation packages.
- `run_manual_validation_report.py`: summarize manual validation after reference masks are available.

## Exploratory Scripts

Several scripts support background searches, histogram panels, threshold visualizations, seeded-region comparators, and sensitivity analysis. Check `docs/ARTIFACT_LEDGER.md` and `docs/methodology/run_logs/` before treating an exploratory script as the current method.
