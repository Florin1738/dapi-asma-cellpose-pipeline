# Source Package Folder

The Python package lives in `src/dapi_norm/`.

## Module Responsibilities

- Image input and channel extraction: `image_arrays.py`, `data_inventory.py`, `data_cli.py`
- Cellpose runtime/counting: `cellpose_runner.py`, `nuclei_outputs.py`
- PI workbook: `pi_simple_summary.py`
- Candidate aSMA region workflows: `cellpose_cell_regions.py`, `seeded_regions.py`, `method_triage.py`
- QC rendering: `qc_preview.py`, `centroid_overlays.py`, `method_review_panels.py`
- Manual validation: `manual_validation_package.py`, `segmentation_validation.py`, `manual_validation_report.py`, annotation helper modules
- Target per-nucleus measurement: `target_normalization.py`, `target_validation.py`

Keep modules small and named by responsibility. Use the output contract in `docs/OUTPUT_CONTRACT.md` for stable column names.
