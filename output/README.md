# Generated Output Folder

This folder contains generated analysis outputs. Many files are ignored by git because they include microscopy-derived masks, figures, tables, and workbooks.

## Active / User-Facing

- `pi_simple_summary/`: PI-requested workbook, plots, DAPI count inputs, and visual QC panels.
- `cellpose_cell_regions/full_plate_cpsam_v2/`: current full-plate Cellpose CH2+CH4 candidate aSMA-associated region run used by the active report and workbook.
- `manual_validation_selection/plate1_ApYYM20AGGSMA_01_stratified_16/`: selection rationale for older manual-validation fields.

## Support / Method Review

- `method_review/`: qualitative comparison panels across segmentation approaches.
- `cellpose_cell_regions/plate1_ApYYM20AGGSMA_01_selected15_cpsam_v2/`: superseded selected-field Cellpose run retained for method-development provenance.
- `cellpose_counts/`: earlier sample count outputs retained for compatibility with existing commands and docs.
- `target_normalization/`: earlier sample target per-nucleus output with legacy directory naming.

## Exploratory / Set Aside

- `seeded_asma_regions/`: watershed, random-walker, and propagation-based candidate-region tests.
- `background_candidate_search/` and `background_exclusion_exploration/`: background-estimation experiments.
- `ch2_absorbance/`, `ch2_raw_intensity/`, `ch2_threshold_retained_area/`, `od_spherical/`: histogram and threshold diagnostics.
- `method_triage/` and `sensitivity/`: method-selection summaries.

See `docs/ARTIFACT_LEDGER.md` before deleting or moving any output folder.
