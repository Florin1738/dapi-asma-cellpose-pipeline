# Manual Validation Folder

This folder contains packages intended for human annotation and validation of candidate segmentation masks.

## Current Package

- `plate1_selected15_propagation_otsu_raw0_asma_region/`: current selected15 manual-validation package. It contains manifests, guide panels, blank reference-mask targets, raw annotation exports, audit outputs, and a validation gate report.

## Superseded Package

- `plate1_xy22_xy23_xy24_xy40_xy41_propagation_otsu_raw0_asma_region/`: earlier 5-field package retained for provenance.

## Important Status

The current manual-validation package is not complete until reference masks are filled by a human reviewer. Placeholder all-zero masks must not be treated as ground truth.

Quantitative validation metrics such as precision, recall, F1, false positives, false negatives, and IoU should only be reported after non-empty manual reference masks exist and the validation gate passes.
