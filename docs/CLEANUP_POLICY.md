# Cleanup Policy

This repository contains source code, private microscopy inputs, generated outputs, exploratory diagnostics, and report artifacts. Cleanup should preserve reproducibility before reducing disk clutter.

## Do Not Delete Without Explicit Approval

- `data/aSMA_DAPI_plates/`
- `ApYYM20AGGSMA_02/`
- `.models/cellpose/`
- `output/pi_simple_summary/`
- `output/cellpose_cell_regions/full_plate_cpsam_v2/`
- `output/cellpose_cell_regions/plate1_ApYYM20AGGSMA_01_selected15_cpsam_v2/`
- `manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/`
- `reports/cellpose_ch2_ch4_full_plate/`
- `reports/cellpose_ch2_ch4_asma_quantification/`

These paths are either private source inputs, current user-facing deliverables,
dependencies of the current full-plate report, or retained provenance for older
method-development work.

## Regenerable But Keep Until Reviewed

- `output/method_review/`
- `output/manual_validation_selection/`
- `output/pi_simple_summary/qc/`
- `output/pi_simple_summary/plots/`

These are reproducible support artifacts, but they are useful for visual review and method decisions.

## Exploratory / Set Aside

The following categories can be archived or removed later after confirming they are not referenced by current docs, reports, or tests:

- `output/background_candidate_search/`
- `output/background_exclusion_exploration/`
- `output/ch2_absorbance/`
- `output/ch2_threshold_retained_area/`
- `output/od_spherical/`
- non-selected runs under `output/seeded_asma_regions/`
- non-selected runs under `output/cellpose_cell_regions/`
- `output/method_triage/`
- `output/sensitivity/`

## Rules for Future Cleanup

1. Update `docs/ARTIFACT_LEDGER.md` before moving or deleting a major output.
2. Check references with `rg "<folder_or_filename>"`.
3. Preserve the command or script needed to regenerate an artifact before deleting it.
4. Do not commit real microscopy images, masks, generated outputs, model weights, or lab-private annotations.
5. If a folder contains manual annotations, treat it as non-regenerable unless the annotation source is separately preserved.
