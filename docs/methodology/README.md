# Methodology Documentation

This folder is the project's methodology memory. It records what is being measured,
why each method exists, which thresholds or parameters were chosen, which sources
support the choices, and which methods are accepted versus exploratory.

## Folder Map

```text
docs/methodology/
  README.md
  CURRENT_METHODS.md
  THRESHOLDS_AND_PARAMETERS.md
  CELL_SEGMENTATION_STRATEGY.md
  MANUAL_MASK_VALIDATION.md
  ADVANCED_ALPHA_SMA_FIBER_METRICS.md
  RUN_LOG_TEMPLATE.md
  run_logs/
    README.md
    YYYY-MM-DD-short-topic.md
```

## Required Logging Rule

Every new analysis method, segmentation rule, threshold, parameter change, or
source-backed methodological claim must get a run log in:

```text
docs/methodology/run_logs/
```

The run log must include:

- status: accepted, exploratory, rejected, or superseded
- input data and output paths
- exact command or script
- exact parameters and thresholds
- why each threshold or parameter was chosen
- source links used to justify the method
- visual QC outputs
- validation status and caveats

## Current Method Status

- Accepted PI-facing deliverable: raw CH2/aSMA integrated intensity, CH4/DAPI
  nuclei count, and raw intensity divided by nuclei count.
- Accepted count method: Cellpose segmentation of CH4/DAPI nuclei, pending
  manual validation.
- Active full-plate candidate-region package: Cellpose CH2+CH4 aSMA-associated
  regions under `output/cellpose_cell_regions/full_plate_cpsam_v2/`, reported
  in `reports/cellpose_ch2_ch4_full_plate/`. In current user-facing outputs,
  Cellpose means retained Cellpose objects with at least one DAPI-positive
  nucleus centroid inside the object.
- Exploratory methods: percentile background correction, alpha-SMA foreground
  thresholding, DAPI-derived cell-neighborhood territories, DAPI-seeded
  aSMA-associated regions, and Cellpose CH2+CH4 candidate aSMA-associated
  regions.
- Exploratory method-selection report: cross-method triage in
  `output/method_triage/plate1_xy22_xy23_xy24_xy40_xy41/`; no
  region-restricted method is currently accepted as validated whole-cell
  segmentation.
- Rejected as final cell segmentation: the current DAPI expansion territory
  prototype. It is useful diagnostically, but it does not correctly capture
  irregular cell bodies in the available images.
- Validation infrastructure ready, awaiting manual labels: IoU-based manual
  mask validation, manual-validation package generation, annotation-readiness
  audit, raw-only editor exports, annotation handoff bundles, safe
  manual-reference-mask commit tooling, and candidate-method validation report
  generation described in `MANUAL_MASK_VALIDATION.md`.

Current manual-validation package:

```text
manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/
```

The current package contains guide panels, a manifest, blank reference masks,
raw-only annotation panels, raw-only editor TIFF exports, annotation handoff
bundles, a bias-aware annotation review gallery, and an annotation audit for a
stratified 15-field set:

```text
XY22, XY23, XY24, XY40, XY41, XY01, XY74, XY66, XY08, XY95,
XY13, XY10, XY11, XY33, XY16
```

The set covers explicit challenge fields, low/mid/high raw aSMA, low/high
DAPI-positive nucleus count, high per-nucleus artifacts, high CH2 saturation,
and method-disagreement fields. The status file starts each image as
`not_started` and must be changed to `complete_non_empty` or `confirmed_empty`
before validation. The reference masks are intentionally empty until a human
fills them; quantitative precision/recall/F1/IoU must not be reported from this
package yet.

Earlier five-field package retained for traceability:

`manual_validation/plate1_xy22_xy23_xy24_xy40_xy41_propagation_otsu_raw0_asma_region/`

Current annotation handoff package:

```text
manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/annotation_handoff/
```

It contains per-image NPZ layer bundles for napari or another label editor. The
candidate and nuclei labels are context layers only; the editable layer is
`manual_reference_labels`, and completed labels must be committed to the
authoritative TIFFs in `reference_masks_to_fill/` with
`scripts/commit_manual_reference_mask.py`.

Current raw-only editor export:

```text
manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/raw_annotation_exports/
```

It contains one folder per selected field with only raw CH2/aSMA, raw CH4/DAPI,
and an editable manual-reference-label scratch TIFF. Candidate masks are
intentionally absent so first-pass annotation can be done without tracing the
candidate segmentation. After editing a scratch label TIFF, commit it with the
per-image command template in `raw_annotation_export_manifest.csv` after
replacing `YOUR_INITIALS`.

Current raw-annotation import summary:

```text
manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/raw_annotation_import/bulk_import_summary.csv
```

It currently reports all 15 fields as `skipped_empty_unconfirmed`, meaning the
scratch TIFFs are still blank and no authoritative manual/reference masks were
changed.

Current annotation review gallery:

```text
manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/annotation_review_gallery/index.html
```

It shows each selected field's raw-only annotation panel, status, audit blocking
reasons, and exact safe commit commands. Candidate guide panels stay hidden
until a field is marked complete because they are post-label QC only.

Current annotation audit:

```text
manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/annotation_audit/
```

It reports 0 of 15 images validation-ready because all status rows remain
`not_started`.

Current validation gate:

```text
manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/validation_gate/
```

Run `scripts/run_manual_validation_pipeline.py` after manual/reference labels
are committed. It reruns the annotation audit first, refuses to score incomplete
packages, and only then runs quantitative candidate-vs-reference validation.
The current gate report is blocked for all 15 fields because labels are still
`not_started`.

Current Cellpose selected-15 visual review report:

```text
output/cellpose_cell_regions/plate1_ApYYM20AGGSMA_01_selected15_cpsam_v2/review_report/index.html
```

This report groups fields by Cellpose candidate-region QC status and links each
full QC panel plus excluded-signal check. It is qualitative review only; it does
not replace the manual/reference-mask validation gate.

## Key Constraint

The project currently has DAPI and alpha-SMA images only. There is no separate
general cytoplasm, membrane, phalloidin, or whole-cell stain. That limits the
ability to produce true cell segmentation without either a geometric
approximation or a biased target-channel-driven segmentation.
