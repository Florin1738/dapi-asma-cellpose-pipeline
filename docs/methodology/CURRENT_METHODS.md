# Current Methods

Last updated: 2026-07-06.

## Biological Question

The assay is quantifying immunofluorescence images of alpha-smooth muscle actin
with DAPI nuclei used as the cell-count reference. The default/original mapping
is `CH2` = aSMA target and `CH4` = DAPI. The nontechnical Cellpose runner now
accepts explicit target and DAPI channel selections because new exports may keep
only `XY##` folder naming consistent while channel assignments differ. In this
experimental context, aSMA is a fibrosis/myofibroblast-activation readout.

The computational goal is not to interpret biology directly. The repository
produces auditable image-derived metrics that can later be grouped by treatment,
dose, plate, and replicate.

## Current PI Workbook Metric

Status: accepted as the PI-requested simple deliverable.

Formula:

```text
aSMA intensity = sum(all CH2 pixels)
Nuclei Count = Cellpose count from CH4/DAPI
Ratio = aSMA intensity / Nuclei Count
```

What it measures:

- whole-field raw aSMA burden
- divided by DAPI-positive nucleus count

What it does not measure:

- DAPI fluorescence normalization
- true per-cell aSMA expression
- cytoplasmic aSMA inside validated cell boundaries
- stress-fiber organization
- background-corrected intensity

Known limitation:

Dense fields can appear falsely low after division by nuclei count. Plate 1
`XY22` demonstrated this: visually high aSMA burden was flattened by a very high
nuclei denominator.

## Current DAPI/Nuclei Count Method

Status: accepted as the current automated count method, pending manual validation.

Method:

- Input channel: `CH4`, user-confirmed DAPI/nuclei.
- Backend: Cellpose.
- Model: `cpsam_v2`.
- Output: label masks, nucleus counts, per-nucleus centroids, and QC overlays.

Parameters are documented in `THRESHOLDS_AND_PARAMETERS.md` and in each run's
`config_resolved.yaml`.

## Background-Corrected Whole-Field aSMA

Status: exploratory/provisional.

Formula:

```text
background_value = percentile_10(CH2 image)
corrected_pixel = max(CH2_pixel - background_value, 0)
corrected_integrated = sum(corrected_pixel across full field)
corrected_per_nucleus = corrected_integrated / DAPI nucleus count
```

Rationale:

This is a simple image-level background estimate when no manual cell-free
background ROI exists. It is not a true cell mask and not local background
correction.

Do not present this as final cell-level expression.

## CellProfiler ECM CH1 Quantification

Status: current reviewable ECM output for the July 2026 download; not manual
ground-truth validated.

Formula:

```text
background_value = median(CH1 pixels outside the CellProfiler ECM-positive mask)
ecm_positive_integrated_background_corrected =
    sum(max(CH1_pixel - background_value, 0) inside the ECM-positive mask)
```

Current run:

```text
reports/cellprofiler_ecm_ch1_2026_07_06_k2_arm64/
```

Current CellProfiler threshold:

```text
module = IdentifyPrimaryObjects
method = Robust Background
strategy = Global
deviations = 2
lower/upper outlier fractions = 0.05 / 0.05
threshold correction factor = 1.0
```

DAPI role:

- `CH4` DAPI is retained as context/QC metadata.
- ECM is not normalized by DAPI-positive nucleus count.
- DAPI fluorescence intensity is not used for ECM normalization.

Selection status:

The earlier `k=3` output was selected by an automatic/qualitative heuristic and
was later judged too aggressive by human QC. The current `k=2` setting is a
reviewable threshold choice supported by stored representative threshold-sweep
panels and aggregate diagnostics, not a ground-truth-validated optimum.

Do not present this as ECM-positive cell counts or true cell segmentation.

## Cell-Restricted aSMA Status

Status: current full-plate Cellpose candidate-region analysis is runnable and
reported; validated whole-cell segmentation remains unresolved.

The current images do not include a separate pan-cell marker. Therefore, true
whole-cell segmentation is not yet validated.

Current active full-plate package:

```text
output/cellpose_cell_regions/full_plate_cpsam_v2/
reports/cellpose_ch2_ch4_full_plate/
```

Current nontechnical batch runner:

```text
Run Cellpose DAPI aSMA Pipeline.command
Run Cellpose DAPI aSMA Pipeline Windows.cmd
scripts/run_user_cellpose_batch.py
docs/NONTECHNICAL_CELLPOSE_BATCH_RUNNER.md
```

Status: accepted production wrapper around the current Cellpose-first method.
It does not change the segmentation or measurement definition; it automates
folder discovery, explicit channel mapping, per-acquisition count/region runs,
merged CSV export, PI-style workbook export, run logs, and QC overlays on Mac
and Windows prepared project folders.

Current active user-facing Cellpose measurements:

- `whole_field_target_integrated_raw`: sum of every selected target/aSMA-channel
  pixel in the image. The active target channel is recorded as
  `target_channel_id`.
- `cellpose_target_integrated_raw`: sum of selected target/aSMA-channel pixels
  inside the Cellpose retained-region mask. In current user-facing outputs,
  Cellpose retained regions are Cellpose objects that contain at least one
  DAPI-positive nucleus centroid.
- Each corresponding per-nucleus endpoint divides by filtered DAPI-positive
  nucleus count. Current nucleus filtering is none, so this equals the Cellpose
  DAPI-channel label count. The active DAPI channel is recorded as
  `dapi_channel_id`. DAPI intensity is not used for normalization.

Interpretation of masking:

```text
background_value_per_px = 0.0
```

The current Cellpose workbook and report do not subtract a scalar background
value. They exclude pixels outside the retained-region mask. Objects with no
DAPI-positive nucleus centroid inside the object are excluded before Cellpose
intensity and area are reported.

Scientific status:

- The output is appropriate for auditing how full-field raw signal changes when
  restricted to DAPI-anchored Cellpose candidate aSMA-associated regions.
- The retained-region rule directly addresses the concern that a separated
  Cellpose-positive object can be counted as positive area even when no DAPI
  nucleus is detected in that object.
- Automated QC status is carried in the current summary table and workbook.
  Current full-plate counts are `reviewable_not_validated: 20`,
  `needs_manual_review: 163`, and `reject_qc_failure: 44`. These are automated
  review flags, not manual validation.
- The method should still be described as candidate-region analysis, not as
  validated whole-cell segmentation, because `CH2`/aSMA is the biological
  endpoint signal rather than an independent pan-cell marker.

Methods evaluated:

- DAPI-derived spatial expansion territories: exploratory, not robust enough
  for final use.
- CH2/aSMA threshold masks: exploratory only; they define foreground using the
  signal being measured.
- DAPI-seeded CH2/aSMA foreground partitioning by marker-controlled watershed:
  implemented exploratory comparator; not validated whole-cell segmentation.
- DAPI-seeded CH2/aSMA foreground partitioning by random walker: implemented
  exploratory comparator. It changed some nucleus-to-pixel assignment inside the
  retained foreground, but did not change the retained foreground union or the
  current union-level intensity endpoint.
- DAPI-seeded CH2/aSMA foreground partitioning by CellProfiler-style
  propagation: implemented exploratory comparator using `centrosome.propagate`.
  It is the closest project-local analogue to CellProfiler
  `IdentifySecondaryObjects` propagation. It preserves the intuitive
  `XY22 > XY23 > XY24` ordering, rejects low-DAPI challenge fields by QC, and
  remains exploratory because the CH2/aSMA endpoint signal defines the
  foreground mask.
- Cellpose CH2+CH4 candidate object segmentation: implemented exploratory
  comparator using a pretrained Cellpose `cpsam_v2` model with CH2/aSMA as the
  candidate object/cytoplasm signal and CH4/DAPI as nuclear context. In the
  first Plate 1 comparison, it produced `XY22 > XY23 > XY24` candidate-region
  intensity per DAPI-positive nucleus. Strict QC rejected `XY23` and `XY24`
  because substantial background-corrected CH2 remained outside the Cellpose
  candidate regions, and it rejected low DAPI-count challenge fields `XY40` and
  `XY41`. It is still not accepted as true whole-cell segmentation because
  CH2/aSMA is the biological endpoint, not a pan-cell stain.
- CellProfiler secondary-object/propagation workflows: candidate reference
  workflow, but they need careful testing because robust use typically requires
  a cell-body, membrane, cytoskeleton, or brightfield-derived boundary image.

Cross-method triage:

- The current triage report compares the PI whole-field metric, DAPI-seeded
  Otsu/watershed regions, DAPI-seeded Otsu/random-walker regions, DAPI-seeded
  CellProfiler-style propagation regions, and Cellpose CH2+CH4 candidate
  regions for `XY22`, `XY23`, `XY24`, `XY40`, and `XY41`.
- Output:
  `output/method_triage/plate1_xy22_xy23_xy24_xy40_xy41/`.
- Conclusion: no region-restricted method is accepted as validated whole-cell
  segmentation. `XY22` remains reviewable but not validated; `XY23` and `XY24`
  are explicitly marked as mixed QC-rejection rows; `XY40` and `XY41` are
  marked as all available region-restricted methods rejected.
- This triage is a method-selection aid only. It does not replace the accepted
  PI workbook metric and does not support biological interpretation.

Side-by-side visual method review:

- Output:
  `output/method_review/plate1_xy22_xy23_xy24_xy40_xy41_propagation_vs_cellpose/`.
- Expanded selected-field output:
  `output/method_review/plate1_selected15_propagation_vs_cellpose/`.
- The package contains matched full-field and crop panels comparing CH4/DAPI
  nuclei, raw CH2/aSMA, CellProfiler-style propagation/Otsu regions, Cellpose
  CH2+CH4 candidate regions, and a method-disagreement map.
- Result: Cellpose retained much larger candidate regions than propagation in
  `XY22`, `XY23`, and `XY24`; the method-region Jaccard values were
  approximately `0.423`, `0.230`, and `0.365`. Both methods remained
  QC-rejected in low-nucleus challenge fields `XY40` and `XY41`.
- On the stratified 15-field manual-validation candidate set, Cellpose was
  `reviewable_not_validated` only for `XY22`; most dense positive fields were
  `needs_manual_review`, and low/low-signal or artifact-risk fields were
  rejected by QC. Propagation/Otsu remained more conservative and was
  reviewable for most non-challenge fields, but is still not validated.
- Interpretation: the panels make the tradeoff visually auditable. Propagation
  is conservative and preserves the expected `XY22 > XY23 > XY24` ordering.
  Cellpose captures broader cell-like regions but retains substantial
  method-only area and is rejected by strict QC in `XY23`, `XY24`, `XY40`, and
  `XY41`. This is qualitative QC only, not validation.

Raw0 threshold/method sensitivity:

- Output:
  `output/sensitivity/plate1_xy22_xy23_xy24_xy40_xy41_region_restricted_raw0/`.
- Review crop panel:
  `output/seeded_asma_regions/plate1_ApYYM20AGGSMA_01_xy22_xy23_xy24_xy40_xy41_propagation_otsu_reg005_raw0/qc_review_crops_green_x_segmentation.png`.
- All sensitivity runs used raw intensity inside retained regions
  (`background_value_per_px = 0`).
- With the same Otsu CH2 foreground mask, watershed, random walker, and
  CellProfiler-style propagation give identical image-level retained-region
  totals. This is expected because the current endpoint sums the retained
  foreground union; the seed-assignment algorithm only matters for future
  per-cell/per-object attribution.
- Per-image Otsu is the only tested automatic threshold family that preserved
  the expected `XY22 > XY23 > XY24` ordering and returned zero retained region
  for low-DAPI challenge fields `XY40` and `XY41`.
- Li preserved `XY22 > XY23 > XY24`, but produced very large rejected challenge
  values if QC is ignored. Triangle failed the expected ordering and also
  produced nonzero challenge values.
- Decision: use Otsu as the current exploratory region-restricted threshold for
  review panels and robustness checks. This is not a promotion to validated
  whole-cell segmentation.

## Current ECM CH1 Quantification

Status: current reviewable CellProfiler-backed ECM output for the July 2026
download; not a manual segmentation validation.

Current output:

```text
reports/cellprofiler_ecm_ch1_2026_07_06_k2_arm64/
```

Input channel mapping:

```text
CH1 = ECM
CH4 = DAPI nuclei
```

Endpoint:

```text
ecm_positive_integrated_background_corrected
  = sum(max(CH1 - background, 0) inside the CellProfiler ECM-positive mask)
```

Background is the per-image median of CH1 pixels outside the CellProfiler
ECM-positive mask. DAPI-positive nucleus count is retained as context/QC
metadata only. DAPI count and DAPI fluorescence brightness are not used to
normalize ECM.

Mask method:

```text
Actual CellProfiler IdentifyPrimaryObjects
Threshold method = Robust Background
# of deviations = 2
selection_source = human_qc_override
```

The earlier `k=3` ECM run is superseded because human QC judged it too
aggressive. The current `k=2` run retained more weaker ECM signal while avoiding
near-full-field masks in the representative sweep. No manual ECM ground truth
masks are available, so do not report precision, recall, F1, false positives,
false negatives, or IoU.

## Preferred Near-Term Reporting

For current review, report these separately:

- total raw aSMA burden per image/well
- DAPI nuclei count
- raw aSMA burden per DAPI-positive nucleus, with caveat
- Cellpose retained-region aSMA burden
- Cellpose retained-region aSMA burden per DAPI-positive nucleus
- saturation fraction
- visual QC overlays
- exploratory background-corrected whole-field metrics, clearly labeled if used

Do not collapse all of these into one biological conclusion.
