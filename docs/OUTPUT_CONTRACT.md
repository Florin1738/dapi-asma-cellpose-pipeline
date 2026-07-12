# Output Contract

All outputs must be auditable. Every CSV row should be traceable to an input image, resolved config, segmentation mask, and QC overlay.

## Directory Layout

```text
output/
  logs/
    run_log.txt
    config_resolved.yaml
  summaries/
    image_level_summary.csv
    per_nucleus_locations.csv
    roi_level_summary.csv
    validation_metrics.csv
  masks/
    <image_id>_nuclei_labels.tif
  qc/
    <image_id>_dapi_nucleus_outlines.png
    <image_id>_target_with_nucleus_outlines.png
    <image_id>_dapi_numbered_centroids.png
    <image_id>_segmentation_montage.png
  validation/
    <image_id>_validation_overlay.png
  debug/
    optional intermediate images
```

## Current Count-Only Output

The current implemented Cellpose command is a count-only pilot, not the full target-intensity endpoint. It writes:

```text
output/cellpose_counts/
  logs/
    config_resolved.yaml
    run_log.txt
  summaries/
    nucleus_counts.csv
    per_nucleus_locations.csv
  masks/
    <position_id>_<channel_id>_<model>_labels.tif
  qc/
    <position_id>_<channel_id>_<model>_montage.png
  qc_contact_sheet.png
```

`nucleus_counts.csv` is intentionally smaller than the planned full `image_level_summary.csv`. It records input image, backend, model, candidate channel, unconfirmed channel warning, mask path, montage path, and nucleus count.

`config_resolved.yaml` records software versions, requested device, model path/hash when available, channel extraction policy, segmentation parameters, filtering state, output paths, and one record per processed image. `run_log.txt` is a short human-readable summary of the same run.

## Current Target Per-Nucleus Output

The current target per-nucleus command is a provisional full-field measurement. It writes to a legacy-named directory:

```text
output/target_normalization/
  logs/
    config_resolved.yaml
    run_log.txt
  summaries/
    image_level_summary.csv
    well_level_summary.csv
  plots/
    normalized_intensity_by_well.png
    target_integrated_vs_nucleus_count.png
  qc/
    <position_id>_CH2_target_with_CH4_nucleus_outlines.png
  qc_contact_sheet.png
```

`image_level_summary.csv` uses the full output contract columns plus `well_id`, `mask_path`, and `qc_overlay_path` for traceability. `well_level_summary.csv` is currently identical to the image-level rows because each `XY##` position has one image pair in the local sample data.

Legacy output paths and plot names may include `target_normalization` or
`normalized_intensity`. In this project, those names mean target CH2/aSMA
intensity divided by DAPI-positive nucleus count. They must not be interpreted
as normalization by DAPI fluorescence brightness.

## Current Full-Plate Cellpose Channel-Selectable Output

The current full-plate Cellpose retained-region workflow writes per-acquisition
outputs under:

```text
output/cellpose_cell_regions/full_plate_cpsam_v2/
  Plate_1/
    ApYYM20AGGSMA_01/
    ApYYM20AGGSMA_02/
  Plate_2/
    APIPIKEALDSMA/
    APIPIKEALDSMA_01/
```

Each acquisition folder contains:

```text
logs/config_resolved.yaml
logs/run_log.txt
summaries/cellpose_cell_region_image_metrics.csv
masks/<location>_cellpose_<target-channel>_<dapi-channel>_cpsam_v2_labels.tif
qc/<location>_cellpose_cell_region_qc.png when internal QC was generated
plots/
```

Some full-plate acquisition runs were executed with internal per-field QC image
writing skipped to reduce runtime. The current report-level visual QC for every
field is therefore:

```text
reports/cellpose_ch2_ch4_full_plate/figures/cellpose_overlay_pages/
reports/cellpose_ch2_ch4_full_plate/figures/cellpose_overlay_pages/overlay_index.csv
```

The merged full-plate table is:

```text
reports/cellpose_ch2_ch4_full_plate/tables/cellpose_full_plate_endpoint_summary.csv
```

Current runs record the selected channel mapping explicitly:

```text
target_channel_id
dapi_channel_id
target_path
dapi_path
```

For backward compatibility, current tables may retain `*_ch2_*`, `ch2_path`,
and `ch4_path` aliases when the target is `CH2` and DAPI is `CH4`. New
user-facing logic should prefer the generic `target_*` and `dapi_*` fields.

Important definitions:

```text
whole_field_ch2_integrated_raw
  sum of all CH2/aSMA pixel intensities in the field.

whole_field_ch2_integrated_raw_per_DAPI_positive_nucleus
  whole_field_ch2_integrated_raw divided by DAPI-positive nucleus count.

cellpose_masked_ch2_integrated_raw
  audit column: sum of CH2/aSMA pixel intensities inside all Cellpose candidate
  regions before DAPI-anchored retention filtering.

cellpose_masked_ch2_integrated_raw_per_DAPI_positive_nucleus
  audit column: cellpose_masked_ch2_integrated_raw divided by DAPI-positive
  nucleus count.

dapi_anchored_cellpose_ch2_integrated_raw
  current user-facing Cellpose intensity: sum of CH2/aSMA pixel intensities
  inside Cellpose objects that contain at least one DAPI-positive nucleus
  centroid.

dapi_anchored_cellpose_ch2_integrated_raw_per_DAPI_positive_nucleus
  current user-facing Cellpose per-nucleus endpoint:
  dapi_anchored_cellpose_ch2_integrated_raw divided by DAPI-positive nucleus
  count.
```

In this run, `background_value_per_px = 0.0`. Therefore the
`*_background_corrected` values are numerically identical to the corresponding
raw retained-pixel sums. "Masked" or "background-removed" in the current report
means pixels outside retained Cellpose regions were excluded; it does not mean
a scalar background value was subtracted.

Key current user-facing full-plate columns:

```text
plate
source_id
location
whole_field_ch2_integrated_raw
whole_field_ch2_integrated_raw_per_DAPI_positive_nucleus
dapi_positive_nucleus_count
dapi_anchored_cellpose_ch2_integrated_raw
dapi_anchored_cellpose_ch2_integrated_raw_per_DAPI_positive_nucleus
dapi_anchored_cellpose_ch2_integrated_background_corrected
dapi_anchored_cellpose_ch2_integrated_background_corrected_per_DAPI_positive_nucleus
dapi_anchored_cellpose_masked_area_px
dapi_anchored_cellpose_masked_area_per_DAPI_positive_nucleus
cellpose_object_count
qc_status
qc_flags
source_warnings
cellpose_mask_path
dapi_nuclei_mask_path
source_qc_panel_path
source_excluded_signal_check_path
ch2_path
ch4_path
```

The machine-readable CSV also carries audit/provenance columns for the raw
candidate label set and excluded no-DAPI objects. Those columns are for
traceability and QC; they should not be used for polished user-facing Cellpose
endpoint plots or workbooks unless a future request explicitly asks for a
method-comparison analysis.

The active polished Excel workbook is:

```text
reports/cellpose_ch2_ch4_full_plate/workbooks/cellpose_background_corrected_pi_style_summary.xlsx
```

It contains one sheet per plate and presents `LOCATION`, Cellpose
retained-region `aSMA intensity`, `Nuclei Count`, and `Ratio`. The ratio is
Cellpose retained-region intensity divided by DAPI-positive nucleus count.

The broader audit workbook
`reports/cellpose_ch2_ch4_full_plate/workbooks/cellpose_full_plate_asma_dapi_summary.xlsx`
is retained locally for provenance and QA, but it is not the polished PI-style
endpoint workbook.

## `image_level_summary.csv`

Required columns:

```text
image_id
input_path
backend
dapi_channel
target_channel
raw_nucleus_count
filtered_nucleus_count
target_area_px
target_integrated_raw
target_mean_raw
background_method
background_value_per_px
target_integrated_background_corrected
target_integrated_intensity_per_DAPI_positive_nucleus
dapi_saturation_fraction
target_saturation_fraction
warnings
```

## `per_nucleus_locations.csv`

Required columns:

```text
image_id
input_path
nucleus_id
x_centroid
y_centroid
area_px
bbox_min_row
bbox_min_col
bbox_max_row
bbox_max_col
touches_border
kept_after_filtering
```

## `roi_level_summary.csv`

Only written when ROI masks are provided.

Required columns:

```text
image_id
roi_id
roi_area_px
dapi_nuclei_count_in_roi
target_integrated_raw_in_roi
target_mean_raw_in_roi
background_method
background_value_per_px
target_integrated_background_corrected_in_roi
target_integrated_intensity_per_DAPI_positive_nucleus_in_roi
warnings
```

## `validation_metrics.csv`

Only written when manual ground truth is provided.

Required columns for centroid validation:

```text
image_id
roi_id
match_radius_px
n_manual
n_predicted
true_positives
false_positives
false_negatives
precision
recall
f1
count_error
count_error_percent
```

Required columns for mask validation:

```text
image_id
roi_id
iou_threshold
n_manual
n_predicted
true_positives
matched_count
false_positives
false_positive_count
false_negatives
false_negative_count
precision
recall
f1
mean_iou_matched
count_error
count_error_percent
```

The current manual instance-mask validator also writes compatibility aliases:

```text
candidate_count
reference_count
```

`candidate_count` is equivalent to `n_predicted`; `reference_count` is
equivalent to `n_manual`.

## `seeded_region_image_metrics.csv`

Exploratory only. Written by `scripts/run_seeded_asma_regions.py`.

Required columns:

```text
image_id
source_id
method
foreground_method
foreground_threshold
background_value_per_px
image_area_px
foreground_area_px
foreground_fraction
seeded_region_area_px
seeded_region_fraction
non_seeded_area_px
non_seeded_area_fraction
dapi_positive_nucleus_count
seeded_region_integrated_raw
seeded_region_mean_raw
seeded_region_integrated_background_corrected
seeded_region_intensity_per_DAPI_positive_nucleus
unseeded_foreground_area_px
unseeded_foreground_fraction
foreground_components
foreground_components_with_seed
qc_status
qc_flags
mask_path
qc_panel_path
warnings
```

Interpretation:

```text
foreground = CH2/aSMA pixels retained by the configured foreground method
seeded_region = foreground components assigned to DAPI nuclei
unseeded_foreground = retained CH2/aSMA foreground components without a DAPI seed
```

`unseeded_foreground` is unresolved target foreground, not proven background.
The seeded-region endpoint is not validated whole-cell expression because CH2 is
both the measured target signal and the foreground-defining signal.

Required warning tokens:

```text
exploratory_asma_associated_region_not_validated_cell_mask
ch2_endpoint_signal_used_to_define_region
```

Common QC flags:

```text
not_validated_whole_cell_mask
low_nucleus_count
sizeable_unseeded_target_foreground
high_unseeded_foreground_fraction
low_fraction_of_foreground_components_with_dapi_seed
near_full_field_seeded_region_fraction
```

## `cellpose_cell_region_image_metrics.csv`

Developer audit table only. Written by `scripts/run_cellpose_cell_regions.py`.

Required columns:

```text
image_id
source_id
method
model_name
background_value_per_px
background_method
image_area_px
cellpose_object_count
candidate_region_area_px
candidate_region_fraction
outside_candidate_region_area_px
outside_candidate_region_fraction
nuclei_mask_path
dapi_positive_nucleus_count
normalization_denominator_count
nuclei_filtering_applied
nuclei_filtering_policy
dapi_nuclei_with_centroid_inside_cellpose_region
dapi_nuclei_centroid_coverage_fraction
cellpose_objects_with_dapi_centroid
cellpose_objects_without_dapi_centroid
cellpose_objects_without_dapi_centroid_fraction
cellpose_objects_with_multiple_dapi_centroids
cellpose_objects_with_multiple_dapi_centroids_fraction
target_integrated_raw_in_cellpose_region
target_mean_raw_in_cellpose_region
target_integrated_background_corrected_in_cellpose_region
target_integrated_intensity_per_DAPI_positive_nucleus
target_integrated_intensity_per_cellpose_object
total_image_integrated_background_corrected
excluded_region_integrated_background_corrected
excluded_region_background_corrected_fraction
excluded_signal_display_threshold
qc_status
qc_flags
mask_path
qc_panel_path
excluded_signal_check_path
warnings
```

Interpretation:

```text
cellpose retained region = Cellpose objects with at least one DAPI-positive
  nucleus centroid inside the object
reported Cellpose intensity = sum(CH2 pixels inside retained regions)
reported target_integrated_intensity_per_DAPI_positive_nucleus =
  retained-region target_integrated_background_corrected /
  dapi_positive_nucleus_count
```

This output must not be presented as validated whole-cell expression unless
manual instance-mask validation later supports that claim. `cellpose_object`
means a candidate Cellpose-labeled object, not a validated biological cell.
Audit columns describing the pre-retention candidate labels and excluded
regions are retained for traceability; they are not the polished Cellpose
endpoint. Current warning tokens must include:

```text
exploratory_cellpose_output_not_validated_whole_cell_mask
ch2_asma_used_as_candidate_cytoplasm_channel
do_not_interpret_as_true_cell_segmentation_without_manual_validation
```

Common QC flags:

```text
not_validated_whole_cell_mask
zero_cellpose_object_count
zero_dapi_positive_nucleus_count
near_full_field_cellpose_region_fraction
very_low_cellpose_region_fraction_with_dapi_nuclei
low_dapi_nuclei_centroid_coverage
candidate_objects_without_dapi_centroid
candidate_objects_with_multiple_dapi_centroids
sizeable_background_corrected_ch2_outside_cellpose_region
majority_background_corrected_ch2_outside_cellpose_region
```

## `watershed_vs_random_walker_*.csv`

Exploratory method comparison table.

Required columns:

```text
source_id
both_empty
watershed_seeded_region_fraction
random_walker_seeded_region_fraction
seeded_region_fraction_delta_random_minus_watershed
watershed_per_nucleus
random_walker_per_nucleus
per_nucleus_delta_random_minus_watershed
union_jaccard
label_agreement_inside_shared_union
watershed_qc_status
random_walker_qc_status
```

`union_jaccard` is `NaN` when both methods produce empty masks; those rows must
not be used as evidence of segmentation robustness.

## `method_triage_summary.csv`

Exploratory cross-method decision table. Written by
`scripts/run_method_triage.py`.

Required columns:

```text
image_id
whole_field_raw_per_nucleus
whole_field_raw_integrated
dapi_nucleus_count
seeded_watershed_per_nucleus
seeded_watershed_region_fraction
seeded_watershed_unseeded_foreground_fraction
seeded_watershed_qc_status
seeded_watershed_qc_flags
seeded_random_walker_per_nucleus
seeded_random_walker_region_fraction
seeded_random_walker_unseeded_foreground_fraction
seeded_random_walker_qc_status
seeded_random_walker_qc_flags
seeded_propagation_per_nucleus
seeded_propagation_region_fraction
seeded_propagation_unseeded_foreground_fraction
seeded_propagation_qc_status
seeded_propagation_qc_flags
cellpose_candidate_per_nucleus
cellpose_candidate_region_fraction
cellpose_excluded_corrected_fraction
cellpose_qc_status
cellpose_qc_flags
missing_sources
any_region_restricted_qc_reject
accepted_region_restricted_method
triage_status
```

Interpretation:

```text
whole_field_raw_per_nucleus = PI workbook full-field CH2 sum / DAPI nuclei count
seeded_*_per_nucleus = CH2 inside exploratory DAPI-seeded aSMA-associated foreground / DAPI nuclei count
seeded_propagation_per_nucleus = CH2 inside CellProfiler-style propagated aSMA-associated foreground / DAPI nuclei count
cellpose_candidate_per_nucleus = CH2 inside exploratory Cellpose candidate regions / DAPI nuclei count
```

`missing_sources` must list absent source tables for a row, using semicolon
separation, rather than filling absent measurements with zero. Region-restricted
methods must remain unaccepted unless manual validation or another explicit
future validation standard supports promotion.

Allowed `triage_status` values:

```text
not_validated_manual_validation_required
manual_review_flagged_not_validated
mixed_region_restricted_qc_rejection_not_validated
all_region_restricted_methods_rejected
region_restricted_sources_missing
manual_validation_available_review_required
```

`manual_validation_available_review_required` does not mean a method is
accepted. It only records that manual validation was declared available; method
promotion still requires explicit review of the validation metrics and
acceptance criteria.

## `sensitivity_long.csv`

Exploratory robustness diagnostic table. Written by
`scripts/run_sensitivity_summary.py`.

Required columns:

```text
run_id
image_id
method
foreground_method
per_nucleus
region_fraction
unseeded_foreground_fraction
qc_status
qc_flags
```

Interpretation:

```text
run_id = user-supplied label for the source seeded-region metrics table
per_nucleus = seeded_region_intensity_per_DAPI_positive_nucleus from that run
region_fraction = seeded_region_fraction from that run
```

This table compares exploratory region-restricted methods only. It must not be
used as manual segmentation validation.

## `sensitivity_run_summary.csv`

Exploratory run-level robustness summary. Written by
`scripts/run_sensitivity_summary.py`.

Required columns:

```text
run_id
n_images
n_reviewable
n_manual_review
n_rejected
expected_order_preserved
challenge_all_rejected
challenge_all_zero_per_nucleus
manual_review_present
hard_reject_present
max_challenge_per_nucleus
```

Interpretation:

```text
expected_order_preserved = all ordered positions are present and strictly descending
challenge_all_rejected = all challenge positions are present and QC-rejected
challenge_all_zero_per_nucleus = all challenge positions are present and report zero retained-region per-nucleus intensity
max_challenge_per_nucleus = largest challenge-field retained-region per-nucleus value
```

`challenge_all_zero_per_nucleus` is intentionally separate from
`challenge_all_rejected`; a rejected method can still emit dangerous numeric
values if QC is ignored.

## `sensitivity_image_summary.csv`

Exploratory image-level robustness summary. Written by
`scripts/run_sensitivity_summary.py`.

Required columns:

```text
image_id
n_runs
min_per_nucleus
max_per_nucleus
mean_per_nucleus
coefficient_of_variation
all_runs_rejected
any_manual_review
any_qc_reject
```

## `manual_validation_manifest.csv`

Written by `scripts/prepare_manual_validation_package.py`.

Required columns:

```text
image_id
source_id
validation_task
ch2_path
ch4_path
candidate_mask_path
nuclei_mask_path
manual_reference_mask_path
annotation_panel_path
guide_panel_path
method
foreground_method
dapi_positive_nucleus_count
candidate_integrated_raw
candidate_intensity_per_DAPI_positive_nucleus
qc_status
qc_flags
```

Interpretation:

```text
candidate_mask_path = automated candidate instance-label mask to validate
manual_reference_mask_path = blank or filled manual/reference instance-label mask
annotation_panel_path = raw-only panel to use while drawing manual/reference labels
guide_panel_path = candidate-overlay panel for post-label QC comparison
validation_task = manual labeling target, e.g. asma_associated_region
```

The manifest does not prove validation. It only records the files required to
produce manual/reference masks and later run the IoU validator.

## `manual_labeling_status.csv`

Written by `scripts/prepare_manual_validation_package.py`.

Required columns:

```text
image_id
manual_reference_mask_path
annotation_panel_path
status
labeler
completed_date
notes
```

Allowed `status` values before validation:

```text
not_started
complete_non_empty
confirmed_empty
```

Interpretation:

```text
not_started = default generated state; do not validate this image
complete_non_empty = manual/reference labels were drawn and positive labels exist
confirmed_empty = human reviewer confirmed no reference objects in this image
```

When this file is present next to a generated validation package, the validator
rejects `not_started` rows and checks that empty masks are marked
`confirmed_empty` while non-empty masks are marked `complete_non_empty`. It also
checks that each row's `manual_reference_mask_path` resolves to the same file as
the reference mask being validated.

## `manual_annotation_audit.csv`

Written by `scripts/audit_manual_annotation_package.py`.

Required columns:

```text
image_id
status
labeler
completed_date
manual_reference_mask_path
annotation_panel_path
reference_mask_exists
annotation_panel_exists
mask_shape
mask_dtype
mask_state
positive_label_count
foreground_area_px
status_mask_consistent
package_has_positive_reference
validation_ready_image
blocking_reasons
```

Interpretation:

```text
validation_ready_image = true only when the status, mask contents, mask path,
annotation panel, full package status set, and package-level positive-reference
state are consistent enough for quantitative validation
blocking_reasons = semicolon-separated reasons preventing validation
```

The audit also writes:

```text
manual_annotation_audit_report.md
manual_annotation_status_contact_sheet.png
```

The audit checks manual annotation readiness only. It does not validate a
candidate segmentation method by itself.

## `annotation_handoff_manifest.csv`

Written by `scripts/prepare_manual_annotation_handoff.py`.

Required columns:

```text
image_id
validation_task
layer_bundle_path
manual_reference_mask_path
annotation_panel_path
guide_panel_path
ch2_path
ch4_path
candidate_mask_path
nuclei_mask_path
status
bundle_shape
instructions
```

Each `layer_bundle_path` points to a compressed NPZ file with:

```text
ch2
ch4
candidate_labels
nuclei_labels
manual_reference_labels
```

The handoff bundle is for annotation convenience only. The authoritative manual
reference label image remains the TIFF recorded in
`manual_reference_mask_path`; after editing, that TIFF and
`manual_labeling_status.csv` must be updated before audit and validation. The
handoff tool rejects status rows whose manual-reference mask path or annotation
panel path does not match the corresponding manifest row for the same image.

## `method_comparison_review_summary.csv`

Written by `scripts/render_method_review_panels.py`.

Required columns:

```text
image_id
dapi_nucleus_count
propagation_region_area_px
cellpose_region_area_px
both_region_area_px
propagation_only_area_px
cellpose_only_area_px
union_area_px
method_region_jaccard
crop_box
propagation_per_DAPI_positive_nucleus
cellpose_per_DAPI_positive_nucleus
propagation_qc_status
cellpose_qc_status
propagation_qc_flags
cellpose_qc_flags
interpretation
```

Interpretation:

```text
propagation = CellProfiler-style DAPI-seeded Otsu/propagation aSMA-associated region
cellpose = Cellpose CH2+CH4 candidate aSMA-associated region
method_region_jaccard = overlap of candidate region masks, not validation IoU
crop_box = y0,x0,y1,x1 crop used in the matched-crop visual panel
interpretation = qualitative_qc_only_not_manual_validation
```

This summary compares two automated candidate methods against each other. It is
not manual validation and must not be reported as segmentation accuracy.

## `validation_candidate_field_features.csv`

Written by `scripts/select_manual_validation_fields.py`.

Required columns:

```text
image_id
source_id
ch2_path
ch4_path
nuclei_mask_path
dapi_positive_nucleus_count
target_integrated_raw
target_integrated_raw_per_DAPI_positive_nucleus
target_saturation_fraction
dapi_saturation_fraction
method_region_jaccard
selection_reasons
```

`selected_manual_validation_fields.csv` uses the same columns, restricted to the
fields chosen for manual/reference annotation.

Interpretation:

```text
target_integrated_raw = sum(all CH2/aSMA pixels)
target_integrated_raw_per_DAPI_positive_nucleus =
  target_integrated_raw / DAPI-positive nucleus count
method_region_jaccard = optional automated method-to-method overlap, not manual validation IoU
selection_reasons = semicolon-separated feature buckets that caused selection
```

The selector defines which images should be manually annotated. It does not
validate any segmentation method.

## `method_validation_summary.csv`

Written by `scripts/run_manual_validation_report.py` after completed
manual/reference masks are available.

Required columns:

```text
candidate_method
n_images
n_manual
n_predicted
true_positives
false_positives
false_negatives
micro_precision
micro_recall
micro_f1
mean_iou_matched
mean_abs_count_error
max_abs_count_error
passes_acceptance_criteria
acceptance_criteria
validation_scope
```

Interpretation:

```text
candidate_method = user-supplied method label from --candidate name=mask_dir
micro_* = aggregate object-level metrics over all evaluated reference images
passes_acceptance_criteria = boolean pass/fail against predeclared thresholds
validation_scope = manual/reference instance-mask task, not biological efficacy
```

## `validation_image_summary.csv`

Written by `scripts/run_manual_validation_report.py`.

Required columns:

```text
candidate_method
image_id
iou_threshold
n_manual
n_predicted
true_positives
false_positives
false_negatives
precision
recall
f1
mean_iou_matched
count_error
count_error_percent
overlay_path
```

`overlay_path` points to a visual QC panel with CH2/aSMA background and
object-match-status colors: yellow for true-positive matched
candidate/reference objects, magenta for false-positive candidate-only object
pixels, green for false-negative manual/reference-only object pixels, and orange
for below-threshold overlap between unmatched candidate and manual/reference
objects. It must be non-empty for every evaluated image. If a candidate mask is
missing for a reference image, the validator counts that case as a
zero-candidate false-negative image and the report still renders an
empty-candidate overlay. These overlays are qualitative review aids only;
quantitative metrics come from the instance-label masks.
