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

## Current Target-Normalization Output

The current target-normalization command is a provisional full-field measurement. It writes:

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
false_positives
false_negatives
precision
recall
f1
mean_iou_matched
count_error
count_error_percent
```
