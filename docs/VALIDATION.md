# Validation Plan

## Key Rule

Visual QC is not precision/recall. Precision, recall, F1, false positives, and false negatives require manual ground truth.

## Recommended Validation Set

Hand-label approximately 10 representative images before trusting batch results:

- 3 clean/typical fields
- 3 dense fields
- 2 dim or high-background fields
- 2 edge cases

The goal is to estimate segmentation/counting error for the intended biological comparison, not to label every image.

## Centroid-Based Validation

Manual input:

```text
manual_centroids.csv
```

Required columns:

```text
image_id
x
y
```

Optional columns:

```text
roi_id
label
```

Matching rule:

```text
A predicted nucleus centroid matches a manual nucleus if distance <= match_radius_px.
Default match_radius_px: 8
```

Metrics:

```text
true_positives
false_positives
false_negatives
precision
recall
f1
count_error
count_error_percent
```

## Mask-Based Validation

Manual input:

```text
manual label-mask folder
```

Executable validator:

```bash
.venv/bin/python scripts/validate_manual_instance_masks.py \
  --candidate-dir output/seeded_asma_regions/plate1_ApYYM20AGGSMA_01_representative_otsu/masks \
  --reference-dir manual_validation/reference_masks \
  --output output/manual_mask_validation/seeded_otsu_vs_manual \
  --iou-threshold 0.5
```

Matching rule:

```text
A predicted mask matches a manual mask if IoU >= iou_threshold.
Default iou_threshold: 0.5
iou_threshold must be > 0 and <= 1.
```

Metrics:

```text
n_manual
n_predicted
true_positives / matched_count
false_positives / false_positive_count
false_negatives / false_negative_count
precision
recall
f1
mean_iou_matched
count_error
count_error_percent
```

Validation is anchored to the manual/reference set. A manual image with no
candidate mask is counted as zero predicted objects, so all manual objects in
that image become false negatives. Candidate images with no manual/reference
mask are logged as not evaluated because there is no ground truth for them.

The validator also rejects non-2-D masks, non-integer label masks, negative
labels, and binary/single-label masks that contain multiple disconnected
objects. For instance validation, each object needs its own positive integer
label.

## Validation Overlays

For each validation image, save:

```text
output/validation/<image_id>_validation_overlay.png
```

Color convention:

- true positives: green
- false positives: magenta
- false negatives: yellow
- unmatched manual labels: yellow outline
- unmatched predicted labels: magenta outline

Each overlay should include a small legend and the matching settings used.

## Acceptance Thresholds

No universal threshold is scientifically correct for all experiments. A practical starting point:

- inspect whether false positives and false negatives are balanced or biased
- compare count error across experimental groups
- require the same segmentation settings across groups unless there is a documented reason
- do not use the per-DAPI-positive-nucleus target endpoint until segmentation error is acceptable for the biological effect size being tested
