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

Matching rule:

```text
A predicted mask matches a manual mask if IoU >= iou_threshold.
Default iou_threshold: 0.5
```

Metrics:

```text
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
- do not use the normalized target endpoint until segmentation error is acceptable for the biological effect size being tested

