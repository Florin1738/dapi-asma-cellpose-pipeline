# Target Per-Nucleus Run

Run date: 2026-06-22.

## Purpose

This run measures provisional candidate-target signal from `CH2` and divides it by the Cellpose candidate-DAPI nucleus count from `CH4`.

Endpoint:

```text
target_integrated_intensity_per_DAPI_positive_nucleus =
    full-field background-corrected CH2 integrated intensity / CH4 Cellpose nucleus count
```

This is provisional because channel identity is inferred from the rendered RGB exports, not confirmed microscope metadata.

## Command

```bash
.venv/bin/python scripts/run_target_normalization.py \
  --input ApYYM20AGGSMA_02 \
  --counts output/cellpose_counts \
  --output output/target_normalization \
  --target-channel CH2 \
  --dapi-channel CH4 \
  --background-percentile 10
```

Validation:

```bash
.venv/bin/python scripts/validate_target_normalization.py --output output/target_normalization
```

Observed validator output:

```text
summary_rows=12
plots_exist=True
qc_overlays_exist=True
formulas_match=True
```

The validator checks that:

- count rows used as denominators are Cellpose `CH4` rows
- referenced masks exist and their nonzero label counts match the reported nucleus counts
- `background_value_per_px` matches the declared percentile method
- corrected integrated intensity and the per-DAPI-positive-nucleus endpoint are recomputed from source images
- image-level and well-level summaries do not diverge
- expected plots and QC overlays exist

## Measurement Details

- Target channel: `CH2`, candidate target stain.
- Nucleus denominator: `CH4` Cellpose labels/counts from `output/cellpose_counts`.
- Measurement area: full image field, not a segmented target ROI.
- Background method: 10th percentile of each `CH2` image.
- Background correction: `clip(CH2 - background, lower=0)`.
- Nucleus filtering: none yet; `raw_nucleus_count == filtered_nucleus_count`.

## Outputs

Tables:

```text
output/target_normalization/summaries/image_level_summary.csv
output/target_normalization/summaries/well_level_summary.csv
```

Plots:

```text
output/target_normalization/plots/normalized_intensity_by_well.png
output/target_normalization/plots/target_integrated_vs_nucleus_count.png
```

The `target_normalization` directory and `normalized_intensity` plot names are
legacy artifact names. They refer to target CH2 intensity divided by
DAPI-positive nucleus count, not normalization by DAPI fluorescence intensity.

QC:

```text
output/target_normalization/qc/<XY##_CH2_target_with_CH4_nucleus_outlines>.png
output/target_normalization/qc_contact_sheet.png
output/cellpose_counts/qc_dapi_centroid_xs/<XY##_CH4_DAPI_with_Cellpose_centroid_Xs>.png
output/cellpose_counts/dapi_centroid_x_contact_sheet.png
```

Run metadata:

```text
output/target_normalization/logs/config_resolved.yaml
output/target_normalization/logs/run_log.txt
```

## Result Ranking

Values are shown in scaled units for readability.

| Position | CH4 nuclei | Corrected CH2 integrated intensity, billions | Corrected CH2 per nucleus, millions |
| --- | ---: | ---: | ---: |
| XY08 | 215 | 11.673 | 54.291 |
| XY07 | 201 | 9.669 | 48.106 |
| XY04 | 345 | 15.959 | 46.259 |
| XY02 | 237 | 10.391 | 43.842 |
| XY09 | 121 | 4.990 | 41.238 |
| XY03 | 221 | 8.667 | 39.216 |
| XY01 | 335 | 11.334 | 33.834 |
| XY06 | 157 | 5.079 | 32.353 |
| XY10 | 65 | 1.044 | 16.064 |
| XY11 | 132 | 1.913 | 14.490 |
| XY05 | 243 | 0.422 | 1.737 |
| XY12 | 251 | 0.401 | 1.598 |

## Visual Checks

The ranking plot shows `XY08`, `XY07`, `XY04`, and `XY02` as the highest provisional target signal per nucleus. `XY05` and `XY12` are much lower, and the target-channel QC contact sheet visually supports that pattern: both fields show weak/diffuse `CH2` signal relative to the high-signal fields.

The count-dependence plot shows that raw corrected full-field `CH2` integrated intensity is moderately related to Cellpose nucleus count across these 12 images:

```text
Pearson r = 0.62
Spearman rho = 0.36
```

The current per-DAPI-positive-nucleus endpoint is much less related to nucleus count:

```text
Pearson r = 0.19
Spearman rho = 0.03
```

This does not prove biological independence from cell density. It only shows that the simple division by nucleus count removed much of the first-order count dependence in this small pilot set.

## Caveats

- `CH2` and `CH4` identities remain unconfirmed.
- The input files are rendered RGB TIFF exports rather than raw grayscale acquisition channels.
- The numerator is full-field `CH2`, not a target-object mask or cell-associated target measurement.
- The denominator uses all Cellpose labels; no area, edge, or manual-QC filtering has been applied.
- Target-channel saturation exists, especially `XY04` at `7.84%` saturated pixels.
