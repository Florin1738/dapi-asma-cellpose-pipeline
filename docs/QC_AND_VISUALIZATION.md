# QC and Visualization

## Goal

The pipeline must make segmentation failures visible. CSV summaries are not enough for microscopy analysis because count errors usually come from missed nuclei, merged nuclei, split nuclei, debris, edge artifacts, or channel mix-ups.

## Required QC Images

For every processed image:

1. DAPI channel with nucleus outlines.
2. Target channel with DAPI-derived nucleus outlines.
3. DAPI channel with numbered centroids.
4. Segmentation montage showing DAPI, target, overlay, and labels.

Planned file names:

```text
output/qc/<image_id>_dapi_nucleus_outlines.png
output/qc/<image_id>_target_with_nucleus_outlines.png
output/qc/<image_id>_dapi_numbered_centroids.png
output/qc/<image_id>_segmentation_montage.png
```

## How I Will Show Outputs Back To You

When the pipeline is implemented and run, I can review outputs in three ways:

1. Attach representative QC PNGs in the Codex response using absolute local paths.
2. Open selected QC PNGs locally for visual inspection.
3. Generate a lightweight HTML QC report with thumbnails, counts, warnings, and links to full-size overlays.

Example response format after a real run:

```markdown
![Segmentation montage](/absolute/path/to/output/qc/sample_segmentation_montage.png)
```

## What To Look For

Review overlays for:

- missed nuclei
- merged nuclei
- split nuclei
- debris counted as nuclei
- nuclei excluded by border filtering
- dense-region failures
- wrong DAPI/target channel selection
- saturated target signal
- target background artifacts

## QC Sampling Strategy

Before trusting batch results, review about 10 representative fields:

- 3 clean/typical fields
- 3 dense fields
- 2 dim or high-background fields
- 2 edge cases

If segmentation quality is poor, tune Cellpose parameters first:

- `diameter`
- `flow_threshold`
- `cellprob_threshold`
- `min_nucleus_area`
- `max_nucleus_area`
- `exclude_border`

## Future HTML Report

A later implementation can write:

```text
output/qc_report.html
```

The report should include image thumbnails, raw and filtered nucleus counts, target normalized endpoint, warnings, and direct links to masks/CSVs.

