# Manual Mask Validation

Last updated: 2026-06-29.

This project can now prepare manual-validation packages and validate candidate
instance masks against manual/reference label masks with an IoU threshold. This
is the required quantitative route for claiming cell-region or whole-cell mask
accuracy.

## Preparing A Manual Validation Package

Use this before drawing manual/reference masks. It writes raw-only annotation
panels, candidate-overlay guide panels, an auditable manifest, and blank
full-field `uint32` TIFFs for manual labels.

Current recommended validation-selection output:

```text
output/manual_validation_selection/plate1_ApYYM20AGGSMA_01_stratified_16/
```

Current recommended manual-validation package:

```text
manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/
```

Selected fields:

```text
XY22,XY23,XY24,XY40,XY41,XY01,XY74,XY66,XY08,XY95,XY13,XY10,XY11,XY33,XY16
```

This stratified set covers explicit challenge fields, low/mid/high raw aSMA,
low/high DAPI-positive nucleus count, high per-nucleus artifact risk, high CH2
saturation, and method-disagreement fields.

Example for recreating the current Plate 1 exploratory propagation/Otsu
aSMA-associated region candidate:

```bash
.venv/bin/python scripts/run_seeded_asma_regions.py \
  --input "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01" \
  --counts-root output/pi_simple_summary/cellpose_counts/Plate_1/ApYYM20AGGSMA_01 \
  --output output/seeded_asma_regions/plate1_ApYYM20AGGSMA_01_selected15_propagation_otsu_reg005_raw0 \
  --positions XY22,XY23,XY24,XY40,XY41,XY01,XY74,XY66,XY08,XY95,XY13,XY10,XY11,XY33,XY16 \
  --foreground-method otsu \
  --segmentation-method propagation \
  --propagation-regularization 0.05 \
  --background-value 0 \
  --min-size 128
```

Then create the manual package:

```bash
.venv/bin/python scripts/prepare_manual_validation_package.py \
  --input "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01" \
  --seeded-run output/seeded_asma_regions/plate1_ApYYM20AGGSMA_01_selected15_propagation_otsu_reg005_raw0 \
  --output manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region \
  --positions XY22,XY23,XY24,XY40,XY41,XY01,XY74,XY66,XY08,XY95,XY13,XY10,XY11,XY33,XY16 \
  --iou-threshold 0.5 \
  --task asma_associated_region \
  --force-overwrite-reference-masks
```

Use `--force-overwrite-reference-masks` only when deliberately recreating blank
placeholders before manual labels have been drawn. Without that flag, the
generator refuses to overwrite existing reference TIFFs so completed manual
labels are not accidentally erased.

Package contents:

```text
manual_validation_manifest.csv
manual_labeling_status.csv
README.md
annotation_panels_raw_only/
guide_panels/
reference_masks_to_fill/
```

Use `annotation_panels_raw_only/` while drawing manual labels. These panels show
only CH2/aSMA, CH4/DAPI with green X nucleus centroids, and written annotation
rules. They intentionally do not show the automated candidate mask.

Use `guide_panels/` only after labeling for qualitative comparison. Those panels
show CH2/aSMA, CH4/DAPI with green X nucleus centroids, the automated candidate
retained region over CH2, and QC status/flags. Tracing the candidate overlay
would bias the reference labels.

The blank masks in `reference_masks_to_fill/` are placeholders; they must be
filled before any quantitative validation metrics are reported.

`manual_labeling_status.csv` starts with every image marked `not_started`.
Before validation, every image must be changed to one of:

```text
complete_non_empty = manual/reference labels were drawn and the mask contains objects
confirmed_empty = the field was reviewed and intentionally has no reference objects
```

The validator auto-detects this status file when it sits next to
`reference_masks_to_fill/`, or it can be passed explicitly with
`--completion-status`.

Each status row is also bound to its `manual_reference_mask_path`. The validator
rejects a status file whose recorded mask path does not match the actual
reference TIFF being validated, so status files should not be copied across
packages.

## Committing Edited Manual Labels

Use `scripts/commit_manual_reference_mask.py` after exporting an edited label
image from napari, Fiji/ImageJ, or another label editor. This command validates
the edited labels, writes them to the authoritative package TIFF, and updates
`manual_labeling_status.csv`.

Example for a non-empty annotated field:

```bash
.venv/bin/python scripts/commit_manual_reference_mask.py \
  --package manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region \
  --image-id XY22 \
  --labels path/to/edited_XY22_manual_reference_labels.tif \
  --labeler FS \
  --completed-date 2026-06-29 \
  --notes "manual aSMA-associated regions traced from CH2 with DAPI context"
```

The default `--status auto` writes `complete_non_empty` only when the edited
mask contains at least one positive label. It refuses empty masks so an
accidental blank export is not silently treated as reviewed.

For a field that is intentionally empty, make that decision explicit:

```bash
.venv/bin/python scripts/commit_manual_reference_mask.py \
  --package manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region \
  --image-id XY40 \
  --labels path/to/empty_XY40_manual_reference_labels.tif \
  --status confirmed_empty \
  --labeler FS \
  --completed-date 2026-06-29 \
  --notes "reviewed; no traceable aSMA-associated region"
```

The commit command rejects masks with the wrong shape, non-integer labels,
boolean exports, negative labels, labels that exceed `uint32`, disconnected
objects sharing one label, and status/mask contradictions such as
`confirmed_empty` with positive labels. Integer label images must use distinct
IDs for distinct objects; a single connected object may use label `1`.

For the raw-only export workflow, bulk-import edited scratch TIFFs after the
annotation pass:

```bash
.venv/bin/python scripts/import_raw_annotation_labels.py \
  --package manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region \
  --raw-export-manifest manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/raw_annotation_exports/raw_annotation_export_manifest.csv \
  --output manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/raw_annotation_import \
  --labeler FS \
  --completed-date 2026-06-29 \
  --notes "bulk import from raw-only annotation export"
```

The bulk importer commits non-empty edited labels through
`scripts/commit_manual_reference_mask.py`. Empty edited label TIFFs are skipped
unless the image ID is passed with `--confirm-empty`; use
`--require-all-decisions` when you want the command to fail before any commit if
any field is still empty and unconfirmed.

Before writing any authoritative mask or status file, the bulk importer validates
the raw export manifest and dry-runs every committable edited label. It rejects
duplicate raw manifest image IDs, stale `authoritative_reference_mask_path`
values that do not match the target package, raw-export manifests that indicate
candidate masks were included, and edited labels that would fail the one-off
commit checks. This prevents a later bad TIFF from leaving the package partially
committed.

## What This Validates

The validator compares:

```text
candidate label mask = automated segmentation output
reference label mask = manual human annotation or accepted reference mask
```

Each object must be a 2-D integer label:

```text
0 = background
1, 2, 3, ... = individual objects/regions
```

The validator rejects float/raw-intensity images, RGB/multichannel arrays,
negative labels, and binary or single-label masks that contain multiple
disconnected objects. Binary masks are allowed only if the intended question is
one foreground object versus background. For instance segmentation, each object
must have a separate label.

The validator also rejects an all-empty reference set. This prevents scoring a
fresh package before manual/reference labels have actually been filled.

When a completion-status file is present, the validator rejects any image still
marked `not_started`, any non-empty reference mask not marked
`complete_non_empty`, and any empty reference mask not marked `confirmed_empty`.
This prevents a partially labeled package from being silently treated as fully
reviewed. It also rejects stale or copied status files whose recorded
`manual_reference_mask_path` values do not match the reference masks being
validated.

The report runner applies this completion check to every row in
`manual_labeling_status.csv`, not only the subset of masks discovered in
`reference_masks_to_fill/`.

## Manual Labeling Rules

For the current `asma_associated_region` task:

```text
primary instance unit =
one contiguous aSMA-positive cellular region associated with one or more nearby
DAPI-positive nuclei
```

Operational rules:

- `0` remains background.
- Each manual object/region gets a separate positive integer label.
- Split touching regions when there is a clear intensity valley or visible cell
  boundary.
- Keep inseparable merged regions as one label and record the ambiguity outside
  the mask if needed.
- Do not label CH2-negative DAPI cells as aSMA-associated regions.
- Include edge objects only if most of the aSMA-associated region is visible in
  the field.
- Annotate rejected or challenge fields when visible aSMA-associated regions
  exist; leave a field empty only when the manual reference truly has no object.

This is not the same as validated whole-cell segmentation. For true whole-cell
validation, the reference labeler would need to trace whole-cell boundaries,
which may not be visually supportable from CH2/aSMA plus CH4/DAPI alone.

## IoU Threshold

The initial package records `--iou-threshold 0.5`. This is a pragmatic first
object-overlap threshold for matching candidate and manual instance masks, not a
universal biological cutoff or final acceptance criterion. The threshold is
logged so validation can be rerun at stricter values before promoting a method.

## File Naming

The command matches files by `XY##` in the filename.

Example:

```text
candidate masks:
output/seeded_asma_regions/.../masks/XY22_seeded_intensity_watershed_labels.tif

manual masks:
manual_validation/reference_masks/XY22_manual_cell_region_labels.tif
```

The `XY##` token must be unique within each directory.

## Command

After manual/reference masks are filled:

```bash
.venv/bin/python scripts/validate_manual_instance_masks.py \
  --candidate-dir output/seeded_asma_regions/plate1_ApYYM20AGGSMA_01_selected15_propagation_otsu_reg005_raw0/masks \
  --reference-dir manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/reference_masks_to_fill \
  --completion-status manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/manual_labeling_status.csv \
  --output manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/validation_results \
  --iou-threshold 0.5
```

## Outputs

```text
manual_mask_validation_summary.csv
manual_mask_validation_matches.csv
```

Summary metrics:

```text
n_manual
n_predicted
candidate_count
reference_count
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

Match table:

```text
candidate_label
reference_label
iou
match_status
matched
```

The match table includes true-positive matches plus unmatched false-positive
candidate objects and unmatched false-negative manual/reference objects.

The run also writes:

```text
logs/config_resolved.yaml
logs/run_log.txt
```

These files record the IoU threshold, exact mask paths, evaluated image IDs,
manual/reference images that had no candidate mask, and candidate images that
were present but not evaluated because no manual/reference mask existed.

## Candidate Method Report

After manual/reference masks are filled and `manual_labeling_status.csv` is
updated, use the report runner to compare one or more candidate segmentation
methods under the same reference set:

```bash
.venv/bin/python scripts/run_manual_validation_report.py \
  --candidate propagation_otsu=output/seeded_asma_regions/plate1_ApYYM20AGGSMA_01_selected15_propagation_otsu_reg005_raw0/masks \
  --candidate cellpose_ch2_ch4=output/cellpose_cell_regions/plate1_ApYYM20AGGSMA_01_selected15_cpsam_v2/masks \
  --reference-dir manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/reference_masks_to_fill \
  --completion-status manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/manual_labeling_status.csv \
  --manifest manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/manual_validation_manifest.csv \
  --output manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/method_validation_report \
  --iou-threshold 0.5
```

Repeat `--candidate name=mask_dir` for additional methods. The report writes:

```text
method_validation_summary.csv
validation_image_summary.csv
manual_validation_report.md
per_candidate/<method>/manual_mask_validation_summary.csv
per_candidate/<method>/manual_mask_validation_matches.csv
overlays/<method>/<XY##_candidate_vs_reference_overlay.png>
overlays/<method>_contact_sheet.png
```

Overlays use CH2/aSMA as the background and the object-match table for color:
yellow marks true-positive matched candidate/reference objects, magenta marks
false-positive candidate-only object pixels, green marks false-negative
manual/reference-only object pixels, and orange marks below-threshold overlap
between unmatched candidate and manual/reference objects. These are visual QC
aids; the quantitative decision is still based on the manual/reference masks
and the predeclared acceptance thresholds.

Every evaluated image must have an overlay path. If a candidate method lacks a
mask for an evaluated reference image, the validator counts the image as a
zero-candidate false-negative case and the report renders an empty-candidate
overlay so the failure is still visually reviewable.

The current package cannot produce this report yet because every status row is
still `not_started`.

## Annotation Completeness Audit

Before running quantitative validation, audit the manual package:

```bash
.venv/bin/python scripts/audit_manual_annotation_package.py \
  --package manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region \
  --output manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/annotation_audit
```

Outputs:

```text
annotation_audit/manual_annotation_audit.csv
annotation_audit/manual_annotation_audit_report.md
annotation_audit/manual_annotation_status_contact_sheet.png
```

The current package audit reports:

```text
Overall validation-ready: false
images: 15
validation-ready images: 0
blocked images: 15
not_started: 15
positive-reference images: 0
```

This audit is a readiness check only. It does not validate a candidate method.

## One-Command Validation Gate

After manual/reference labels are committed, use the validation gate instead of
running metrics directly:

```bash
.venv/bin/python scripts/run_manual_validation_pipeline.py \
  --package manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region \
  --candidate propagation=output/seeded_asma_regions/plate1_ApYYM20AGGSMA_01_selected15_propagation_otsu_reg005_raw0/masks \
  --candidate cellpose_ch2_ch4=output/cellpose_cell_regions/plate1_ApYYM20AGGSMA_01_selected15_cpsam_v2/masks \
  --output manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/validation_gate \
  --iou-threshold 0.5
```

Outputs:

```text
validation_gate/validation_gate_report.md
validation_gate/annotation_audit/
validation_gate/validation_report/   # written only if the audit passes
```

The gate always reruns the annotation audit first. If any field is still
`not_started`, has an invalid mask, has a status/mask contradiction, or the
whole package is empty, it writes a blocker table and does not report precision,
recall, F1, false positives, false negatives, or IoU. Once the audit passes, the
same command compares the configured candidate methods and renders validation
overlays.

## Annotation Handoff Bundles

To make manual labeling easier in napari or a similar layer-based editor, create
per-image NPZ layer bundles:

```bash
.venv/bin/python scripts/prepare_manual_annotation_handoff.py \
  --package manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region \
  --output manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/annotation_handoff
```

Outputs:

```text
annotation_handoff/README.md
annotation_handoff/annotation_handoff_manifest.csv
annotation_handoff/layers_npz/<XY##_annotation_layers.npz>
```

Each bundle contains:

```text
ch2
ch4
candidate_labels
nuclei_labels
manual_reference_labels
```

The candidate and nuclei labels are context layers only. The manual annotator
should edit only `manual_reference_labels`, export the completed integer
instance-label image to a scratch TIFF or NPZ, and then run
`scripts/commit_manual_reference_mask.py`. Do not manually overwrite
`reference_masks_to_fill/` or hand-edit `manual_labeling_status.csv`; the commit
tool performs the structural checks and status update. After that, rerun the
annotation audit, and only then run quantitative validation.

The handoff builder treats `manual_validation_manifest.csv` as the package
contract. It refuses to create bundles if `manual_labeling_status.csv` points to
a different manual-reference mask or annotation panel for the same image. This
prevents stale or copied status files from silently directing edits to the wrong
TIFF.

## Raw-Only Editor Export

For first-pass blinded manual/reference annotation, use the raw-only TIFF export
instead of the NPZ handoff:

```bash
.venv/bin/python scripts/prepare_raw_annotation_export.py \
  --package manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region \
  --output manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/raw_annotation_exports
```

Outputs:

```text
raw_annotation_exports/README.md
raw_annotation_exports/raw_annotation_export_manifest.csv
raw_annotation_exports/<XY##>/<XY##_CH2_raw.tif>
raw_annotation_exports/<XY##>/<XY##_CH4_raw.tif>
raw_annotation_exports/<XY##>/<XY##_editable_manual_reference_labels.tif>
```

Candidate masks are intentionally not exported here. This prevents first-pass
manual tracing from being anchored to the candidate segmentation. The scratch
label TIFF can be edited in Fiji, napari, or another label editor, then committed
with the per-image command template in `raw_annotation_export_manifest.csv`
after replacing `YOUR_INITIALS`.

The raw-only export treats `manual_validation_manifest.csv` as the source of
truth. It refuses to export if `manual_labeling_status.csv` points to a
different manual-reference mask for the same image. This prevents stale status
rows from directing edits to the wrong authoritative mask.

After editing the scratch label TIFFs, use
`scripts/import_raw_annotation_labels.py` to commit all non-empty labels and any
explicitly confirmed empty fields in one run. This writes
`raw_annotation_import/bulk_import_summary.csv`, so skipped empty fields remain
visible before the validation gate is run.

## Annotation Review Gallery

Render the static worklist gallery after preparing handoff bundles and running
the annotation audit:

```bash
.venv/bin/python scripts/render_manual_annotation_gallery.py \
  --package manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region \
  --output manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/annotation_review_gallery
```

Open:

```text
manual_validation/plate1_selected15_propagation_otsu_raw0_asma_region/annotation_review_gallery/index.html
```

The gallery is not a validation result. It is a labeling worklist. It requires a
current annotation audit, embeds the raw-only annotation panel for each field,
lists current audit blockers, links to the annotation handoff bundle, and
provides exact safe commit commands for non-empty and intentionally empty
labels. Candidate guide-panel links stay hidden until a field is marked complete
because guide panels are post-label QC only.

## Interpretation Rules

Do not report precision, recall, F1, false positives, false negatives, or cell
boundary accuracy unless this validator has been run against manual/reference
masks.

Validation is anchored to the manual/reference set. If a manual/reference image
has no candidate mask, the validator evaluates that image with a zero-valued
candidate mask and counts all reference objects as false negatives. Candidate
images with no manual/reference mask are logged as not evaluated because they
lack ground truth.

For the current CH2/DAPI-only aSMA images:

```text
validated whole-cell segmentation = no
validated cell/background exclusion = no
exploratory DAPI-seeded aSMA-associated ROI = implemented
manual validation package = generated, awaiting filled reference masks
```

The current blocker is not code. It is the absence of manual/reference labels or
an independent pan-cell/cell-body channel.

## Suggested Manual Label Set

Minimum representative set before considering any method validated:

```text
high aSMA / dense: XY22-like
low aSMA / sparse: XY23-like
moderate aSMA: XY24-like
artifact/low-nucleus field: XY40/XY41-like
high-density positive field: XY95-like
```

For aSMA-associated ROI validation, label the intended aSMA-associated cellular
regions. For true whole-cell segmentation validation, manual labels must trace
whole-cell boundaries, which may not be possible from CH2/DAPI alone.
