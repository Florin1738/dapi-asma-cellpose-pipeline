# Thresholds And Parameters

Last updated: 2026-07-06.

This file records every analysis threshold or parameter that can affect
quantification. Parameters are grouped by status.

## Accepted Current Deliverable

### PI Workbook Raw aSMA Intensity

Status: accepted PI-facing deliverable.

Thresholds:

```text
none
```

Formula:

```text
sum(all CH2 pixels)
```

Rationale:

This directly implements the PI request for "aSMA intensity" in the red/CH2
channel. It is intentionally simple and does not attempt background subtraction
or cell segmentation.

Limitations:

- Includes background, empty field, and non-cell pixels.
- Sensitive to saturation and acquisition/export settings.
- Not a final per-cell expression metric.

## Accepted Current Segmentation

### Cellpose DAPI Nuclei Count

Status: accepted automated count method, pending manual validation.

Parameters:

```text
default channel = CH4
run-selected channel = recorded in channel_id / dapi_channel_id
channel_identity_confirmed = true
model = cpsam_v2
channels = [0, 0]
diameter = None
flow_threshold = 0.4
cellprob_threshold = 0.0
filtering = none
```

Rationale:

- CH4 was user-confirmed as DAPI for the original and 2026-07-04 downloaded
  datasets. The nontechnical runner now also accepts an explicit DAPI channel
  selection and records it in the run logs, CSVs, and workbook metadata.
- Cellpose `cpsam_v2` is the project v1 segmentation backend.
- `diameter=None` lets Cellpose estimate object scale.
- `flow_threshold=0.4` and `cellprob_threshold=0.0` are the current Cellpose
  call defaults used by the project scripts.

Validation status:

- Visual QC overlays exist.
- Manual centroid/mask ground truth has not been supplied, so precision, recall,
  F1, false positive rate, and false negative rate must not be claimed.

## Exploratory/Provisional

### CellProfiler ECM CH1 Positive Mask

Status: current reviewable ECM output for the July 2026 download; not manual
ground-truth validated.

Current output:

```text
reports/cellprofiler_ecm_ch1_2026_07_06_k2_arm64/
```

Parameters:

```text
input ECM channel = CH1
DAPI context channel = CH4
CellProfiler module = IdentifyPrimaryObjects
threshold method = Robust Background
threshold deviations = 2
lower outlier fraction = 0.05
upper outlier fraction = 0.05
averaging method = Mean
variance method = Standard deviation
threshold strategy = Global
threshold correction factor = 1.0
discard objects outside diameter range = No
discard border objects = No
fill holes = Never
background value = median of mask-negative CH1 pixels per image
fallback background = whole-field 10th percentile if no mask-negative pixels exist
normalization denominator = none
DAPI count role = context/QC only; not used for ECM normalization
```

Primary endpoint:

```text
ecm_positive_integrated_background_corrected
  = sum(max(CH1 - background, 0) inside the CellProfiler ECM-positive mask)
```

Parameter selection:

```text
candidate deviations tested in full current run = 2, 3, 5, 8
expanded threshold sensitivity ladder = 1, 1.5, 2, 2.5, 3, 4, 5, 8
representative fields = 24
selected deviations = 2
selection_source = human_qc_override
```

Expanded threshold sensitivity diagnostics:

```text
output root = reports/cellprofiler_ecm_ch1_2026_07_06_threshold_sensitivity/
k=1:   median area fraction 0.1687; median integrated endpoint 2.02e9
k=1.5: median area fraction 0.1344; median integrated endpoint 1.82e9
k=2:   median area fraction 0.1108; median integrated endpoint 1.63e9
k=2.5: median area fraction 0.0923; median integrated endpoint 1.43e9
k=3:   median area fraction 0.0771; median integrated endpoint 1.22e9
k=4:   median area fraction 0.0557; median integrated endpoint 8.95e8
k=5:   median area fraction 0.0160; median integrated endpoint 2.20e8
k=8:   median area fraction 0.0099; median integrated endpoint 1.64e8
```

Rationale:

The first full run selected `k=3` by an automatic/qualitative heuristic that
targeted moderate mask area and low automated QC flags. That was not a
ground-truth validation. Human QC judged `k=3` too aggressive for ECM because
it removed too much weaker CH1 structure. On the representative sweep, `k=2`
increased median ECM-positive area fraction from `0.077` to `0.111` compared
with `k=3` and still avoided near-full-field masks. The `k=2` run therefore
supersedes the July 5 `k=3` output.

Validation status:

Visual QC overlays exist for every image. Manual ECM ground-truth masks are not
available, so precision, recall, F1, false positives, false negatives, and IoU
must not be reported.

Other CellProfiler threshold settings to run before promoting an ECM method to
accepted/final:

```text
Robust Background deviation ladder = 1, 1.5, 2, 2.5, 3, 4, 5, 8
Robust Background correction factor around k=2 = 0.8, 1.0, 1.2
Otsu global threshold
Minimum Cross-Entropy global threshold
Background threshold
Adaptive versions of the strongest global candidates, if illumination varies
Manual threshold only if justified by controls or manual ground truth
```

Rationale for alternatives:

- Robust Background is appropriate when most pixels are background; it can be
  brittle if the ECM-positive fraction varies substantially across images.
- Otsu and Minimum Cross-Entropy are standard CellProfiler alternatives when the
  image histogram better supports separating foreground/background peaks.
- Threshold correction factors test whether the chosen automatic threshold is
  systematically too permissive or too stringent without changing method family.

### Cross-Method Region-Restricted Triage

Status: exploratory method-selection report; not an endpoint.

Current comparison scope:

```text
positions = XY22, XY23, XY24, XY40, XY41
manual_validation_available = false
accepted_region_restricted_method = false for all current rows
missing measurements = blank/NA, not zero-filled
```

Decision rule:

```text
If no region-restricted source rows exist for an image, triage status is
region_restricted_sources_missing.

If all available region-restricted methods for an image are QC-rejected, triage
status is all_region_restricted_methods_rejected.

If at least one available region-restricted method is QC-rejected and at least
one is not, triage status is
mixed_region_restricted_qc_rejection_not_validated.

If at least one available region-restricted method is flagged
needs_manual_review and none are QC-rejected, triage status is
manual_review_flagged_not_validated.

If available region-restricted methods are not QC-rejected but manual validation
is unavailable, triage status is not_validated_manual_validation_required.

If manual validation is declared available, triage status is
manual_validation_available_review_required; this still does not auto-accept a
method.
```

Rationale:

The comparison table is meant to prevent method drift. It keeps the accepted
PI-facing whole-field metric separate from exploratory region-restricted
methods, and it requires manual validation before any region-restricted method
can be promoted.

Representative output:

```text
output/method_triage/plate1_xy22_xy23_xy24_xy40_xy41/
```

### Region-Restricted Raw0 Sensitivity

Status: exploratory robustness diagnostic; not an endpoint.

Current comparison scope:

```text
positions = XY22, XY23, XY24, XY40, XY41
ordered_positions = XY22, XY23, XY24
challenge_positions = XY40, XY41
background_value_per_px = 0 for all compared runs
```

Run-level checks:

```text
expected_order_preserved = true only when XY22 > XY23 > XY24
challenge_all_rejected = true only when every challenge field has qc_status reject_qc_failure
challenge_all_zero_per_nucleus = true only when every challenge field reports 0 retained-region intensity per DAPI-positive nucleus
max_challenge_per_nucleus = maximum retained-region intensity per DAPI-positive nucleus among challenge fields
```

Rationale:

`challenge_all_rejected` alone is not sufficient because a method can correctly
mark challenge fields as QC failures while still emitting large numeric values.
`challenge_all_zero_per_nucleus` records whether the method is numerically safe
if a downstream user accidentally ignores QC status.

Representative output:

```text
output/sensitivity/plate1_xy22_xy23_xy24_xy40_xy41_region_restricted_raw0/
```

### Cellpose Target+DAPI Retained aSMA-Associated Object Segmentation

Status: active full-plate candidate-region analysis; not accepted as
whole-cell segmentation.

Current full-plate parameters:

```text
default method = cellpose_ch2_ch4_candidate_asma_associated_region
model = cpsam_v2
default input channel order = [CH2/aSMA, CH4/DAPI]
run-selected target_channel_id = recorded in run logs and output tables
run-selected dapi_channel_id = recorded in run logs and output tables
channel_axis = 0
diameter = None
flow_threshold = 0.4
cellprob_threshold = 0.0
background_value_per_px = 0.0
whole-cell claim = false
```

Retained-object rule:

```text
For each Cellpose candidate object:
  retain the object if at least one DAPI-positive nucleus centroid falls inside it
  exclude the object if no DAPI-positive nucleus centroid falls inside it
```

This rule is deterministic and is not an arbitrary intensity cutoff. In current
user-facing outputs, `Cellpose` means this retained object set.

QC scoring:

```text
zero Cellpose candidate-object count = reject
zero DAPI-positive nucleus count = reject
Cellpose candidate-region fraction >= 0.90 = reject as near-full-field mask
Cellpose candidate-region fraction <= 0.01 when DAPI nuclei exist = reject
DAPI centroid coverage < 0.50 = reject
candidate objects without DAPI centroid fraction > 0.10 = needs manual review
candidate objects with multiple DAPI centroids fraction > 0.05 = needs manual review
background-corrected CH2 outside candidate regions > 0.20 = needs manual review
background-corrected CH2 outside candidate regions > 0.50 = reject
otherwise status = reviewable_not_validated
```

Rationale:

This tests whether a pretrained learned segmentation model can produce more
plausible aSMA-associated object masks than hand-built threshold/radius methods.
It is practical with the project-local Cellpose environment and the already
cached `cpsam_v2` weights.

For the 2026-07-04 downloaded dataset, visual inspection of representative
fields and the provided overlay mapping showed `CH4` as DAPI, `CH2` as the
red/aSMA target channel, and `CH1` as the green ECM channel not measured in the
current endpoint.

Reason not accepted as final:

- CH2/aSMA is the endpoint signal, not a pan-cell cytoplasm or membrane marker.
- The masks can therefore represent aSMA-positive objects/regions rather than
  all cells.
- Manual instance-mask validation is still required before precision, recall,
  F1, IoU, or whole-cell segmentation accuracy can be reported.

Representative outputs:

```text
output/cellpose_cell_regions/full_plate_cpsam_v2/
reports/cellpose_ch2_ch4_full_plate/
output/cellpose_cell_regions/plate1_ApYYM20AGGSMA_01_xy22_xy23_xy24_cpsam_v2/
output/cellpose_cell_regions/plate1_ApYYM20AGGSMA_01_xy40_xy41_challenge_cpsam_v2/
```

First-pass result:

```text
XY22 target_integrated_intensity_per_DAPI_positive_nucleus = 5.64e7
XY23 target_integrated_intensity_per_DAPI_positive_nucleus = 3.66e7
XY24 target_integrated_intensity_per_DAPI_positive_nucleus = 2.65e7
XY22 qc_status = reviewable_not_validated
XY23 qc_status = reject_qc_failure
XY24 qc_status = reject_qc_failure
XY40 qc_status = reject_qc_failure
XY41 qc_status = reject_qc_failure
```

Full-plate result summary:

```text
fields processed = 227
Plate 1 fields = 111
Plate 2 fields = 116
median whole-field raw CH2/aSMA per DAPI-positive nucleus = 6.00e7
median Cellpose retained-region raw CH2/aSMA per DAPI-positive nucleus = 3.77e7
median fraction of whole-field raw CH2 retained inside Cellpose regions = 66.7%
automated QC status = reviewable_not_validated: 20, needs_manual_review: 163, reject_qc_failure: 44
```

### Nontechnical Cellpose Batch Runner

Status: accepted production wrapper around the current Cellpose DAPI nuclei
count and Cellpose CH2+CH4 retained-region workflow.

Default parameters:

```text
dapi_channel = CH4
target_channel = CH2
model = cpsam_v2
model_cache = .models/cellpose/cpsam_v2
nuclei_channel_identity_confirmed = true
region_input_channel_order = [CH2/aSMA, CH4/DAPI]
region_channel_axis = 0
diameter = None
flow_threshold = 0.4
cellprob_threshold = 0.0
background_value_per_px = 0.0
```

Exported primary endpoint:

```text
target_integrated_intensity_per_DAPI_positive_nucleus =
  raw_retained_CH2_integrated_intensity /
  filtered_DAPI_positive_nucleus_count
```

Rationale:

The runner does not introduce a new segmentation threshold. It packages the
current Cellpose method into one folder-picking workflow and writes the same
auditable logs, masks, CSVs, and QC overlays.

The generated workbook filename retains `background_corrected` for compatibility
with the current report package. Current runner defaults use
`background_value_per_px = 0.0`, so no scalar background is subtracted.

### DAPI-Seeded aSMA-Associated Region Segmentation

Status: implemented exploratory comparator; not accepted as whole-cell
segmentation.

Current reviewable seeded-region parameters:

```text
method family = seeded_intensity_watershed, seeded_intensity_random_walker, or seeded_intensity_propagation
input nucleus labels = existing CH4/DAPI Cellpose labels
CH2 foreground method = per-image Otsu threshold
CH2 foreground minimum connected object size = 128 px
CH2 foreground holes = filled
foreground components without any DAPI seed = excluded
background sensitivity value = 0 raw CH2 counts per pixel
background correction = none for current raw0 sensitivity runs
fixed-radius expansion = false
whole-cell claim = false
```

Earlier Li and fixed-value/background-subtracted runs remain on disk as
exploratory history. They are not the current reviewable threshold setting.

Implemented random-walker comparator:

```text
method = seeded_intensity_random_walker
input nucleus labels = existing CH4/DAPI Cellpose labels
CH2 foreground method = per-image Otsu threshold
random_walker_beta = 30
random_walker_mode = cg_j
random_walker_tol = 1e-5
random_walker_prob_tol = 0.01
flat multi-seed foreground component fallback = nearest-marker watershed on a flat elevation image
fixed-radius expansion = false
whole-cell claim = false
```

Implemented CellProfiler-style propagation comparator:

```text
method = seeded_intensity_propagation
backend = centrosome.propagate
input nucleus labels = existing CH4/DAPI Cellpose labels
CH2 foreground method = per-image Otsu threshold
CH2 foreground minimum connected object size = 128 px
CH2 foreground holes = filled
propagation mask = CH2 foreground OR DAPI nucleus labels
reported/measured labels = propagated labels clipped back to CH2 foreground
reported/measured components = CH2 foreground components that contain a DAPI label
propagation_regularization = 0.05
propagation input image = CH2 scaled by 1st to 99.8th percentile and Gaussian-smoothed
fixed-radius expansion = false
whole-cell claim = false
```

Representative propagation output:

```text
output/seeded_asma_regions/plate1_ApYYM20AGGSMA_01_representative_propagation_otsu_reg005/
```

First-pass result:

```text
XY22 seeded propagation CH2/aSMA per DAPI-positive nucleus = 3.55e7
XY23 seeded propagation CH2/aSMA per DAPI-positive nucleus = 1.47e7
XY24 seeded propagation CH2/aSMA per DAPI-positive nucleus = 1.30e7
XY40 seeded propagation CH2/aSMA per DAPI-positive nucleus = 0.00e+00
XY41 seeded propagation CH2/aSMA per DAPI-positive nucleus = 0.00e+00
XY40 qc_status = reject_qc_failure
XY41 qc_status = reject_qc_failure
```

Additional QC scoring added on 2026-06-28:

```text
low nucleus count warning threshold = 10 nuclei
sizeable unseeded target foreground threshold = 0.10
high unseeded foreground fraction threshold = 0.15
low seeded component coverage threshold = 0.50
large seeded area with low nuclei threshold = seeded area fraction >= 0.20 and nuclei < 10
near-full-field seeded region threshold = seeded area fraction >= 0.80
```

`unseeded_foreground_fraction` means CH2/aSMA foreground that was not assigned
to a DAPI-seeded object because the foreground component did not contain a DAPI
label. It is unresolved target foreground, not proven background. Values >= 0.10
trigger manual review; values >= 0.15 also get the stronger high-unseeded flag.

Sensitivity thresholds tested:

```text
Li threshold
Otsu threshold
Triangle threshold
fixed value threshold = 8100 raw CH2 counts
```

Rationale:

This tests an accepted family of seeded secondary-object approaches without
using a constant radius around each nucleus. It asks: "what CH2/aSMA-associated
region is connected to DAPI nuclei and how much CH2 signal lies inside it?"

Reason not accepted as final:

- CH2/aSMA is the endpoint signal, so using it to define the ROI is circular for
  whole-cell expression.
- Visual QC shows plausible behavior for `XY22`, `XY23`, and `XY24`, but
  artifact-like regions in `XY40` and `XY41` are also segmented.
- Manual labels or a pan-cell/body marker would be required before reporting
  cell-boundary accuracy or using it as final background exclusion.

Sensitivity conclusion:

Otsu was the best compromise among tested thresholds for the representative
fields: it preserved `XY22 > XY23 > XY24` while avoiding the severe over-inclusion
seen with the fixed `8100` threshold and the XY22 undercall seen with triangle.
This makes Otsu the best exploratory seeded-ROI threshold among the tested
options, but it still does not validate whole-cell segmentation because the ROI
remains defined by CH2/aSMA.

Watershed-vs-random-walker conclusion:

Using the same Otsu CH2 foreground mask, the watershed and random-walker
comparators produced identical retained foreground unions and identical
union-level intensity-per-DAPI-positive-nucleus values for the representative
fields. The nucleus-to-pixel labels inside the shared union differed modestly in
dense fields, but the current endpoint does not use per-object label assignment.
Rejected empty-mask fields are marked as `both_empty=true` and
`union_jaccard=NaN` in the comparison table rather than being treated as
evidence of segmentation agreement.
Therefore the practical sensitivity is dominated by the CH2 foreground rule, not
by the seeded partitioning backend.

Representative output:

```text
output/seeded_asma_regions/plate1_ApYYM20AGGSMA_01_representative/
```

### P10 Background Correction

Status: exploratory.

Parameter:

```text
background_percentile = 10
```

Formula:

```text
background_value = percentile_10(CH2 image)
corrected_pixel = max(CH2_pixel - background_value, 0)
```

Rationale:

The 10th percentile is a simple per-image low-intensity baseline when no
cell-free background ROI exists. It avoids choosing a bright/positive threshold,
but it is still a scalar image-wide approximation.

Limitations:

- Does not correct local uneven illumination.
- Does not identify cells.
- Can be distorted if most of an image is occupied by cells or artifacts.

### Alpha-SMA Global Foreground Threshold

Status: exploratory only; do not use as final endpoint without explicit approval.

Parameter:

```text
global Otsu threshold on P10-corrected CH2 pixels sampled across the batch
```

Rationale:

This was used to explore whether excluding dim CH2 pixels changed result
patterns. A global batch threshold is more comparable than per-image thresholds.

Reason not accepted:

It defines aSMA-positive area using CH2 intensity itself. That can be useful for
"aSMA-positive burden" after control-based validation, but it is not the same as
all-cell segmentation and should not be used to exclude background by default.

### DAPI-Derived Cell-Neighborhood Territory

Status: rejected as final method; exploratory diagnostic only.

Prototype parameter:

```text
territory expansion distance = 2 * median DAPI nucleus diameter per image
```

Rationale:

This was a geometry-only attempt to exclude pixels far from DAPI nuclei without
using a CH2 intensity cutoff.

Reason not accepted:

It does not correctly capture the irregular, elongated fibroblast-like cell
bodies seen in the aSMA images. It can over-keep dense regions and under-keep
long cytoplasmic/stress-fiber processes.

### Superseded XY22 Fluorescence Pseudo-OD Diagnostic

Status: superseded exploratory visualization.

Initial transform used for the first fluorescence overlay panel:

```text
RGB_normalized = RGB / dtype_max
RGB_density = -ln(1 - RGB_normalized + eps)
rho = sqrt(R_density^2 + G_density^2 + B_density^2)
theta = atan2(G_density, R_density)
phi = arccos(B_density / rho)
```

Parameters:

```text
eps = 1e-6
theta-phi bins = 256 for full panel, 128 for theta-zoom panel
rho ranges = [0, 0.5), [0.5, 1), [1, 2), [2, 4), [4, 8), [8, infinity)
histogram display = log10(pixel_count + 1)
```

Rationale:

This adaptation made bright fluorescence correspond to higher rho, which made
the 1D rho histogram interpretable.

Reason superseded:

- It was not faithful to the NLTD 2.0 plugin implementation, which uses
  brightfield optical density.
- The theta-phi map was visually uninformative because the microscope overlay
  has red and blue signal but no green channel.
- Keep only as a historical note explaining why the first panel failed.

### Corrected XY22 NLTD-Style OD Debug Panel

Status: exploratory visualization only; not an accepted aSMA quantification
endpoint.

Correction logged on 2026-06-28.

The first `XY22` theta-phi panel was useful only for its 1D rho histogram. It
was not a faithful NLTD 2.0 implementation because it used a fluorescence
pseudo-OD transform rather than the plugin's brightfield OD transform. The
corrected debug panel now:

```text
maps CH2/CH4 fluorescence to a white-background OD-compatible image
applies OD = -ln((RGB8 + 1) / 256)
uses theta = atan2(G_OD, R_OD)
uses phi = arccos(B_OD / rho)
```

Parameters:

```text
fluorescence density eps = 1e-6
theta/phi histogram bins = 256
source image display percentiles = [1, 99.8]
theta display range = 0-5 degrees
phi display range = 0-90 degrees
OD transform = -ln((RGB8 + 1) / 256)
near-zero OD cleanup threshold = 1e-12
angle-valid pixel rule = rho > 1e-12 and finite(theta) and finite(phi)
rho histogram y-axis = log scale
2D histogram display = log10(pixel_count + 1)
```

The corrected panel still shows theta collapse because the image has no green
axis. For this two-channel fluorescence case, `rho` and `phi-versus-rho` are the
interpretable color-space diagnostics; theta-phi TPOMs should not be treated as
useful assay outputs.

Limitations:

- This is not an accepted aSMA quantification endpoint.
- The microscope overlay is pseudocolored RGB, not a raw brightfield image.
- The fixed white-background overlay is a diagnostic transformation, not a
  biological measurement.

### CH2-Only Absorbance-Style Density Histograms

Status: exploratory visualization only; not an accepted aSMA quantification
endpoint.

Use case:

This view is for comparing the CH2/alpha-smooth muscle actin intensity
distribution across wells without using CH4/DAPI or color-combination
theta-phi maps.

Formula:

```text
CH2_fraction = CH2_raw / dtype_max
fixed_red_8bit = round((1 - CH2_fraction) * 255)
CH2_density = -ln((fixed_red_8bit + 1) / 256)
```

Parameters used for the Plate 1 `XY22`, `XY23`, `XY24` panel:

```text
histogram bins = 256
histogram x-range = 0 to -ln(1/256)
histogram y-axis = log scale
image display percentiles = [1, 99.8]
image display scaling = shared across wells in the panel
near-zero density cleanup = 1e-12
background exclusion = none
cell segmentation = none
```

Interpretation:

Higher density corresponds to brighter CH2/aSMA fluorescence after an
OD-compatible diagnostic transform. This should be described as
`absorbance-style density`, not physical absorbance measured by a brightfield
microscope.

### CH2 Raw Intensity Histograms

Status: requested raw-intensity visualization.

Use case:

This view directly shows the microscope-exported CH2/alpha-smooth muscle actin
pixel-value distribution. It should be used when the goal is to inspect raw
intensity rather than a transformed density or absorbance-style view.

Formula:

```text
histogram x-axis = raw CH2 pixel value
histogram y-axis = raw pixel count, linear scale
raw integrated intensity = sum(all CH2 pixels)
```

Parameters used for the Plate 1 `XY22`, `XY23`, `XY24` panel:

```text
histogram bins = 256
histogram x-range = 0 to dtype maximum
histogram y-axis = linear
all pixels included = true
background exclusion = none
cell segmentation = none
CH4/DAPI normalization = none
image display percentiles = [1, 99.8]
image display scaling = shared across wells in the panel
```

Important caveat:

The image thumbnails are display-scaled so they are visible. The histogram,
raw integrated intensity, mean, median, p95, p99, and saturation fraction are
computed from raw CH2 pixel values.

### CH2 Background Candidate Search

Status: exploratory candidate screen; not a confirmed negative-control method.

Use case:

Search all available CH2/aSMA images for fields with low CH2 signal that might
serve as candidate background/negative-like distributions after plate-map
review.

Candidate categories:

```text
blank_like_low_CH2_low_DAPI:
  low CH2 score + low DAPI score

low_aSMA_with_DAPI_present:
  CH4 p95 >= global median CH4 p95
  CH4 saturated fraction <= 0.01
  sorted by lowest CH2 p95, then CH2 p99, then CH2 mean
```

Histogram normalization:

```text
histogram_y = pixels_in_bin / total_pixels_in_image
```

This reports fraction of pixels per bin. It is used for comparing histogram
shapes across images. It does not alter the raw CH2 values on the x-axis.

Current pooled candidate estimates:

```text
blank_like_low_CH2_low_DAPI:
  mode bin center ~8,320
  p50 = 8,445
  p95 = 12,107
  p99 = 18,344

low_aSMA_with_DAPI_present:
  mode bin center ~8,320
  p50 = 8,413
  p95 = 13,445
  p99 = 24,066
```

Interpretation:

The shared mode near `8.3k` raw CH2 counts is a plausible baseline/black-level
estimate. The p95/p99 values are candidate upper background-like ranges, not
proof that all pixels below those values are background.

Curated correction after visual review:

Original candidate panel columns `3` and `4` were rejected because they showed
visible CH2-positive structures. The retained candidates were original columns
`2`, `5`, `6`, and `1`, with original column `2` judged best overall.

```text
best_single_XY12:
  mode bin center = 7,808
  p50 = 7,710
  p95 = 8,453
  p99 = 8,718

top3_XY12_XY41_XY40:
  mode bin center = 8,064
  p50 = 8,156
  p95 = 9,357
  p99 = 10,268

accepted4_XY12_XY41_XY40_XY01:
  mode bin center = 8,320
  p50 = 8,252
  p95 = 9,473
  p99 = 11,144
```

Preferred next-step background level:

```text
CH2 scalar background baseline ~= 8.1k raw counts
```

Use the top-3 mode/median region for background subtraction sensitivity testing.
Use p95/p99 as upper-range sensitivity markers, not as the scalar value to
subtract.

### CH2 8,100 Threshold Retained-Area QC

Status: exploratory threshold QC; not accepted as a segmentation method.

Threshold:

```text
raw CH2 threshold = 8,100
retained mask = CH2 >= 8,100
```

Derived measurements:

```text
retained_area_fraction = count(CH2 >= 8,100) / total_pixels
retained_raw_sum = sum(CH2 where CH2 >= 8,100)
background_subtracted_sum = sum(max(CH2 - 8,100, 0))
```

Observed in representative low/medium/high examples:

```text
low examples retained ~19.6% and ~54.3% of pixels
medium examples retained ~97.2% and ~98.2% of pixels
high examples retained ~99.8% and ~99.9% of pixels
```

Interpretation:

The threshold is useful as an empirical scalar background floor, but it is not a
selective foreground/cell mask. Use `sum(max(CH2 - 8,100, 0))` for
background-subtraction sensitivity testing; do not use `CH2 >= 8,100` alone as a
validated CH2-positive area endpoint.

## Display-Only Parameters

These affect QC image display only, not measurement:

```text
display scaling often uses percentiles such as [1, 99.8]
red/blue/green colors are pseudocolor display choices
green X marks centroid positions
yellow overlays may indicate illustrative high-intensity cues
```

Display parameters must not be confused with measurement thresholds.
