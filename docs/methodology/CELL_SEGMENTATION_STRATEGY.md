# Cell Segmentation Strategy

Last updated: 2026-06-28.

## Current Constraint

The available images contain:

```text
CH4 = DAPI nuclei
CH2 = alpha-SMA target signal
```

There is no separate pan-cell, membrane, cytoplasm, phalloidin, brightfield, or
phase-contrast channel in the current analysis set.

This matters because robust cell segmentation needs either:

- a reliable cell-body/cell-edge image, or
- a validated model that can infer cells from available channels without
  circularly using the target readout as the cell-definition rule.

## Why The DAPI Expansion Prototype Is Not Final

The DAPI territory prototype expanded each nucleus by a distance based on median
nuclear diameter. It helped explain why whole-field ratios were misleading, but
it does not correctly trace fibroblast-like shapes:

- cells are elongated and irregular
- aSMA stress fibers can extend far from nuclei
- dense fields cause expanded territories to merge or cover too much area
- sparse fields can miss long cell processes

Therefore it is logged as exploratory/rejected for final segmentation.

## CellProfiler Reference Workflow

CellProfiler has the type of workflow we are looking for:

```text
IdentifyPrimaryObjects(DAPI nuclei)
IdentifySecondaryObjects(cells seeded from nuclei)
MeasureObjectIntensity(aSMA image inside cell/cytoplasm objects)
ExportToSpreadsheet
SaveImages/OverlayOutlines for QC
```

The key module is `IdentifySecondaryObjects`. It identifies secondary objects
such as cells using primary objects such as nuclei as seeds. It can use:

- a cell-body or cell-edge image, when available
- a fixed-distance method, when no cell-body image is available
- propagation/watershed-like approaches that combine seed objects and image
  intensity information

## Critical Caveat For Our Images

Using CH2/aSMA as the secondary-object cell-body image is scientifically risky.
It would tend to grow/segment regions where aSMA is bright and fail to represent
aSMA-low cells. That makes the cell mask partly defined by the endpoint we are
trying to measure.

Consequences:

- aSMA-high cells may be overrepresented
- aSMA-low cells may be missed or undersegmented
- per-cell aSMA intensity may be biased upward
- "percent aSMA-positive cells" may be circular unless anchored to controls or
  manual labels

## Candidate CellProfiler Tests

These can be tested visually, but none should become default without review:

### Option A: Distance-N Secondary Objects

Start from DAPI nuclei and expand by a fixed pixel distance.

Pros:

- no CH2 intensity threshold
- reproducible
- available directly in CellProfiler

Cons:

- similar limitation to our DAPI expansion prototype
- poor fit for elongated fibroblast morphology

Status: not preferred as final method.

### Option B: Propagation From DAPI Using CH2/aSMA

Start from DAPI nuclei and use CH2 intensity to help find secondary object
boundaries.

Pros:

- may follow aSMA-positive cell regions better than fixed expansion
- produces object masks that can be visually inspected

Cons:

- circular: aSMA defines the cell ROI and is then measured inside that ROI
- likely misses aSMA-low cells
- should be called "aSMA-associated object measurement," not all-cell
  segmentation

Status: candidate exploratory workflow only.

## Implemented Exploratory Seeded-Region Test

On 2026-06-28 we implemented a project-local analogue of nucleus-seeded,
image-guided secondary-object segmentation:

```text
1. Segment/count CH4 DAPI nuclei with the existing Cellpose nuclei pipeline.
2. Build a CH2/aSMA foreground mask with a data-driven image threshold.
3. Keep only CH2 foreground components that contain at least one DAPI label.
4. Partition those retained CH2 foreground components by marker-controlled
   watershed seeded from the DAPI nucleus labels.
5. Measure CH2 intensity inside the resulting labels.
```

Important limitation: CH2/aSMA foreground components without direct DAPI-label
overlap are not assigned to a seeded region. These pixels are now logged as
`unseeded_foreground` / unresolved target foreground. They are not proven
background, because valid aSMA structures can sit adjacent to nuclei without
overlapping nuclear pixels.

On the same date we added a second comparator using scikit-image random walker:

```text
segmentation_method = random_walker
random_walker_beta = 30
foreground_method = otsu
```

The random-walker comparator is also marker-based and non-radius. It starts from
DAPI labels and uses CH2 image gradients/similarity inside retained CH2
foreground components. Flat foreground components fall back to deterministic
nearest-marker watershed because there is no local intensity information for a
random-walker solver to use.

Important naming rule:

```text
output object = DAPI-seeded aSMA-associated region
output object != whole-cell mask
```

This is not a fixed-radius expansion. It is closer in spirit to
CellProfiler-style seeded secondary-object analysis, but because CH2/aSMA is
the biological endpoint, it remains circular as a cell-definition rule.

Representative Plate 1 results for `XY22`, `XY23`, and `XY24` were useful for
diagnosis:

```text
XY22 seeded-region bg-corrected CH2 per nucleus ~= 3.76e7
XY23 seeded-region bg-corrected CH2 per nucleus ~= 1.65e7
XY24 seeded-region bg-corrected CH2 per nucleus ~= 9.11e6
```

That matches the visual impression that `XY22` has more aSMA signal than
`XY23` or `XY24`. However, the same method segmented obvious artifact-like
regions in low/blank-like fields such as `XY40` and `XY41`. Therefore the method
is useful as a QC/exploratory aSMA-associated ROI comparator, but it is not
robust enough to promote as final cell segmentation.

Watershed-vs-random-walker comparison:

```text
output/seeded_asma_regions/plate1_ApYYM20AGGSMA_01_watershed_vs_random_walker_otsu_beta30_summary.csv
output/seeded_asma_regions/plate1_ApYYM20AGGSMA_01_watershed_vs_random_walker_otsu_beta30_plot.png
```

Result:

```text
retained-region union Jaccard = 1.0 for all non-empty representative fields
union-level intensity-per-nucleus endpoint delta = 0 for all representative fields
label agreement inside the shared union = ~0.84-1.00 for reviewable fields
empty rejected fields = both_empty=true and union_jaccard=NaN
```

Interpretation: random walker changes some nucleus-to-pixel assignment inside
dense retained foreground, but the current measurement sums the whole retained
foreground union. Therefore random walker does not solve the scientific problem
of deciding which CH2 pixels are true cell/aSMA signal versus background. That
decision is controlled by the CH2 foreground definition, not by watershed versus
random-walker partitioning.

### Option C: aSMA Burden And Stress-Fiber Metrics Without All-Cell Masks

Measure aSMA as a field/well phenotype:

- total background-corrected aSMA burden
- aSMA-positive area fraction after a control-defined threshold
- stress-fiber texture/structure metrics
- DAPI count reported separately

Pros:

- matches the available channels
- avoids pretending we have true cell boundaries
- can be scientifically meaningful for an anti-fibrotic assay

Cons:

- not true per-cell expression
- threshold or texture rules need validation against controls/manual review

Status: likely most honest path if no new imaging channels are available.

## Recommendation

Do not claim true cell segmentation from the current images.

Use CellProfiler as a reference/test workflow for secondary objects, but treat
any CH2-guided segmentation as an exploratory aSMA-associated ROI method. For
decision-grade reporting, keep total aSMA burden, DAPI nuclei count, saturation,
and QC visuals separate until a segmentation method is visually and manually
accepted.

The current best decision-grade path remains:

```text
primary: full-field or control/background-corrected CH2/aSMA burden
normalizer/context: DAPI-positive nucleus count
exploratory comparator: DAPI-seeded aSMA-associated regions
not accepted: whole-cell segmentation from CH2+DAPI only
```

## Pending QC Panel Requirement

When method testing resumes, generate an intuitive visual QC panel for each
candidate segmentation method before considering it reviewable. The panel should
include, at minimum:

- CH4/DAPI image with green X markers at the counted nucleus centroids.
- CH2/aSMA image with candidate segmented/retained regions overlaid in light
  green with enough transparency that excluded CH2 signal remains visible.
- A combined CH2 + DAPI centroid view so reviewers can tell whether excluded
  CH2-looking structures lack nearby DAPI nuclei or were missed by the method.
- Per-field labels showing location, method name, raw CH2 integrated intensity,
  retained-region CH2 integrated intensity if applicable, DAPI-positive nucleus
  count, and target integrated intensity per DAPI-positive nucleus.
- Explicit visual distinction between excluded CH2 signal, retained seeded
  regions, DAPI nuclei, and unseeded foreground.

The purpose is qualitative error detection: reviewers should be able to scan for
obvious false exclusions, artifact retention, missed nuclei, and cases where
cell-like/aSMA-like structures were excluded only because the current method
could not associate them with a DAPI-positive nucleus.
