# Local Data Inventory

Inventory date: 2026-06-22.

## Current Plate Data

The active plate data are under:

```text
data/aSMA_DAPI_plates/
```

This folder contains Plate 1 and Plate 2 microscopy images organized by plate, acquisition folder, `XY##` field, and channel TIFF. CH2 is confirmed as alpha-smooth muscle actin (aSMA), and CH4 is confirmed as DAPI for the current workflow.

## Legacy Sample Folder

The repository also contains an older local sample folder:

```text
ApYYM20AGGSMA_02/
```

Both data locations are ignored by git and should remain uncommitted.

## Basic Counts

Observed files:

```text
36 tif
13 txt
12 lnk
3 md navigation READMEs added locally
```

Approximate folder size:

```text
87 MB
```

## Image Pattern

The image folder appears to contain 12 XY positions:

```text
XY01
XY02
...
XY12
```

Each observed XY folder contains:

```text
<sample>_XY##_CH2.tif
<sample>_XY##_CH4.tif
<sample>_XY##_Overlay.tif
```

Initial metadata check on `XY01`:

```text
CH2: 960 x 720 x 3 RGB pseudocolor, 16-bit TIFF, plus 120 x 160 thumbnail
CH4: 960 x 720 x 3 RGB pseudocolor, 16-bit TIFF, plus 120 x 160 thumbnail
```

The current TIFF files are RGB rendered exports, not plain grayscale single-channel scientific images. See `docs/DATA_ORGANIZATION.md` for the full organization and risk assessment.

## Sidecar and Error Files

The `.lnk` files and `_Error.txt` files appear to be download/shortcut sidecars, not microscopy images. The pipeline should ignore them during image discovery. The `.md` files are local navigation READMEs added for this project.

`___All_Errors.txt` reports that several files/folders were not downloaded, with a UTC timestamp of `06/22/2026 19:40:58`.

## Open Question

The pipeline still needs confirmation of which channel is DAPI and which channel is the target stain. Current evidence suggests `CH4` is DAPI and `CH2` is target, but that is based on visual morphology and pseudocolor, not embedded metadata.
