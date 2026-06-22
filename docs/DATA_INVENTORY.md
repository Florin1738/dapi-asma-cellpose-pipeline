# Local Data Inventory

Inventory date: 2026-06-22.

This repository currently contains a local data folder:

```text
ApYYM20AGGSMA_02/
```

The folder is ignored by git and should remain uncommitted.

## Basic Counts

Observed files:

```text
36 tif
13 txt
12 lnk
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
CH2: 960 x 720, 16-bit TIFF
CH4: 960 x 720, 16-bit TIFF
```

## Sidecar and Error Files

The `.lnk` files and `_Error.txt` files appear to be download/shortcut sidecars, not microscopy images. The pipeline should ignore them during image discovery.

`___All_Errors.txt` reports that several files/folders were not downloaded, with a UTC timestamp of `06/22/2026 19:40:58`.

## Open Question

The pipeline still needs confirmation of which channel is DAPI and which channel is the target stain. The file names indicate `CH2` and `CH4`, but channel identity cannot be inferred safely from the names alone.

