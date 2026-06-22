# Data Organization: `ApYYM20AGGSMA_02`

Inventory date: 2026-06-22.

## Summary

The local dataset is stored at:

```text
ApYYM20AGGSMA_02/
```

The usable image files are in the nested folder:

```text
ApYYM20AGGSMA_02/ApYYM20AGGSMA_02/
```

There are 12 positions, named `XY01` through `XY12`. Each position has one `CH2` TIFF, one `CH4` TIFF, and one `Overlay` TIFF.

## Current Structure

```text
ApYYM20AGGSMA_02/
  README.md
  ___All_Errors.txt

  ApYYM20AGGSMA_02/
    README.md
    H02.lnk
    H03.lnk
    ...
    H12.lnk
    H06(2).lnk

    XY01/
      ApYYM20AGGSMA_XY01_CH2.tif
      ApYYM20AGGSMA_XY01_CH4.tif
      ApYYM20AGGSMA_XY01_Overlay.tif
    ...
    XY12/
      ApYYM20AGGSMA_XY12_CH2.tif
      ApYYM20AGGSMA_XY12_CH4.tif
      ApYYM20AGGSMA_XY12_Overlay.tif

  __ApYYM20AGGSMA_02/
    README.md
    XY01/_H02_Error.txt
    ...
    XY12/_H06_Error.txt
```

## File Counts

Generated inventory at `output/data_inventory/dataset_summary.md` found:

```text
36 .tif
13 .txt
12 .lnk
3 .md local navigation READMEs
```

The image positions are complete in the exported TIFF tree: 12 positions with 3 TIFF files each.

## TIFF Format

These TIFFs are not clean single-channel grayscale scientific images. They are RGB pseudocolor exports:

- Primary image series: `960 x 720 x 3`, axes `YXS`, `uint16`.
- Secondary thumbnail series: `120 x 160 x 3`, `uint8`.
- Compression: LZW.
- No OME metadata, ImageJ metadata, channel names, wavelength labels, DAPI strings, target-stain names, or filter labels were found.

Observed RGB organization:

- `CH2` contains signal in the red RGB component.
- `CH4` contains signal in the blue RGB component.
- `Overlay` combines red and blue.

## Candidate Channel Map

Based on visual morphology and RGB export colors:

| Filename channel | RGB component | Candidate role | Confidence |
|---|---|---|---|
| `CH2` | red | target stain | likely, not metadata-confirmed |
| `CH4` | blue | DAPI nuclei | likely, not metadata-confirmed |
| `Overlay` | red + blue | visual overlay only | confirmed as derived display image |

The candidate DAPI assignment is based on CH4 showing compact punctate nuclear-looking objects. CH2 shows broad fibrous/cytoplasmic-looking signal and is therefore the candidate target channel.

Do not treat this as auditable channel identity until a human confirms the acquisition channel map.

## Sidecars and Download State

The `.lnk` files point back to a Keyence-style Windows source path:

```text
D:\Keyence_Users_D\963\6-22-2026\ApYYM20AGGSMA_02\XY##
```

The `__ApYYM20AGGSMA_02/` folder is not image data. It contains download error text files. `___All_Errors.txt` reports that multiple `_H##` folders/files were not downloaded because of `WebMeTAException`, correlation ID `daec379b-29f5-490d-8ab8-832e81a9b70d`, timestamp `06/22/2026 19:40:58 UTC`.

There is an `H06` duplication:

- `H06.lnk` points to `XY05`.
- `H06(2).lnk` points to `XY12`.
- `XY12/_H06_Error.txt` also uses `_H06`, so the source export/download naming likely reused an H index.

## Quantification Risk

The current files may be display-rendered RGB exports rather than raw fluorescence intensity images. That matters because rendered exports can be clipped, contrast-scaled, pseudocolored, or otherwise transformed.

For serious quantification, prefer raw grayscale channel exports or OME-TIFFs with channel metadata. These current TIFFs can still be useful for orientation, Cellpose smoke tests, and visual QC development, but normalized intensity results from them should be treated as provisional until raw/export settings are confirmed.

## Generated Navigation Outputs

The reproducible inventory command is:

```bash
.venv/bin/python scripts/inspect_dataset.py \
  --root ApYYM20AGGSMA_02 \
  --output output/data_inventory \
  --preview-position XY01 \
  --preview-position XY04 \
  --preview-position XY08 \
  --preview-position XY12
```

Generated files:

```text
output/data_inventory/dataset_summary.md
output/data_inventory/image_manifest.csv
output/data_inventory/channel_interpretation_manifest.csv
output/data_inventory/previews/XY01_channels_preview.png
output/data_inventory/previews/XY04_channels_preview.png
output/data_inventory/previews/XY08_channels_preview.png
output/data_inventory/previews/XY12_channels_preview.png
```
