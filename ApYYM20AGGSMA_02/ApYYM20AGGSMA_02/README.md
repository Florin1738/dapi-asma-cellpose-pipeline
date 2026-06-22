# Exported TIFF Image Folder

This folder contains the usable exported TIFF images for the current local dataset.

Observed organization:

```text
XY01/
  ApYYM20AGGSMA_XY01_CH2.tif
  ApYYM20AGGSMA_XY01_CH4.tif
  ApYYM20AGGSMA_XY01_Overlay.tif
...
XY12/
  ApYYM20AGGSMA_XY12_CH2.tif
  ApYYM20AGGSMA_XY12_CH4.tif
  ApYYM20AGGSMA_XY12_Overlay.tif
```

Interpretation:

- `XY##` is a field/position identifier, likely related to well/acquisition position.
- `CH2` is red-only in the RGB export and is the candidate target-stain channel.
- `CH4` is blue-only in the RGB export and is the candidate DAPI/nuclei channel.
- `Overlay` combines red and blue and should be treated as a visual display image, not a measurement input.

Do not assume the channel map is confirmed until acquisition metadata or human notes verify it.

