# Data QC Review

Generated from local data on 2026-06-22.

These preview montages are for channel-role review and navigation. They are not final segmentation outputs. The generated PNGs come from real microscopy data and are intentionally kept under ignored `output/` paths rather than committed to git.

## Interpretation

Across reviewed fields, `CH4` shows compact punctate objects consistent with nuclei, while `CH2` shows broader fibrous/cytoplasmic-looking signal. This supports the provisional map:

```text
CH4 = candidate DAPI/nuclei channel
CH2 = candidate target-stain channel
Overlay = display/QC image only
```

This remains unconfirmed because the TIFFs do not contain channel names, wavelength labels, OME metadata, or ImageJ metadata.

## Representative Fields

### XY01

```text
output/data_inventory/previews/XY01_channels_preview.png
```

### XY04

```text
output/data_inventory/previews/XY04_channels_preview.png
```

### XY08

```text
output/data_inventory/previews/XY08_channels_preview.png
```

### XY12

```text
output/data_inventory/previews/XY12_channels_preview.png
```

## Logic Checks

- If DAPI was expected to be blue in the exported overlay, `CH4` is the correct candidate.
- If the target stain was expected to be red in the exported overlay, `CH2` is the correct candidate.
- If the microscope/acquisition notes say DAPI was not channel 4, stop and update the channel map before analysis.
- If raw grayscale exports are available, prefer those over the current RGB rendered exports for intensity quantification.
