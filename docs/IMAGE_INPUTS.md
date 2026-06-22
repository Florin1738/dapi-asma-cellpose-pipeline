# Image Inputs

## v1 Supported Inputs

The first implementation should support recursive processing of:

- `.tif`
- `.tiff`
- `.ome.tif`
- `.ome.tiff`

Proprietary formats such as `.czi`, `.nd2`, and `.lif` are useful but should not block v1. If those files are central, the safest first step is usually conversion to OME-TIFF using trusted acquisition/export software or a dedicated bioformats conversion tool.

## Axes and Channel Handling

The loader must log the raw array shape, inferred axes, selected DAPI channel, selected target channel, and any Z projection.

Supported layout targets:

- `CYX`
- `YXC`
- `ZCYX`
- `CZYX`
- OME metadata-derived axes when available

If axes are ambiguous, the pipeline must stop with a clear error rather than silently selecting the wrong channel.

## Z-Stack Policy

Default behavior:

```text
z_mode: max
```

Supported planned modes:

- `--z-mode max`
- `--z-mode mean`
- `--z-index <integer>`

The resolved Z policy must be written to `output/logs/config_resolved.yaml`.

## Channel Policy

DAPI channel:

- Used only for nuclei segmentation/counting.
- Must not be used as the signal-normalization intensity.

Target channel:

- Used for integrated raw intensity.
- Used for background-corrected integrated intensity.
- Used for the normalized endpoint after division by filtered DAPI-positive nucleus count.

