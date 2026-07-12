# Data Folder

This folder contains local microscopy data copied into the project workspace. Treat it as private lab data.

## Current Layout

```text
data/aSMA_DAPI_plates/
  plate 1/
    ApYYM20AGGSMA_01/
    ApYYM20AGGSMA_02/
  plate 2/
    APIPIKEALDSMA/
    APIPIKEALDSMA_01/
```

Within each acquisition folder, fields are organized by `XY##` location, with channel TIFFs such as `*_CH2.tif` and `*_CH4.tif`.

## Channel Assignments

- `CH2`: alpha-smooth muscle actin (aSMA), target stain.
- `CH4`: DAPI, nuclear stain.

## Handling Rules

- Do not commit real microscopy images.
- Do not move data folders without updating configs, docs, and scripts that reference them.
- If new data are added from an external download or shared drive, keep the same plate/acquisition/XY/channel organization whenever possible.
