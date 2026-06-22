# Cellpose Smoke Test

Attempt date: 2026-06-22.

## Goal

Run a minimal Cellpose compatibility test on one candidate DAPI image:

```text
ApYYM20AGGSMA_02/ApYYM20AGGSMA_02/XY01/ApYYM20AGGSMA_XY01_CH4.tif
```

## Environment Result

The local project environment was created with `uv` and Python 3.12. Cellpose and PyTorch import successfully.

Observed:

```text
cellpose: 4.2.1.1
torch: 2.12.1
torch.backends.mps.is_available(): True
torch.cuda.is_available(): False
```

## Image Extraction Result

The source TIFF is an RGB pseudocolor export. The one-off smoke-test attempt extracted the active blue RGB component from the primary TIFF series and wrote:

```text
output/cellpose_smoke/XY01_CH4_candidate_dapi_intensity.tif
```

Observed extracted plane:

```text
shape: 720 x 960
dtype: uint16
min: 0
max: 65535
```

## Cellpose Result

Cellpose segmentation did not complete.

Cellpose v4.2.1.1 uses the `cpsam_v2` model by default. That model download is approximately 1.15 GB. The download reached about 739 MB and then failed with:

```text
BrokenPipeError: [Errno 32] Broken pipe
```

No completed model file was left in `.models/cellpose`, so no Cellpose mask or segmentation overlay was produced.

Cellpose v4.2.1.1 built-in model names in the installed package are:

```text
cpsam_v2
cpdino
cpdino-vitb
cpsam
```

The old `nuclei` model name from earlier Cellpose workflows is not exposed as a built-in model in this installed Cellpose 4 package.

## Retry Command

Use this command to retry the smoke test through the reusable script:

```bash
mkdir -p .models/cellpose
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_cellpose_smoke.py \
  --input ApYYM20AGGSMA_02/ApYYM20AGGSMA_02/XY01/ApYYM20AGGSMA_XY01_CH4.tif \
  --output output/cellpose_smoke \
  --model cpsam_v2 \
  --gpu
```

Expected successful outputs:

```text
output/cellpose_smoke/ApYYM20AGGSMA_XY01_CH4_candidate_dapi_intensity.tif
output/cellpose_smoke/ApYYM20AGGSMA_XY01_CH4_cpsam_v2_masks.tif
output/cellpose_smoke/ApYYM20AGGSMA_XY01_CH4_cpsam_v2_smoke_montage.png
```

## Practical Consequence

The Python/Cellpose environment setup is mostly in place, but the model cache is incomplete. Before running Cellpose segmentation, one of these needs to happen:

1. Retry the `cpsam_v2` model download on a stable network.
2. Download the Cellpose model weights manually into `.models/cellpose`.
3. Revisit whether the project should pin Cellpose 3.x to use the older nuclei-oriented workflow described in the original handoff.
4. Run the Cellpose setup on the Windows NVIDIA workstation, where larger model downloads and GPU execution may be more practical.

## Important Caveat

Even after Cellpose runs, this dataset is currently based on RGB rendered TIFF exports. Segmentation QC can be useful, but intensity quantification should prefer raw grayscale channel exports or OME-TIFFs if available.
