# Cellpose Smoke Test

Run date: 2026-06-22.

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
numpy: 2.5.0
scikit-image: 0.26.0
tifffile: 2026.6.1
torch.backends.mps.is_available(): True
torch.cuda.is_available(): False
```

## Model Cache

Cellpose v4.2.1.1 uses `cpsam_v2` for this workflow. The model was downloaded outside the Cellpose runtime with `curl` and cached locally:

```text
.models/cellpose/cpsam_v2
size: 1.1G
sha256: 0f1cc3f7ecdd8a037a57c6c48d9d8921391be4cbce3fa9f13c3e3a2e1253c667
```

The model file is not committed to git.

Rebuild the local model cache if needed:

```bash
mkdir -p .models/cellpose
curl -L --fail --retry 8 --retry-all-errors --retry-delay 5 \
  --connect-timeout 30 \
  -o .models/cellpose/cpsam_v2.download \
  https://huggingface.co/mouseland/cellpose-sam/resolve/main/cpsam_v2
mv .models/cellpose/cpsam_v2.download .models/cellpose/cpsam_v2
shasum -a 256 .models/cellpose/cpsam_v2
```

## Image Extraction Result

The source TIFF is an RGB pseudocolor export. The smoke-test script extracted the active blue RGB component from the primary TIFF series and wrote:

```text
output/cellpose_smoke/ApYYM20AGGSMA_XY01_CH4_candidate_dapi_intensity.tif
```

Observed extracted plane:

```text
shape: 720 x 960
dtype: uint16
min: 0
max: 65535
```

## Command

```bash
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_cellpose_smoke.py \
  --input ApYYM20AGGSMA_02/ApYYM20AGGSMA_02/XY01/ApYYM20AGGSMA_XY01_CH4.tif \
  --output output/cellpose_smoke \
  --model cpsam_v2 \
  --gpu
```

## Result

Cellpose segmentation completed successfully.

```text
label_count: 335
```

Successful outputs:

```text
output/cellpose_smoke/ApYYM20AGGSMA_XY01_CH4_candidate_dapi_intensity.tif
output/cellpose_smoke/ApYYM20AGGSMA_XY01_CH4_cpsam_v2_masks.tif
output/cellpose_smoke/ApYYM20AGGSMA_XY01_CH4_cpsam_v2_smoke_montage.png
```

Visual review: the smoke-test montage shows most bright nuclei outlined, with no obvious whole-field failure. This is a compatibility and logic check, not a validated biological result.

## Important Caveat

This dataset is currently based on RGB rendered TIFF exports. Segmentation QC is useful, but intensity quantification should prefer raw grayscale channel exports or OME-TIFFs if available. `CH4` remains candidate DAPI until channel metadata confirms the stain identity.
