# Cellpose Nuclei Count Run

Run date: 2026-06-22.

## Purpose

This run counts nuclei from the candidate DAPI channel in the local `ApYYM20AGGSMA_02` sample data. It is a count-only pipeline check: it writes masks, centroid tables, and QC montages, but does not yet measure target-channel intensity.

## Inputs

Dataset root:

```text
ApYYM20AGGSMA_02
```

Processed images:

```text
ApYYM20AGGSMA_02/ApYYM20AGGSMA_02/XY01/ApYYM20AGGSMA_XY01_CH4.tif
...
ApYYM20AGGSMA_02/ApYYM20AGGSMA_02/XY12/ApYYM20AGGSMA_XY12_CH4.tif
```

`CH4` is treated as `candidate_DAPI`, not confirmed DAPI, because acquisition metadata has not yet been found in the folder.

## Command

```bash
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_cellpose_counts.py \
  --input ApYYM20AGGSMA_02 \
  --output output/cellpose_counts \
  --channel CH4 \
  --model cpsam_v2 \
  --gpu
```

Windows PowerShell equivalent:

```powershell
$env:CELLPOSE_LOCAL_MODELS_PATH = "$PWD\.models\cellpose"
python scripts\run_cellpose_counts.py `
  --input ApYYM20AGGSMA_02 `
  --output output\cellpose_counts `
  --channel CH4 `
  --model cpsam_v2 `
  --gpu
```

Runtime notes:

```text
cellpose: 4.2.1.1
torch: 2.12.1
torch.backends.mps.is_available(): True
torch.cuda.is_available(): False
model: .models/cellpose/cpsam_v2
model_sha256: 0f1cc3f7ecdd8a037a57c6c48d9d8921391be4cbce3fa9f13c3e3a2e1253c667
```

## Outputs

Summary tables:

```text
output/cellpose_counts/summaries/nucleus_counts.csv
output/cellpose_counts/summaries/per_nucleus_locations.csv
```

Run metadata:

```text
output/cellpose_counts/logs/config_resolved.yaml
output/cellpose_counts/logs/run_log.txt
```

Masks:

```text
output/cellpose_counts/masks/XY01_CH4_cpsam_v2_labels.tif
...
output/cellpose_counts/masks/XY12_CH4_cpsam_v2_labels.tif
```

QC:

```text
output/cellpose_counts/qc/XY01_CH4_cpsam_v2_montage.png
...
output/cellpose_counts/qc/XY12_CH4_cpsam_v2_montage.png
output/cellpose_counts/qc_contact_sheet.png
```

## Counts

| Position | Nucleus count |
| --- | ---: |
| XY01 | 335 |
| XY02 | 237 |
| XY03 | 221 |
| XY04 | 345 |
| XY05 | 243 |
| XY06 | 157 |
| XY07 | 201 |
| XY08 | 215 |
| XY09 | 121 |
| XY10 | 65 |
| XY11 | 132 |
| XY12 | 251 |

All rows in `nucleus_counts.csv` include:

```text
channel_identity_confirmed: False
warnings: channel_identity_unconfirmed
```

## Visual QC

The quickest review image is:

```text
output/cellpose_counts/qc_contact_sheet.png
```

Representative full-size montages reviewed:

```text
output/cellpose_counts/qc/XY01_CH4_cpsam_v2_montage.png
output/cellpose_counts/qc/XY04_CH4_cpsam_v2_montage.png
output/cellpose_counts/qc/XY10_CH4_cpsam_v2_montage.png
```

The reviewed montages show Cellpose outlines on most bright nuclei. `XY10` has the lowest count and also visibly has fewer candidate nuclei, which is a useful sanity check. This is still not formal accuracy validation because no manual ground-truth annotations exist yet.

## Validation Checks Completed

Primary command:

```bash
.venv/bin/python scripts/validate_cellpose_counts.py --output output/cellpose_counts
```

Windows PowerShell:

```powershell
python scripts\validate_cellpose_counts.py --output output\cellpose_counts
```

Observed validator output:

```text
summary_rows=12
total_nucleus_count=2523
per_nucleus_rows=2523
mask_counts_match_csv=True
```

Additional checks:

- All 12 `CH4` positions produced a mask and QC montage.
- `per_nucleus_locations.csv` includes `kept_after_filtering`; all current values are `True` because no filtering is implemented in the count-only run.
- `config_resolved.yaml` records Cellpose/PyTorch versions, requested device, model path/hash, channel extraction policy, segmentation parameters, and per-image input/output paths.
- The output tables keep the channel identity warning visible.

## Caveats

- These are rendered RGB TIFF exports, not raw grayscale acquisition channels.
- `CH4` likely corresponds to the blue/DAPI-like channel, but this is not confirmed by metadata.
- Counts have not been filtered by area, edge contact, shape, or manual QC acceptance.
- Target-channel intensity measurement has not been implemented in this count-only run.
