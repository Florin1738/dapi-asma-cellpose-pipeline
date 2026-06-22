#!/usr/bin/env bash
set -euo pipefail

python scripts/run_pipeline.py \
  --input "/path/to/images" \
  --output "/path/to/output" \
  --dapi-channel 0 \
  --target-channel 1 \
  --backend cellpose \
  --diameter auto \
  --background percentile \
  --background-percentile 10 \
  --save-qc

