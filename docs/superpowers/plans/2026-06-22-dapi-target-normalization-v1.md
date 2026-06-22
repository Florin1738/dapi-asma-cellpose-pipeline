# DAPI Target Normalization V1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Cellpose-first CLI that counts DAPI-positive nuclei and reports background-corrected target integrated intensity per DAPI-positive nucleus.

**Architecture:** Implement a small Python package under `src/dapi_norm` with separate modules for I/O, channel extraction, segmentation, measurement, QC, validation, and CLI orchestration. Use synthetic data to test the pipeline before real microscopy data.

**Tech Stack:** Python 3.12, tifffile, numpy, scipy, scikit-image, pandas, matplotlib, typer, pyyaml, pytest, Cellpose, PyTorch.

---

## File Structure

- Create `src/dapi_norm/io.py` for discovery, TIFF loading, and image metadata.
- Create `src/dapi_norm/channels.py` for axis resolution and DAPI/target extraction.
- Create `src/dapi_norm/segmentation.py` for the Cellpose adapter and label-mask filtering.
- Create `src/dapi_norm/measure.py` for background correction, saturation metrics, and target endpoint math.
- Create `src/dapi_norm/qc.py` for overlay and montage generation.
- Create `src/dapi_norm/validation.py` for centroid and mask validation.
- Create `src/dapi_norm/cli.py` for command-line parsing and run orchestration.
- Create `scripts/run_pipeline.py` as a thin CLI entry point.
- Create `scripts/generate_synthetic_test_data.py` for private-data-free test images.
- Create `tests/` for unit and end-to-end synthetic tests.

## Task 1: Synthetic Measurement Baseline

**Files:**
- Create: `scripts/generate_synthetic_test_data.py`
- Create: `tests/test_measurement_math.py`
- Create: `src/dapi_norm/measure.py`

- [ ] **Step 1: Write failing measurement tests**

```python
import numpy as np

from dapi_norm.measure import measure_target_signal


def test_target_signal_per_nucleus_with_percentile_background():
    target = np.array([[10, 10, 10], [10, 50, 50], [10, 50, 50]], dtype=np.float32)
    result = measure_target_signal(
        target=target,
        mask=np.ones_like(target, dtype=bool),
        nucleus_count=4,
        background_method="percentile",
        background_percentile=0,
    )

    assert result.target_integrated_raw == 250.0
    assert result.background_value_per_px == 10.0
    assert result.target_integrated_background_corrected == 160.0
    assert result.target_integrated_intensity_per_DAPI_positive_nucleus == 40.0
```

- [ ] **Step 2: Run test to verify it fails**

Run:

```bash
pytest tests/test_measurement_math.py -v
```

Expected: fails because `dapi_norm.measure` is not implemented.

- [ ] **Step 3: Implement measurement dataclass and function**

Implement `TargetMeasurement` and `measure_target_signal()` in `src/dapi_norm/measure.py` with explicit handling for zero nuclei by returning `nan` and warning text.

- [ ] **Step 4: Run the test to verify it passes**

Run:

```bash
pytest tests/test_measurement_math.py -v
```

Expected: pass.

## Task 2: TIFF Discovery and Channel Extraction

**Files:**
- Create: `src/dapi_norm/io.py`
- Create: `src/dapi_norm/channels.py`
- Create: `tests/test_channel_extraction.py`

- [ ] **Step 1: Write tests for CYX, YXC, ZCYX, and CZYX extraction**
- [ ] **Step 2: Implement explicit axis-aware extraction**
- [ ] **Step 3: Add clear errors for ambiguous axes**
- [ ] **Step 4: Run channel tests**

Run:

```bash
pytest tests/test_channel_extraction.py -v
```

Expected: pass for explicit axes and fail clearly for ambiguous axes.

## Task 3: Cellpose Segmentation Adapter

**Files:**
- Create: `src/dapi_norm/segmentation.py`
- Create: `tests/test_segmentation_filtering.py`

- [ ] **Step 1: Test object filtering on synthetic label masks**
- [ ] **Step 2: Implement filtering by min area, max area, and border contact**
- [ ] **Step 3: Implement Cellpose adapter behind a function named `segment_nuclei_cellpose`**
- [ ] **Step 4: Keep Cellpose import lazy so base tests can run without Cellpose installed**
- [ ] **Step 5: Run filtering tests without requiring Cellpose**

Run:

```bash
pytest tests/test_segmentation_filtering.py -v
```

Expected: pass.

## Task 4: CSV Output Contract

**Files:**
- Create: `src/dapi_norm/results.py`
- Create: `tests/test_output_contract.py`

- [ ] **Step 1: Encode required CSV column lists from `docs/OUTPUT_CONTRACT.md`**
- [ ] **Step 2: Test exact image-level and per-nucleus column order**
- [ ] **Step 3: Implement CSV writing helpers**
- [ ] **Step 4: Run output contract tests**

Run:

```bash
pytest tests/test_output_contract.py -v
```

Expected: pass.

## Task 5: QC Overlay Rendering

**Files:**
- Create: `src/dapi_norm/qc.py`
- Create: `tests/test_qc_outputs.py`

- [ ] **Step 1: Test that QC functions write non-empty PNGs**
- [ ] **Step 2: Implement DAPI outlines, target outlines, numbered centroids, and montage**
- [ ] **Step 3: Run QC tests**

Run:

```bash
pytest tests/test_qc_outputs.py -v
```

Expected: pass and create temporary PNGs during test execution.

## Task 6: Centroid Validation

**Files:**
- Create: `src/dapi_norm/validation.py`
- Create: `tests/test_centroid_validation.py`

- [ ] **Step 1: Write centroid matching tests with known true positives, false positives, and false negatives**
- [ ] **Step 2: Implement one-to-one nearest-neighbor matching within `match_radius_px`**
- [ ] **Step 3: Compute precision, recall, F1, and count error**
- [ ] **Step 4: Run validation tests**

Run:

```bash
pytest tests/test_centroid_validation.py -v
```

Expected: pass.

## Task 7: CLI Orchestration

**Files:**
- Create: `src/dapi_norm/cli.py`
- Create: `scripts/run_pipeline.py`
- Create: `tests/test_synthetic_pipeline.py`

- [ ] **Step 1: Generate synthetic images in a temporary folder**
- [ ] **Step 2: Run the CLI against synthetic images with a lightweight test segmentation mode or mocked Cellpose adapter**
- [ ] **Step 3: Verify masks, CSVs, logs, and QC PNGs exist**
- [ ] **Step 4: Run end-to-end synthetic test**

Run:

```bash
pytest tests/test_synthetic_pipeline.py -v
```

Expected: pass without private microscopy data.

## Task 8: Real Cellpose Smoke Test

**Files:**
- Modify: `docs/INSTALLATION.md`
- Modify: `README.md`

- [ ] **Step 1: Install Cellpose in the project environment after user approval**
- [ ] **Step 2: Run the synthetic pipeline with `--backend cellpose`**
- [ ] **Step 3: Confirm model download path and record package versions**
- [ ] **Step 4: Inspect one QC montage visually**
- [ ] **Step 5: Document exact environment and command used**

Run:

```bash
python scripts/run_pipeline.py \
  --input output/synthetic/input \
  --output output/synthetic/cellpose_run \
  --dapi-channel 0 \
  --target-channel 1 \
  --backend cellpose \
  --diameter auto \
  --background percentile \
  --background-percentile 10 \
  --save-qc
```

Expected: writes summaries, masks, and QC overlays.

