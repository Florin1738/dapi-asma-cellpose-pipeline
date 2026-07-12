# Roadmap

## Phase 0: Groundwork

- Document scientific metric and naming.
- Document software choice and install matrix.
- Define output CSV contract.
- Define QC and validation protocol.
- Create project metadata and example config.

## Phase 1: Synthetic Data and Measurement Math

- Generate synthetic two-channel images with known nuclei and target signal.
- Test background correction.
- Test target integrated intensity per DAPI-positive nucleus calculation.
- Test CSV output schema.

## Phase 2: Image I/O

- Discover TIFF/OME-TIFF files recursively.
- Load image arrays and metadata.
- Resolve axes and extract DAPI/target channels.
- Implement Z projection and Z-index selection.
- Fail clearly on ambiguous axes.

## Phase 3: Cellpose Segmentation

- Add Cellpose adapter.
- Save label masks.
- Extract counts, centroids, areas, bounding boxes, and border-touching flags.
- Add area and border filters.
- Log raw and filtered counts.

## Phase 4: Target Measurement

- Measure target integrated raw intensity.
- Add percentile background correction.
- Add saturation checks.
- Compute `target_integrated_intensity_per_DAPI_positive_nucleus`.

## Phase 5: QC Outputs

- Save DAPI outline overlays.
- Save target outline overlays.
- Save numbered centroid images.
- Save segmentation montages.
- Add optional HTML QC report.

## Phase 6: Validation

- Add centroid-based manual validation.
- Add mask-based validation when manual masks are available.
- Save validation overlays and metrics.

## Phase 7: Optional Backends and GUI Interop

- Add StarDist comparison backend after Cellpose is stable.
- Document CellProfiler reference workflow if needed.
- Document QuPath/ROI workflow if whole-slide or tissue-section images become central.
