# Advanced Alpha-SMA Fiber Metrics

Last updated: 2026-06-28.

## Why This Exists

Alpha-SMA biology is not always captured by total fluorescence intensity alone.
In myofibroblast activation, organized alpha-SMA stress fibers can be more
phenotypically meaningful than diffuse signal.

## Key Paper To Revisit

Hillsley, Santos, and Rosales 2021:

```text
A deep learning approach to identify and segment alpha-smooth muscle actin
stress fiber positive cells
Scientific Reports 11, Article 21855
https://www.nature.com/articles/s41598-021-01304-4
```

Why it matters:

- The paper frames alpha-SMA stress-fiber organization as a myofibroblast
  phenotype.
- It notes that simple average alpha-SMA intensity can overlook stress-fiber
  organization.
- It describes a computer-vision approach for classifying alpha-SMA
  stress-fiber-positive and stress-fiber-negative cells.

## Candidate Future Metrics

These are future ideas, not current accepted outputs:

- aSMA-positive area fraction, using a pre-specified control-derived threshold
- stress-fiber orientation/coherence
- filamentness or vesselness-style features on CH2
- texture features inside aSMA-positive regions
- fraction of cells or fields with organized aSMA fibers
- deep-learning classification of aSMA stress-fiber-positive objects

## Required Validation Before Use

Before any fiber metric becomes part of the assay:

- define positive and negative control conditions
- create representative manual annotations or review labels
- define field-level exclusion criteria
- test sensitivity to saturation and focus artifacts
- produce visual overlays showing detected fibers/regions
- record all parameters in `THRESHOLDS_AND_PARAMETERS.md`
- log each run in `docs/methodology/run_logs/`

## Current Status

This is a research direction to preserve. It should not affect the current PI
workbook or the current simple quantification outputs until explicitly promoted.
