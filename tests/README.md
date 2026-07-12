# Tests Folder

This folder contains pytest tests for the image-analysis package.

## Test Data Policy

Tests use synthetic arrays and temporary files. They should not require private microscopy data.

## Coverage Areas

- data inventory and image extraction
- Cellpose command wrappers and output validation
- PI workbook generation
- target per-DAPI-positive-nucleus calculations
- candidate-region workflows
- QC panel generation
- manual-validation package and report tooling
- terminology guardrails for the scientific contract

Before claiming a code path is complete, run the relevant targeted tests or document why they could not be run.
