# Methodology Run Log Template

Copy this template into `docs/methodology/run_logs/YYYY-MM-DD-short-topic.md`.

````text
# <Short Methodology Log Title>

Date:
Status: accepted | exploratory | rejected | superseded

## Purpose

What question or problem does this method/run address?

## Input Data

- input root:
- plates/runs/images:
- channel assumptions:

## Command Or Script

```bash
<exact command>
```

## Outputs

- tables:
- plots:
- QC overlays:
- logs/configs:

## Parameters And Thresholds

| Name | Value | Rationale | Status |
| --- | ---: | --- | --- |
| example_parameter | 1 | why this value was chosen | accepted/exploratory |

## Sources Consulted

- source title and URL:
- how it affected the method:

## Visual QC

Describe what was inspected and link output paths.

## Validation Status

What has been verified? What requires manual ground truth?

## Caveats

List known failure modes and interpretation limits.

## Decision

Should this method be accepted, revised, kept exploratory, or rejected?
````
