from __future__ import annotations

import re
from pathlib import Path

from scripts.run_target_normalization import app as target_measurement_app
from scripts.validate_target_normalization import app as target_validation_app

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ACTIVE_TERMINOLOGY_PATHS = [
    PROJECT_ROOT / "README.md",
    PROJECT_ROOT / "docs" / "DATA_ORGANIZATION.md",
    PROJECT_ROOT / "docs" / "IMAGE_INPUTS.md",
    PROJECT_ROOT / "docs" / "METHODOLOGY_AUDIT.md",
    PROJECT_ROOT / "docs" / "OUTPUT_CONTRACT.md",
    PROJECT_ROOT / "docs" / "QC_AND_VISUALIZATION.md",
    PROJECT_ROOT / "docs" / "ROADMAP.md",
    PROJECT_ROOT / "docs" / "TARGET_NORMALIZATION_RUN.md",
    PROJECT_ROOT / "docs" / "VALIDATION.md",
    PROJECT_ROOT / "scripts" / "run_target_normalization.py",
    PROJECT_ROOT / "scripts" / "validate_target_normalization.py",
    PROJECT_ROOT / "src" / "dapi_norm" / "__init__.py",
    PROJECT_ROOT / "src" / "dapi_norm" / "target_validation.py",
    PROJECT_ROOT / "pyproject.toml",
]
RISKY_TERMINOLOGY = re.compile(
    "|".join(
        [
            r"DAPI[- ]normalized",
            r"normalized intensity",
            r"normalizes it",
            r"normalized value",
            r"normalized endpoint",
            r"target normalized endpoint",
            r"target-channel intensity normalization",
            r"target-channel normalization",
            r"target-normalization command",
            r"target-normalization CLI",
            r"Validate target-normalization",
            r"Validate the generated normalized",
            r"Required target-normalization artifact",
        ]
    ),
    re.IGNORECASE,
)


def test_target_measurement_cli_describes_per_nucleus_endpoint_without_dapi_normalized_wording():
    help_text = target_measurement_app.info.help or ""

    assert "per DAPI-positive nucleus" in help_text
    assert "normalized" not in help_text.lower()


def test_target_validation_cli_describes_per_nucleus_endpoint_without_dapi_normalized_wording():
    help_text = target_validation_app.info.help or ""

    assert "per-DAPI-positive-nucleus" in help_text
    assert "normalization" not in help_text.lower()


def test_active_docs_and_user_facing_code_avoid_misleading_normalized_endpoint_wording():
    matches: list[str] = []
    for path in ACTIVE_TERMINOLOGY_PATHS:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        for line_number, line in enumerate(text.splitlines(), start=1):
            if RISKY_TERMINOLOGY.search(line):
                matches.append(f"{path.relative_to(PROJECT_ROOT)}:{line_number}:{line.strip()}")

    assert matches == []


def test_current_target_run_docs_do_not_list_unproduced_extra_correlation_artifacts():
    for path in [PROJECT_ROOT / "README.md", PROJECT_ROOT / "docs" / "TARGET_NORMALIZATION_RUN.md"]:
        text = path.read_text(encoding="utf-8")
        current_section = text.split("Current target", maxsplit=1)[-1]
        assert "nucleus_count_vs_raw_and_normalized_intensity.png" not in current_section
        assert "count_intensity_correlations.csv" not in current_section
