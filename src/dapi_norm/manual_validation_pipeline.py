from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from dapi_norm.manual_annotation_audit import run_manual_annotation_audit
from dapi_norm.manual_validation_report import run_manual_validation_report


def run_manual_validation_pipeline(
    *,
    package_dir: Path,
    candidate_dirs: dict[str, Path],
    output_dir: Path,
    iou_threshold: float,
    min_precision: float = 0.8,
    min_recall: float = 0.8,
    min_f1: float = 0.8,
    min_mean_iou: float = 0.5,
) -> dict[str, Any]:
    package_path = Path(package_dir)
    output_path = Path(output_dir)
    audit_outputs = run_manual_annotation_audit(
        package_dir=package_path,
        output_dir=output_path / "annotation_audit",
    )
    audit_rows = _read_csv(audit_outputs["csv"])
    blocked_rows = [
        row
        for row in audit_rows
        if str(row.get("validation_ready_image", "")).strip().lower() != "true"
    ]
    gate_report = output_path / "validation_gate_report.md"
    if blocked_rows:
        _write_gate_report(
            gate_report,
            validation_ready=False,
            audit_report=audit_outputs["report"],
            blocked_rows=blocked_rows,
            stale_validation_report=output_path / "validation_report",
            validation_outputs=None,
        )
        return {
            "validation_ready": False,
            "audit_csv": audit_outputs["csv"],
            "audit_report": audit_outputs["report"],
            "audit_contact_sheet": audit_outputs["contact_sheet"],
            "gate_report": gate_report,
        }

    validation_outputs = run_manual_validation_report(
        candidate_dirs=candidate_dirs,
        reference_dir=package_path / "reference_masks_to_fill",
        completion_status_path=package_path / "manual_labeling_status.csv",
        manifest_path=package_path / "manual_validation_manifest.csv",
        output_dir=output_path / "validation_report",
        iou_threshold=iou_threshold,
        min_precision=min_precision,
        min_recall=min_recall,
        min_f1=min_f1,
        min_mean_iou=min_mean_iou,
    )
    _write_gate_report(
        gate_report,
        validation_ready=True,
        audit_report=audit_outputs["report"],
        blocked_rows=[],
        stale_validation_report=None,
        validation_outputs=validation_outputs,
    )
    return {
        "validation_ready": True,
        "audit_csv": audit_outputs["csv"],
        "audit_report": audit_outputs["report"],
        "audit_contact_sheet": audit_outputs["contact_sheet"],
        "gate_report": gate_report,
        "method_summary": validation_outputs["method_summary"],
        "image_summary": validation_outputs["image_summary"],
        "validation_report": validation_outputs["report"],
        "overlay_dir": validation_outputs["overlay_dir"],
    }


def _write_gate_report(
    path: Path,
    *,
    validation_ready: bool,
    audit_report: Path,
    blocked_rows: list[dict[str, str]],
    stale_validation_report: Path | None,
    validation_outputs: dict[str, Path] | None,
) -> None:
    lines = [
        "# Manual Validation Gate",
        "",
        f"Validation-ready after audit: `{validation_ready}`",
        f"Audit report: `{audit_report}`",
        "",
    ]
    if blocked_rows:
        lines.extend(
            [
                "Validation report was not run because the manual/reference package is not ready.",
                "",
                "| Image | Status | Blocking reasons |",
                "| --- | --- | --- |",
            ]
        )
        for row in blocked_rows:
            lines.append(
                f"| {row['image_id']} | {row['status']} | {row['blocking_reasons']} |"
            )
        if stale_validation_report is not None and stale_validation_report.exists():
            lines.extend(
                [
                    "",
                    f"Stale validation report directory exists: `{stale_validation_report}`.",
                    "The current gate is blocked; do not use those metrics.",
                ]
            )
    else:
        lines.append("Validation report was run against completed manual/reference masks.")
        if validation_outputs:
            lines.extend(
                [
                    "",
                    f"Method summary: `{validation_outputs['method_summary']}`",
                    f"Image summary: `{validation_outputs['image_summary']}`",
                    f"Validation report: `{validation_outputs['report']}`",
                    f"Overlays: `{validation_outputs['overlay_dir']}`",
                ]
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))
