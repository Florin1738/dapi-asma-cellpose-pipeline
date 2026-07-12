from __future__ import annotations

import csv
from pathlib import Path

from dapi_norm.manual_annotation_gallery import render_manual_annotation_gallery


def test_render_manual_annotation_gallery_writes_status_cards_and_commit_commands(
    tmp_path: Path,
):
    package_dir = _write_gallery_package(tmp_path)

    outputs = render_manual_annotation_gallery(
        package_dir=package_dir,
        output_dir=package_dir / "annotation_review_gallery",
    )

    assert outputs["index"].exists()
    html = outputs["index"].read_text(encoding="utf-8")
    assert "Manual Annotation Gallery" in html
    assert "This gallery is not a validation result" in html
    assert "not_started: 1" in html
    assert "complete_non_empty: 1" in html
    assert "XY01" in html
    assert "status_not_complete" in html
    assert "annotation_panels_raw_only/XY01_manual_annotation_panel.png" in html
    assert "Guide panel is post-label QC only" in html
    assert "guide_panels/XY01_manual_validation_guide.png" not in html
    assert "guide_panels/XY02_manual_validation_guide.png" in html
    assert "scripts/commit_manual_reference_mask.py" in html
    assert "--image-id XY01" in html
    assert "--status confirmed_empty" in html
    assert "&lt;initials&gt;" not in html
    assert "YOUR_INITIALS" in html
    assert "Do not manually overwrite reference_masks_to_fill/" in html


def test_render_manual_annotation_gallery_requires_matching_status_rows(tmp_path: Path):
    package_dir = _write_gallery_package(tmp_path)
    status_path = package_dir / "manual_labeling_status.csv"
    rows = _read_rows(status_path)
    _write_csv(status_path, [row for row in rows if row["image_id"] != "XY02"])

    try:
        render_manual_annotation_gallery(
            package_dir=package_dir,
            output_dir=package_dir / "annotation_review_gallery",
        )
    except ValueError as exc:
        assert "missing status rows" in str(exc)
    else:
        raise AssertionError("expected missing status row failure")


def test_render_manual_annotation_gallery_rejects_non_raw_annotation_panel_path(
    tmp_path: Path,
):
    package_dir = _write_gallery_package(tmp_path)
    manifest_path = package_dir / "manual_validation_manifest.csv"
    rows = _read_rows(manifest_path)
    rows[0]["annotation_panel_path"] = rows[0]["guide_panel_path"]
    _write_csv(manifest_path, rows)

    try:
        render_manual_annotation_gallery(
            package_dir=package_dir,
            output_dir=package_dir / "annotation_review_gallery",
        )
    except ValueError as exc:
        assert "annotation_panels_raw_only" in str(exc)
    else:
        raise AssertionError("expected non-raw annotation panel failure")


def test_render_manual_annotation_gallery_requires_audit_rows(tmp_path: Path):
    package_dir = _write_gallery_package(tmp_path)
    (package_dir / "annotation_audit" / "manual_annotation_audit.csv").unlink()

    try:
        render_manual_annotation_gallery(
            package_dir=package_dir,
            output_dir=package_dir / "annotation_review_gallery",
        )
    except FileNotFoundError as exc:
        assert "manual_annotation_audit.csv" in str(exc)
    else:
        raise AssertionError("expected missing audit failure")


def _write_gallery_package(tmp_path: Path) -> Path:
    package_dir = tmp_path / "manual_validation" / "package"
    annotation_dir = package_dir / "annotation_panels_raw_only"
    guide_dir = package_dir / "guide_panels"
    reference_dir = package_dir / "reference_masks_to_fill"
    audit_dir = package_dir / "annotation_audit"
    for path in [annotation_dir, guide_dir, reference_dir, audit_dir]:
        path.mkdir(parents=True)

    manifest_rows = []
    status_rows = []
    audit_rows = []
    specs = [
        ("XY01", "not_started", "status_not_complete;package_all_reference_masks_empty"),
        ("XY02", "complete_non_empty", ""),
    ]
    for image_id, status, blocking in specs:
        annotation_path = annotation_dir / f"{image_id}_manual_annotation_panel.png"
        guide_path = guide_dir / f"{image_id}_manual_validation_guide.png"
        reference_path = reference_dir / f"{image_id}_manual_reference_labels.tif"
        annotation_path.write_bytes(b"annotation")
        guide_path.write_bytes(b"guide")
        reference_path.write_bytes(b"reference")
        manifest_rows.append(
            {
                "image_id": image_id,
                "source_id": image_id,
                "validation_task": "asma_associated_region",
                "ch2_path": f"data/{image_id}_CH2.tif",
                "ch4_path": f"data/{image_id}_CH4.tif",
                "candidate_mask_path": f"masks/{image_id}_candidate.tif",
                "nuclei_mask_path": f"masks/{image_id}_nuclei.tif",
                "manual_reference_mask_path": str(reference_path),
                "annotation_panel_path": str(annotation_path),
                "guide_panel_path": str(guide_path),
                "method": "propagation",
                "foreground_method": "otsu",
                "dapi_positive_nucleus_count": "12",
                "candidate_integrated_raw": "1234",
                "candidate_intensity_per_DAPI_positive_nucleus": "102.8",
                "qc_status": "reviewable_not_validated",
                "qc_flags": "not_validated_whole_cell_mask",
            }
        )
        status_rows.append(
            {
                "image_id": image_id,
                "manual_reference_mask_path": str(reference_path.resolve()),
                "annotation_panel_path": str(annotation_path.resolve()),
                "status": status,
                "labeler": "tester" if status != "not_started" else "",
                "completed_date": "2026-06-29" if status != "not_started" else "",
                "notes": "",
            }
        )
        audit_rows.append(
            {
                "image_id": image_id,
                "status": status,
                "labeler": status_rows[-1]["labeler"],
                "completed_date": status_rows[-1]["completed_date"],
                "manual_reference_mask_path": str(reference_path.resolve()),
                "annotation_panel_path": str(annotation_path.resolve()),
                "reference_mask_exists": "True",
                "annotation_panel_exists": "True",
                "mask_shape": "8x9",
                "mask_dtype": "uint32",
                "mask_state": "empty" if status == "not_started" else "non_empty",
                "positive_label_count": "0" if status == "not_started" else "1",
                "foreground_area_px": "0" if status == "not_started" else "9",
                "status_mask_consistent": "True",
                "package_has_positive_reference": "True",
                "validation_ready_image": "False" if status == "not_started" else "True",
                "blocking_reasons": blocking,
            }
        )
    _write_csv(package_dir / "manual_validation_manifest.csv", manifest_rows)
    _write_csv(package_dir / "manual_labeling_status.csv", status_rows)
    _write_csv(audit_dir / "manual_annotation_audit.csv", audit_rows)
    return package_dir


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))
