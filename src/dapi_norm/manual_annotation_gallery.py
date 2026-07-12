from __future__ import annotations

from collections import Counter
import csv
from html import escape
import os
from pathlib import Path
import shlex


def render_manual_annotation_gallery(*, package_dir: Path, output_dir: Path) -> dict[str, Path]:
    package_path = Path(package_dir)
    output_path = Path(output_dir)
    manifest_rows = _read_csv(package_path / "manual_validation_manifest.csv")
    status_rows = _read_csv(package_path / "manual_labeling_status.csv")
    audit_rows = _read_csv(package_path / "annotation_audit" / "manual_annotation_audit.csv")
    _validate_manifest_status_match(manifest_rows, status_rows)
    _validate_manifest_audit_match(manifest_rows, audit_rows)

    status_by_image = {_image_id(row): row for row in status_rows}
    audit_by_image = {_image_id(row): row for row in audit_rows}
    output_path.mkdir(parents=True, exist_ok=True)
    index_path = output_path / "index.html"
    index_path.write_text(
        _gallery_html(
            package_dir=package_path,
            output_dir=output_path,
            manifest_rows=manifest_rows,
            status_by_image=status_by_image,
            audit_by_image=audit_by_image,
        ),
        encoding="utf-8",
    )
    return {"index": index_path}


def _gallery_html(
    *,
    package_dir: Path,
    output_dir: Path,
    manifest_rows: list[dict[str, str]],
    status_by_image: dict[str, dict[str, str]],
    audit_by_image: dict[str, dict[str, str]],
) -> str:
    status_counts = Counter(row.get("status", "blank") or "blank" for row in status_by_image.values())
    cards = [
        _image_card(
            row=row,
            status_row=status_by_image[_image_id(row)],
            audit_row=audit_by_image.get(_image_id(row), {}),
            package_dir=package_dir,
            output_dir=output_dir,
        )
        for row in manifest_rows
    ]
    return "\n".join(
        [
            "<!doctype html>",
            '<html lang="en">',
            "<head>",
            '<meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width, initial-scale=1">',
            "<title>Manual Annotation Gallery</title>",
            "<style>",
            _css(),
            "</style>",
            "</head>",
            "<body>",
            "<header>",
            "<h1>Manual Annotation Gallery</h1>",
            f"<p><strong>Package:</strong> {escape(str(package_dir))}</p>",
            "<p class=\"warning\">This gallery is not a validation result. It is a manual-labeling worklist. Quantitative validation remains unavailable until reference masks are completed and audited.</p>",
            "<div class=\"counts\">"
            + " ".join(
                f"<span>{escape(status)}: {count}</span>"
                for status, count in sorted(status_counts.items())
            )
            + "</div>",
            "</header>",
            "<main>",
            "\n".join(cards),
            "</main>",
            "</body>",
            "</html>",
        ]
    )


def _image_card(
    *,
    row: dict[str, str],
    status_row: dict[str, str],
    audit_row: dict[str, str],
    package_dir: Path,
    output_dir: Path,
) -> str:
    image_id = _image_id(row)
    status = status_row.get("status", "blank") or "blank"
    annotation_path = _resolve_path(row.get("annotation_panel_path", ""), package_dir=package_dir)
    guide_path = _resolve_path(row.get("guide_panel_path", ""), package_dir=package_dir)
    _validate_raw_annotation_panel_path(image_id, annotation_path)
    layer_bundle = package_dir / "annotation_handoff" / "layers_npz" / f"{image_id}_annotation_layers.npz"
    blocking = audit_row.get("blocking_reasons", "")
    positive_labels = audit_row.get("positive_label_count", "")
    foreground_area = audit_row.get("foreground_area_px", "")
    return "\n".join(
        [
            f'<section class="card status-{escape(status)}">',
            "<div class=\"card-head\">",
            f"<h2>{escape(image_id)}</h2>",
            f"<span class=\"status\">{escape(status)}</span>",
            "</div>",
            "<div class=\"meta\">",
            _meta_row("task", row.get("validation_task", "")),
            _meta_row("QC status", row.get("qc_status", "")),
            _meta_row("nuclei", row.get("dapi_positive_nucleus_count", "")),
            _meta_row("labels", positive_labels),
            _meta_row("foreground px", foreground_area),
            _meta_row("blocking", blocking or "none recorded"),
            "</div>",
            "<div class=\"panels\">",
            _linked_image(annotation_path, output_dir=output_dir, label="raw-only annotation panel"),
            "</div>",
            "<p class=\"note\">Guide panel is post-label QC only; use the raw-only annotation panel while drawing reference labels.</p>",
            "<div class=\"commands\">",
            "<p><strong>Safe commit command after exporting edited labels:</strong></p>",
            f"<pre>{escape(_commit_command(package_dir, image_id))}</pre>",
            "<p><strong>Confirm intentionally empty field:</strong></p>",
            f"<pre>{escape(_confirm_empty_command(package_dir, image_id))}</pre>",
            "</div>",
            "<div class=\"paths\">",
            _path_link("annotation panel", annotation_path, output_dir=output_dir),
            _guide_link_for_status(status, guide_path, output_dir=output_dir),
            _path_link("layer bundle", layer_bundle, output_dir=output_dir),
            _path_link(
                "reference mask",
                _resolve_path(status_row.get("manual_reference_mask_path", ""), package_dir=package_dir),
                output_dir=output_dir,
            ),
            "</div>",
            "<p class=\"note\">Do not manually overwrite reference_masks_to_fill/ or hand-edit manual_labeling_status.csv. Export edited labels to a scratch file, then run the commit command.</p>",
            "</section>",
        ]
    )


def _linked_image(path: Path, *, output_dir: Path, label: str) -> str:
    rel = _relative_href(path, output_dir=output_dir)
    return (
        f'<figure><a href="{escape(rel)}"><img src="{escape(rel)}" alt="{escape(label)}">'
        f"</a><figcaption>{escape(label)}</figcaption></figure>"
    )


def _path_link(label: str, path: Path, *, output_dir: Path) -> str:
    rel = _relative_href(path, output_dir=output_dir)
    return f'<p><strong>{escape(label)}:</strong> <a href="{escape(rel)}">{escape(rel)}</a></p>'


def _guide_link_for_status(status: str, guide_path: Path, *, output_dir: Path) -> str:
    if status not in {"complete_non_empty", "confirmed_empty"}:
        return (
            "<p><strong>guide panel:</strong> hidden until this field is completed; "
            "candidate guide panels are post-label QC only.</p>"
        )
    return _path_link("guide panel", guide_path, output_dir=output_dir)


def _meta_row(label: str, value: str) -> str:
    return f"<p><strong>{escape(label)}:</strong> {escape(str(value))}</p>"


def _commit_command(package_dir: Path, image_id: str) -> str:
    return "\n".join(
        [
            ".venv/bin/python scripts/commit_manual_reference_mask.py \\",
            f"  --package {shlex.quote(str(package_dir))} \\",
            f"  --image-id {image_id} \\",
            f"  --labels path/to/edited_{image_id}_manual_reference_labels.tif \\",
            "  --labeler YOUR_INITIALS \\",
            '  --notes "brief annotation note"',
        ]
    )


def _confirm_empty_command(package_dir: Path, image_id: str) -> str:
    return "\n".join(
        [
            ".venv/bin/python scripts/commit_manual_reference_mask.py \\",
            f"  --package {shlex.quote(str(package_dir))} \\",
            f"  --image-id {image_id} \\",
            f"  --labels path/to/empty_{image_id}_manual_reference_labels.tif \\",
            "  --status confirmed_empty \\",
            "  --labeler YOUR_INITIALS \\",
            '  --notes "reviewed; no traceable aSMA-associated region"',
        ]
    )


def _validate_manifest_status_match(
    manifest_rows: list[dict[str, str]],
    status_rows: list[dict[str, str]],
) -> None:
    manifest_ids = [_image_id(row) for row in manifest_rows]
    status_ids = [_image_id(row) for row in status_rows]
    missing = sorted(set(manifest_ids) - set(status_ids))
    extra = sorted(set(status_ids) - set(manifest_ids))
    if missing:
        raise ValueError("missing status rows for manifest images: " + ", ".join(missing))
    if extra:
        raise ValueError("status rows not present in manifest: " + ", ".join(extra))
    _require_unique("manual_validation_manifest.csv", manifest_ids)
    _require_unique("manual_labeling_status.csv", status_ids)


def _validate_manifest_audit_match(
    manifest_rows: list[dict[str, str]],
    audit_rows: list[dict[str, str]],
) -> None:
    manifest_ids = [_image_id(row) for row in manifest_rows]
    audit_ids = [_image_id(row) for row in audit_rows]
    missing = sorted(set(manifest_ids) - set(audit_ids))
    extra = sorted(set(audit_ids) - set(manifest_ids))
    if missing:
        raise ValueError("missing audit rows for manifest images: " + ", ".join(missing))
    if extra:
        raise ValueError("audit rows not present in manifest: " + ", ".join(extra))
    _require_unique("manual_annotation_audit.csv", audit_ids)


def _validate_raw_annotation_panel_path(image_id: str, path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{image_id} annotation panel does not exist: {path}")
    if "annotation_panels_raw_only" not in path.parts:
        raise ValueError(
            f"{image_id} annotation_panel_path must point under annotation_panels_raw_only: {path}"
        )
    if "guide_panels" in path.parts:
        raise ValueError(f"{image_id} annotation_panel_path must not point to guide_panels: {path}")


def _require_unique(source: str, image_ids: list[str]) -> None:
    duplicates = sorted({image_id for image_id in image_ids if image_ids.count(image_id) > 1})
    if duplicates:
        raise ValueError(f"duplicate image_id rows in {source}: " + ", ".join(duplicates))


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"required CSV does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _resolve_path(value: str, *, package_dir: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    candidates = [path, package_dir / path, package_dir.parent / path, Path.cwd() / path]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return package_dir / path


def _relative_href(path: Path, *, output_dir: Path) -> str:
    return os.path.relpath(path.resolve(), start=output_dir.resolve())


def _image_id(row: dict[str, str]) -> str:
    return row.get("image_id", "").strip().upper().replace(" ", "")


def _css() -> str:
    return """
body { font-family: Arial, sans-serif; margin: 0; color: #202124; background: #f5f6f7; }
header { padding: 24px 28px; background: #1f2933; color: white; }
h1 { margin: 0 0 8px; font-size: 30px; }
.warning { max-width: 980px; line-height: 1.45; }
.counts span { display: inline-block; margin: 6px 8px 0 0; padding: 4px 8px; border: 1px solid #9fb3c8; border-radius: 4px; }
main { padding: 20px; }
.card { background: white; border: 1px solid #d9e2ec; border-radius: 8px; padding: 16px; margin: 0 0 18px; }
.card-head { display: flex; align-items: center; justify-content: space-between; gap: 16px; border-bottom: 1px solid #e4e7eb; padding-bottom: 10px; }
h2 { margin: 0; font-size: 22px; }
.status { font-weight: bold; padding: 4px 8px; border-radius: 4px; background: #e4e7eb; }
.status-not_started .status { background: #ffe8cc; }
.status-complete_non_empty .status, .status-confirmed_empty .status { background: #d3f9d8; }
.meta { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 6px 14px; margin: 12px 0; }
.meta p, .paths p { margin: 3px 0; }
.panels { display: grid; grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); gap: 14px; }
figure { margin: 0; }
img { max-width: 100%; border: 1px solid #cbd2d9; background: #111; }
figcaption { font-size: 13px; color: #52606d; margin-top: 4px; }
pre { white-space: pre-wrap; background: #111827; color: #f9fafb; padding: 12px; border-radius: 6px; overflow-x: auto; }
.note { color: #7c2d12; font-weight: bold; }
""".strip()
