from __future__ import annotations

import csv
import os
import re
import shutil
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

import yaml

from dapi_norm.cellpose_cell_regions import run_cellpose_cell_region_batch
from dapi_norm.cellpose_endpoint_figures import render_cellpose_endpoint_figures
from dapi_norm.cellpose_plate_summary import (
    build_cellpose_plate_summary,
    write_cellpose_plate_summary_csv,
    write_cellpose_plate_summary_markdown,
)
from dapi_norm.cellpose_runner import run_nuclei_count_batch
from dapi_norm.pi_simple_summary import ImagePair, PiSummaryRow, find_image_pairs, write_pi_workbook
from dapi_norm.pi_simple_summary import _normalize_channel_id
from dapi_norm.seeded_regions import load_mask_lookup_from_counts_root

PRIMARY_ENDPOINT = "target_integrated_intensity_per_DAPI_positive_nucleus"
PRIMARY_SOURCE_COLUMN = (
    "dapi_anchored_cellpose_target_integrated_background_corrected_per_DAPI_positive_nucleus"
)
USER_SUMMARY_COLUMNS = [
    "plate",
    "location",
    "source_id",
    "target_integrated_intensity",
    "normalization_denominator_count",
    PRIMARY_ENDPOINT,
    "qc_status",
    "qc_flags",
    "dapi_anchored_cellpose_object_count",
    "dapi_anchored_cellpose_masked_area_px",
    "cellpose_mask_path",
    "dapi_nuclei_mask_path",
    "target_channel_id",
    "dapi_channel_id",
    "target_path",
    "dapi_path",
    "ch2_path",
    "ch4_path",
]
SKIP_DIR_NAMES = {
    ".git",
    ".models",
    ".pytest_cache",
    ".ruff_cache",
    ".venv",
    "__pycache__",
    "output",
    "reports",
}


@dataclass(frozen=True)
class AcquisitionRun:
    input_root: Path
    plate_name: str
    plate_token: str
    run_name: str
    display_name: str
    image_count: int

    @property
    def output_key(self) -> str:
        return f"{self.plate_token}/{self.run_name}"

    @property
    def output_key_tuple(self) -> tuple[str, str]:
        return (self.plate_token, self.run_name)


@dataclass(frozen=True)
class BatchOutputs:
    output_root: Path
    counts_root: Path
    regions_root: Path
    final_root: Path
    full_summary_csv: Path
    full_summary_markdown: Path
    user_summary_csv: Path
    workbook_path: Path
    run_summary_markdown: Path
    run_summary_html: Path
    run_log_yaml: Path
    figures_root: Path | None


@dataclass(frozen=True)
class BatchRunResult:
    acquisitions: list[AcquisitionRun]
    rows_processed: int
    qc_status_counts: dict[str, int]
    outputs: BatchOutputs


def discover_acquisitions(
    input_root: Path | str,
    *,
    target_channel_id: str = "CH2",
    dapi_channel_id: str = "CH4",
) -> list[AcquisitionRun]:
    root = Path(input_root).expanduser().resolve()
    target_channel = _normalize_channel_id(target_channel_id)
    dapi_channel = _normalize_channel_id(dapi_channel_id)
    if not root.exists():
        raise FileNotFoundError(f"Input folder does not exist: {root}")
    if not root.is_dir():
        raise NotADirectoryError(f"Input path is not a folder: {root}")

    all_pairs = _find_image_pairs_without_output_dirs(
        root,
        target_channel_id=target_channel,
        dapi_channel_id=dapi_channel,
    )
    if not all_pairs:
        raise ValueError(
            f"No {target_channel}/{dapi_channel} TIFF image pairs were found. Expected folders like "
            f"Plate 1/AcquisitionName/XY01/*_{target_channel}.tif and *_{dapi_channel}.tif."
        )

    acquisition_roots = sorted(
        {pair.ch2_path.parent.parent.resolve() for pair in all_pairs},
        key=lambda path: _natural_path_sort_key(path.relative_to(root) if path != root else Path(".")),
    )
    raw_runs = [
        _build_acquisition(
            root,
            acquisition_root,
            target_channel_id=target_channel,
            dapi_channel_id=dapi_channel,
        )
        for acquisition_root in acquisition_roots
    ]
    return _deduplicate_run_names(raw_runs)


def run_user_cellpose_batch(
    *,
    input_root: Path | str,
    output_root: Path | str,
    project_root: Path | str | None = None,
    model_name: str = "cpsam_v2",
    gpu: bool = True,
    background_value: float = 0.0,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    diameter: float | None = None,
    target_channel_id: str = "CH2",
    dapi_channel_id: str = "CH4",
    max_images_per_acquisition: int | None = None,
    render_figures: bool = True,
    max_overlay_images: int | None = None,
    overwrite: bool = False,
) -> BatchRunResult:
    project_path = Path(project_root).resolve() if project_root is not None else _project_root()
    input_path = Path(input_root).expanduser().resolve()
    output_path = Path(output_root).expanduser().resolve()
    target_channel = _normalize_channel_id(target_channel_id)
    dapi_channel = _normalize_channel_id(dapi_channel_id)
    if target_channel == dapi_channel:
        raise ValueError("Target and DAPI channels must be different.")
    _prepare_output_root(output_path, overwrite=overwrite)
    _configure_model_cache(project_path, model_name)

    acquisitions = discover_acquisitions(
        input_path,
        target_channel_id=target_channel,
        dapi_channel_id=dapi_channel,
    )
    if max_images_per_acquisition is not None and max_images_per_acquisition <= 0:
        raise ValueError("--max-images-per-acquisition must be positive when provided.")

    counts_root = output_path / "cellpose_counts"
    regions_root = output_path / "cellpose_cell_regions"
    final_root = output_path / "final"
    logs_root = output_path / "logs"
    logs_root.mkdir(parents=True, exist_ok=True)

    processed_runs: list[dict[str, Any]] = []
    for acquisition in acquisitions:
        pairs = _pairs_for_acquisition(
            acquisition.input_root,
            max_images=max_images_per_acquisition,
            target_channel_id=target_channel,
            dapi_channel_id=dapi_channel,
        )
        counts_dir = counts_root / acquisition.plate_token / acquisition.run_name
        regions_dir = regions_root / acquisition.plate_token / acquisition.run_name

        count_rows = run_nuclei_count_batch(
            input_root=acquisition.input_root,
            output_dir=counts_dir,
            channel_id=dapi_channel,
            model_name=model_name,
            gpu=gpu,
            max_images=max_images_per_acquisition,
            channel_identity_confirmed=True,
        )
        region_rows = run_cellpose_cell_region_batch(
            image_pairs=pairs,
            mask_lookup=load_mask_lookup_from_counts_root(counts_dir),
            output_dir=regions_dir,
            model_name=model_name,
            gpu=gpu,
            background_value=background_value,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
            diameter=diameter,
            target_channel_id=target_channel,
            dapi_channel_id=dapi_channel,
            write_internal_qc=True,
        )
        processed_runs.append(
            {
                **asdict(acquisition),
                "input_root": str(acquisition.input_root),
                "counts_output": str(counts_dir),
                "cell_region_output": str(regions_dir),
                "nuclei_images_processed": len(count_rows),
                "cell_region_images_processed": len(region_rows),
            }
        )

    table_root = final_root / "tables"
    workbook_root = final_root / "workbooks"
    table_root.mkdir(parents=True, exist_ok=True)
    workbook_root.mkdir(parents=True, exist_ok=True)

    rows = build_cellpose_plate_summary(regions_root)
    full_summary_csv = write_cellpose_plate_summary_csv(
        table_root / "cellpose_full_plate_endpoint_summary.csv",
        rows,
    )
    full_summary_markdown = write_cellpose_plate_summary_markdown(
        table_root / "cellpose_full_plate_endpoint_summary.md",
        rows,
    )
    user_summary_csv = write_user_summary_csv(
        table_root / "cellpose_user_friendly_endpoint_summary.csv",
        rows,
    )
    workbook_path = write_pi_style_workbook_from_rows(
        workbook_root / "cellpose_background_corrected_pi_style_summary.xlsx",
        rows,
    )

    figures_root: Path | None = None
    if render_figures and rows:
        figures_root = final_root / "figures"
        render_cellpose_endpoint_figures(
            summary_csv=full_summary_csv,
            output_dir=figures_root,
            panel_page_size=12,
            selected_fields_per_plate=14,
            max_overlay_images=max_overlay_images,
        )

    qc_counts = dict(sorted(Counter(str(row.get("qc_status", "")) for row in rows).items()))
    outputs = BatchOutputs(
        output_root=output_path,
        counts_root=counts_root,
        regions_root=regions_root,
        final_root=final_root,
        full_summary_csv=full_summary_csv,
        full_summary_markdown=full_summary_markdown,
        user_summary_csv=user_summary_csv,
        workbook_path=workbook_path,
        run_summary_markdown=final_root / "START_HERE_RUN_SUMMARY.md",
        run_summary_html=final_root / "START_HERE_RUN_SUMMARY.html",
        run_log_yaml=logs_root / "user_cellpose_batch_run.yaml",
        figures_root=figures_root,
    )
    _write_run_summary_files(
        input_root=input_path,
        project_root=project_path,
        model_name=model_name,
        gpu=gpu,
        background_value=background_value,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        diameter=diameter,
        target_channel_id=target_channel,
        dapi_channel_id=dapi_channel,
        max_images_per_acquisition=max_images_per_acquisition,
        acquisitions=acquisitions,
        processed_runs=processed_runs,
        rows=rows,
        qc_status_counts=qc_counts,
        outputs=outputs,
    )
    return BatchRunResult(
        acquisitions=acquisitions,
        rows_processed=len(rows),
        qc_status_counts=qc_counts,
        outputs=outputs,
    )


def write_user_summary_csv(path: Path | str, rows: Iterable[dict[str, Any]]) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=USER_SUMMARY_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(_user_summary_row(row))
    return output_path


def write_pi_style_workbook_from_rows(path: Path | str, rows: Iterable[dict[str, Any]]) -> Path:
    output_path = Path(path)
    row_list = list(rows)
    grouped: dict[str, list[PiSummaryRow]] = {}
    for row in row_list:
        plate = str(row["plate"])
        grouped.setdefault(plate, []).append(
            PiSummaryRow(
                location=str(row["source_id"]),
                asma_intensity=float(
                    _row_value(
                        row,
                        "dapi_anchored_cellpose_target_integrated_background_corrected",
                        "dapi_anchored_cellpose_ch2_integrated_background_corrected",
                    )
                ),
                nuclei_count=int(row["dapi_positive_nucleus_count"]),
                source_id=str(row["source_id"]),
            )
        )
    write_pi_workbook(output_path, grouped, metadata=_workbook_metadata(row_list))
    return output_path


def _workbook_metadata(rows: list[dict[str, Any]]) -> dict[str, str] | None:
    if not rows:
        return None
    first = rows[0]
    return {
        "target_channel_id": str(first.get("target_channel_id", "CH2")),
        "target_channel_name": "aSMA",
        "dapi_channel_id": str(first.get("dapi_channel_id", "CH4")),
        "dapi_channel_name": "DAPI",
        "primary_endpoint": PRIMARY_ENDPOINT,
        "normalization_denominator": (
            "filtered count of DAPI-positive nuclei; DAPI fluorescence brightness is not used"
        ),
        "target_path_column": "target_path",
        "dapi_path_column": "dapi_path",
    }


def _find_image_pairs_without_output_dirs(
    root: Path,
    *,
    target_channel_id: str = "CH2",
    dapi_channel_id: str = "CH4",
) -> list[ImagePair]:
    return [
        pair
        for pair in find_image_pairs(
            root,
            target_channel=target_channel_id,
            dapi_channel=dapi_channel_id,
        )
        if not _pair_is_in_skipped_dir(pair, root)
    ]


def _build_acquisition(
    input_root: Path,
    acquisition_root: Path,
    *,
    target_channel_id: str = "CH2",
    dapi_channel_id: str = "CH4",
) -> AcquisitionRun:
    pairs = _find_image_pairs_without_output_dirs(
        acquisition_root,
        target_channel_id=target_channel_id,
        dapi_channel_id=dapi_channel_id,
    )
    plate_name = _plate_name_for_path(input_root, acquisition_root)
    run_name = _safe_name(acquisition_root.name)
    try:
        relative_name = acquisition_root.relative_to(input_root)
        display_name = acquisition_root.name if relative_name == Path(".") else str(relative_name)
    except ValueError:
        display_name = acquisition_root.name
    return AcquisitionRun(
        input_root=acquisition_root,
        plate_name=plate_name,
        plate_token=plate_name.replace(" ", "_"),
        run_name=run_name,
        display_name=display_name,
        image_count=len(pairs),
    )


def _deduplicate_run_names(runs: list[AcquisitionRun]) -> list[AcquisitionRun]:
    counts: dict[tuple[str, str], int] = {}
    deduplicated: list[AcquisitionRun] = []
    for run in runs:
        key = (run.plate_token, run.run_name)
        counts[key] = counts.get(key, 0) + 1
        if counts[key] == 1:
            deduplicated.append(run)
            continue
        suffix = _safe_name(run.input_root.parent.name)
        deduplicated.append(
            AcquisitionRun(
                input_root=run.input_root,
                plate_name=run.plate_name,
                plate_token=run.plate_token,
                run_name=f"{run.run_name}_{suffix}_{counts[key]}",
                display_name=run.display_name,
                image_count=run.image_count,
            )
        )
    return deduplicated


def _plate_name_for_path(input_root: Path, acquisition_root: Path) -> str:
    parts = [input_root.name]
    try:
        parts.extend(acquisition_root.relative_to(input_root).parts)
    except ValueError:
        parts.extend(acquisition_root.parts)
    for part in parts:
        match = re.fullmatch(r"plate[\s_-]*(\d+)", part.strip(), flags=re.IGNORECASE)
        if match:
            return f"Plate {int(match.group(1))}"
    return "Plate 1"


def _pairs_for_acquisition(
    acquisition_root: Path,
    *,
    max_images: int | None,
    target_channel_id: str = "CH2",
    dapi_channel_id: str = "CH4",
) -> list[ImagePair]:
    target_channel = _normalize_channel_id(target_channel_id)
    dapi_channel = _normalize_channel_id(dapi_channel_id)
    pairs = _find_image_pairs_without_output_dirs(
        acquisition_root,
        target_channel_id=target_channel,
        dapi_channel_id=dapi_channel,
    )
    if max_images is not None:
        pairs = pairs[:max_images]
    if not pairs:
        raise ValueError(f"No {target_channel}/{dapi_channel} image pairs found under {acquisition_root}")
    duplicates = [location for location, count in Counter(pair.location for pair in pairs).items() if count > 1]
    if duplicates:
        raise ValueError(
            f"Duplicate XY locations inside one acquisition folder {acquisition_root}: "
            + ", ".join(sorted(duplicates))
        )
    return pairs


def _prepare_output_root(output_root: Path, *, overwrite: bool) -> None:
    if output_root.exists() and any(output_root.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Output folder already exists and is not empty: {output_root}. "
            "Choose a new output folder or pass --overwrite."
        )
    output_root.mkdir(parents=True, exist_ok=True)
    if overwrite:
        for managed_name in ["cellpose_counts", "cellpose_cell_regions", "final", "logs"]:
            managed_path = output_root / managed_name
            if managed_path.is_symlink() or managed_path.is_file():
                managed_path.unlink()
            elif managed_path.is_dir():
                shutil.rmtree(managed_path)


def _configure_model_cache(project_root: Path, model_name: str) -> None:
    model_root = project_root / ".models" / "cellpose"
    os.environ.setdefault("CELLPOSE_LOCAL_MODELS_PATH", str(model_root))
    model_path = model_root / model_name
    if not model_path.exists() and not Path(model_name).exists():
        raise FileNotFoundError(
            f"Cellpose model is not available at {model_path}. "
            "Run the project setup once or copy the prepared .models/cellpose folder."
        )


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _user_summary_row(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "plate": row["plate"],
        "location": row["location"],
        "source_id": row["source_id"],
        "target_integrated_intensity": _row_value(
            row,
            "dapi_anchored_cellpose_target_integrated_background_corrected",
            "dapi_anchored_cellpose_ch2_integrated_background_corrected",
        ),
        "normalization_denominator_count": row["dapi_positive_nucleus_count"],
        PRIMARY_ENDPOINT: _row_value(
            row,
            PRIMARY_SOURCE_COLUMN,
            "dapi_anchored_cellpose_ch2_integrated_background_corrected_per_DAPI_positive_nucleus",
        ),
        "qc_status": row["qc_status"],
        "qc_flags": row["qc_flags"],
        "dapi_anchored_cellpose_object_count": _anchored_object_count(row),
        "dapi_anchored_cellpose_masked_area_px": row["dapi_anchored_cellpose_masked_area_px"],
        "cellpose_mask_path": row["cellpose_mask_path"],
        "dapi_nuclei_mask_path": row["dapi_nuclei_mask_path"],
        "target_channel_id": row.get("target_channel_id", "CH2"),
        "dapi_channel_id": row.get("dapi_channel_id", "CH4"),
        "target_path": row.get("target_path", row["ch2_path"]),
        "dapi_path": row.get("dapi_path", row["ch4_path"]),
        "ch2_path": row["ch2_path"],
        "ch4_path": row["ch4_path"],
    }


def _row_value(row: dict[str, Any], preferred: str, fallback: str) -> Any:
    if preferred in row:
        return row[preferred]
    return row[fallback]


def _pair_is_in_skipped_dir(pair: ImagePair, root: Path) -> bool:
    try:
        relative_parts = pair.ch2_path.relative_to(root).parts
    except ValueError:
        relative_parts = pair.ch2_path.parts
    return any(part in SKIP_DIR_NAMES for part in relative_parts)


def _anchored_object_count(row: dict[str, Any]) -> int:
    if "dapi_anchored_cellpose_object_count" in row:
        return _int_value(row["dapi_anchored_cellpose_object_count"])
    return _int_value(row["cellpose_object_count"]) - _int_value(
        row.get("no_dapi_cellpose_object_count_excluded_in_anchored_variant", 0)
    )


def _int_value(value: Any) -> int:
    return int(float(value))


def _write_run_summary_files(
    *,
    input_root: Path,
    project_root: Path,
    model_name: str,
    gpu: bool,
    background_value: float,
    flow_threshold: float,
    cellprob_threshold: float,
    diameter: float | None,
    target_channel_id: str,
    dapi_channel_id: str,
    max_images_per_acquisition: int | None,
    acquisitions: list[AcquisitionRun],
    processed_runs: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    qc_status_counts: dict[str, int],
    outputs: BatchOutputs,
) -> None:
    generated_at = datetime.now(UTC).isoformat()
    outputs.final_root.mkdir(parents=True, exist_ok=True)
    run_log = {
        "generated_at_utc": generated_at,
        "input_root": str(input_root),
        "project_root": str(project_root),
        "output_root": str(outputs.output_root),
        "method": {
            "name": "user_cellpose_dapi_ch2_batch",
            "target_channel_id": target_channel_id,
            "dapi_channel_id": dapi_channel_id,
            "dapi_nuclei_backend": "Cellpose",
            "cell_region_backend": "Cellpose",
            "primary_endpoint": PRIMARY_ENDPOINT,
            "normalization_denominator": (
                "filtered count of DAPI-positive nuclei labels; DAPI fluorescence brightness is "
                "not used as the denominator"
            ),
            "cellpose_region_policy": (
                "Cellpose CH2+CH4 objects retained only when at least one DAPI-positive nucleus "
                "centroid falls inside the object"
            ),
            "whole_cell_segmentation_validated": False,
        },
        "parameters": {
            "model_name": model_name,
            "gpu": gpu,
            "background_value": background_value,
            "flow_threshold": flow_threshold,
            "cellprob_threshold": cellprob_threshold,
            "diameter": diameter,
            "target_channel_id": target_channel_id,
            "dapi_channel_id": dapi_channel_id,
            "max_images_per_acquisition": max_images_per_acquisition,
        },
        "acquisitions": [_acquisition_log_record(acquisition) for acquisition in acquisitions],
        "processed_runs": processed_runs,
        "rows_processed": len(rows),
        "qc_status_counts": qc_status_counts,
        "outputs": {
            "full_summary_csv": str(outputs.full_summary_csv),
            "user_summary_csv": str(outputs.user_summary_csv),
            "workbook": str(outputs.workbook_path),
            "figures_root": "" if outputs.figures_root is None else str(outputs.figures_root),
            "counts_root": str(outputs.counts_root),
            "regions_root": str(outputs.regions_root),
        },
        "validation_status": {
            "automated_smoke_validated": True,
            "manual_ground_truth_metrics_available": False,
            "visual_qc_required": True,
        },
        "python": sys.version.split()[0],
    }
    outputs.run_log_yaml.write_text(yaml.safe_dump(run_log, sort_keys=False), encoding="utf-8")

    md_text = _run_summary_markdown(
        generated_at=generated_at,
        input_root=input_root,
        acquisitions=acquisitions,
        processed_runs=processed_runs,
        rows=rows,
        qc_status_counts=qc_status_counts,
        outputs=outputs,
        target_channel_id=target_channel_id,
        dapi_channel_id=dapi_channel_id,
    )
    outputs.run_summary_markdown.write_text(md_text, encoding="utf-8")
    outputs.run_summary_html.write_text(_markdown_to_simple_html(md_text), encoding="utf-8")


def _run_summary_markdown(
    *,
    generated_at: str,
    input_root: Path,
    acquisitions: list[AcquisitionRun],
    processed_runs: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    qc_status_counts: dict[str, int],
    outputs: BatchOutputs,
    target_channel_id: str,
    dapi_channel_id: str,
) -> str:
    processed_lookup = {
        (str(run["plate_token"]), str(run["run_name"])): int(run["cell_region_images_processed"])
        for run in processed_runs
    }
    acquisition_lines = "\n".join(
        "- "
        + _acquisition_summary_line(acquisition, processed_lookup.get(acquisition.output_key_tuple))
        for acquisition in acquisitions
    )
    qc_lines = "\n".join(f"- {status or 'blank'}: {count}" for status, count in qc_status_counts.items())
    if not qc_lines:
        qc_lines = "- none"
    figures_line = (
        "- QC figures: not rendered for this run"
        if outputs.figures_root is None
        else f"- QC figures: `{outputs.figures_root}`"
    )
    return (
        "# START HERE: Cellpose DAPI/aSMA Batch Run\n\n"
        f"Generated at UTC: `{generated_at}`\n\n"
        "## What Was Measured\n\n"
        f"- Primary endpoint: `{PRIMARY_ENDPOINT}`\n"
        "- DAPI was used to segment and count DAPI-positive nuclei.\n"
        "- The denominator is the filtered DAPI-positive nucleus count, not DAPI fluorescence "
        "brightness. Current nucleus filtering is none.\n"
        "- The Cellpose cell-region output is DAPI-anchored: Cellpose objects are retained only when "
        "at least one DAPI-positive nucleus centroid is inside the object.\n"
        "- These masks are auditable image-analysis regions, not a validated biological interpretation.\n\n"
        "## Inputs\n\n"
        f"- Input folder: `{input_root}`\n"
        f"- Acquisitions processed: {len(acquisitions)}\n"
        f"- Image rows in final table: {len(rows)}\n\n"
        "## Channel Mapping\n\n"
        f"- Target/aSMA measurement channel: `{target_channel_id}`\n"
        f"- DAPI nuclei segmentation/count channel: `{dapi_channel_id}`\n\n"
        f"{acquisition_lines}\n\n"
        "## Main Outputs\n\n"
        f"- User-friendly CSV: `{outputs.user_summary_csv}`\n"
        f"- Full audit CSV: `{outputs.full_summary_csv}`\n"
        f"- PI-style workbook: `{outputs.workbook_path}`\n"
        f"{figures_line}\n"
        f"- Nuclei masks and QC: `{outputs.counts_root}`\n"
        f"- Cellpose region masks and QC: `{outputs.regions_root}`\n\n"
        "## Automated QC Status Counts\n\n"
        f"{qc_lines}\n\n"
        "Visual QC overlays should be inspected before treating results as reviewable. Precision, recall, "
        "F1, and IoU require separate manual ground-truth annotations.\n"
    )


def _acquisition_summary_line(acquisition: AcquisitionRun, processed_count: int | None) -> str:
    if processed_count is None or processed_count == acquisition.image_count:
        count_text = f"{acquisition.image_count} image pairs"
    else:
        count_text = f"{processed_count} of {acquisition.image_count} discovered image pairs"
    return f"{acquisition.plate_name} / {acquisition.display_name}: {count_text}"


def _acquisition_log_record(acquisition: AcquisitionRun) -> dict[str, Any]:
    return {
        "input_root": str(acquisition.input_root),
        "plate_name": acquisition.plate_name,
        "plate_token": acquisition.plate_token,
        "run_name": acquisition.run_name,
        "display_name": acquisition.display_name,
        "image_count": acquisition.image_count,
    }


def _markdown_to_simple_html(markdown_text: str) -> str:
    lines = markdown_text.splitlines()
    html_lines = [
        "<!doctype html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        "<title>Cellpose Batch Run Summary</title>",
        "<style>",
        "body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;line-height:1.45;max-width:960px;margin:32px auto;padding:0 20px;color:#172033}",
        "code{background:#f2f4f8;padding:2px 4px;border-radius:4px}",
        "li{margin:4px 0}",
        "</style>",
        "</head><body>",
    ]
    in_list = False
    for line in lines:
        if line.startswith("# "):
            if in_list:
                html_lines.append("</ul>")
                in_list = False
            html_lines.append(f"<h1>{_html_escape(line[2:])}</h1>")
        elif line.startswith("## "):
            if in_list:
                html_lines.append("</ul>")
                in_list = False
            html_lines.append(f"<h2>{_html_escape(line[3:])}</h2>")
        elif line.startswith("- "):
            if not in_list:
                html_lines.append("<ul>")
                in_list = True
            html_lines.append(f"<li>{_inline_markdown_to_html(line[2:])}</li>")
        elif line.strip():
            if in_list:
                html_lines.append("</ul>")
                in_list = False
            html_lines.append(f"<p>{_inline_markdown_to_html(line)}</p>")
    if in_list:
        html_lines.append("</ul>")
    html_lines.append("</body></html>")
    return "\n".join(html_lines) + "\n"


def _inline_markdown_to_html(text: str) -> str:
    parts = re.split(r"(`[^`]*`)", text)
    rendered = []
    for part in parts:
        if part.startswith("`") and part.endswith("`"):
            rendered.append(f"<code>{_html_escape(part[1:-1])}</code>")
        else:
            rendered.append(_html_escape(part))
    return "".join(rendered)


def _html_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _safe_name(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return safe.strip("._") or "Run"


def _natural_path_sort_key(path: Path) -> tuple[Any, ...]:
    return tuple(_natural_text_sort_key(part) for part in path.parts)


def _natural_text_sort_key(text: str) -> tuple[Any, ...]:
    return tuple(int(token) if token.isdigit() else token.lower() for token in re.split(r"(\d+)", text))
