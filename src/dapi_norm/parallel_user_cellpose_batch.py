"""Bounded concurrent orchestration for isolated Cellpose acquisition runs.

Each acquisition is executed through the established single-acquisition CLI so
its workbook, QC figures, masks, and run log remain self-contained.  This
module only schedules those independent runs and records their provenance; it
does not change the scientific measurement workflow.
"""

from __future__ import annotations

import json
import subprocess
import time
from collections import deque
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, TextIO

from dapi_norm.user_cellpose_batch import discover_acquisitions


STANDARD_WORKBOOK_RELATIVE_PATH = Path(
    "final/workbooks/cellpose_background_corrected_pi_style_summary.xlsx"
)
ORCHESTRATION_DIRECTORY_NAME = "_orchestration"


@dataclass(frozen=True)
class PerAcquisitionRun:
    """An immutable invocation of the existing validated single-run CLI."""

    input_root: Path
    output_root: Path
    output_key: str
    display_name: str
    image_count: int
    command: tuple[str, ...]
    log_path: Path


@dataclass
class PerAcquisitionRunRecord:
    """Persistent status for one scheduled acquisition."""

    run: PerAcquisitionRun
    status: str = "pending"
    return_code: int | None = None
    started_at_utc: str | None = None
    completed_at_utc: str | None = None
    duration_seconds: float | None = None
    error: str | None = None
    _started_monotonic: float | None = None


@dataclass(frozen=True)
class ParallelBatchResult:
    """Final outcome of a bounded parallel orchestration run."""

    output_root: Path
    manifest_path: Path
    summary_path: Path
    records: tuple[PerAcquisitionRunRecord, ...]

    @property
    def succeeded(self) -> bool:
        return all(record.status == "succeeded" for record in self.records)


def build_per_acquisition_runs(
    *,
    input_root: Path | str,
    output_root: Path | str,
    python_executable: Path | str,
    runner_script: Path | str,
    target_channel_id: str = "CH2",
    dapi_channel_id: str = "CH4",
    model_name: str = "cpsam_v2",
    gpu: bool = True,
    background_value: float = 0.0,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    diameter: float | None = None,
    max_images_per_acquisition: int | None = None,
    render_figures: bool = True,
    max_overlay_images: int | None = None,
) -> list[PerAcquisitionRun]:
    """Discover input acquisitions and construct one non-overlapping run each.

    This is intentionally a planning function: it neither creates output
    folders nor starts Cellpose.  A caller can use it for a safe dry run.
    """

    input_path = Path(input_root).expanduser().resolve()
    output_path = Path(output_root).expanduser().resolve()
    python_path = Path(python_executable).expanduser().resolve()
    runner_path = Path(runner_script).expanduser().resolve()
    acquisitions = discover_acquisitions(
        input_path,
        target_channel_id=target_channel_id,
        dapi_channel_id=dapi_channel_id,
    )
    if not acquisitions:
        raise ValueError(
            "No valid acquisition folders were found. Expected XY folders containing "
            f"paired {target_channel_id} and {dapi_channel_id} TIFF files."
        )

    logs_root = output_path / ORCHESTRATION_DIRECTORY_NAME / "logs"
    runs: list[PerAcquisitionRun] = []
    output_roots: set[Path] = set()
    for index, acquisition in enumerate(acquisitions, start=1):
        acquisition_output = output_path / Path(acquisition.output_key)
        if acquisition_output in output_roots:
            raise ValueError(
                "Two discovered acquisitions would write to the same output directory: "
                f"{acquisition_output}"
            )
        output_roots.add(acquisition_output)
        log_stem = acquisition.output_key.replace("/", "__")
        command = _single_acquisition_command(
            python_path=python_path,
            runner_path=runner_path,
            input_root=acquisition.input_root,
            output_root=acquisition_output,
            target_channel_id=target_channel_id,
            dapi_channel_id=dapi_channel_id,
            model_name=model_name,
            gpu=gpu,
            background_value=background_value,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
            diameter=diameter,
            max_images_per_acquisition=max_images_per_acquisition,
            render_figures=render_figures,
            max_overlay_images=max_overlay_images,
        )
        runs.append(
            PerAcquisitionRun(
                input_root=acquisition.input_root,
                output_root=acquisition_output,
                output_key=acquisition.output_key,
                display_name=acquisition.display_name,
                image_count=acquisition.image_count,
                command=tuple(command),
                log_path=logs_root / f"{index:03d}__{log_stem}.log",
            )
        )
    return runs


def validate_new_output_root(output_root: Path | str, runs: list[PerAcquisitionRun]) -> None:
    """Refuse to mix a new benchmark/production run with existing outputs."""

    output_path = Path(output_root).expanduser().resolve()
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(
            f"Output root must be new or empty to prevent mixed results: {output_path}"
        )
    for run in runs:
        if run.output_root.exists() and any(run.output_root.iterdir()):
            raise FileExistsError(
                "An acquisition output folder already contains files; choose a new output root: "
                f"{run.output_root}"
            )


def execute_parallel_runs(
    *,
    runs: list[PerAcquisitionRun],
    output_root: Path | str,
    project_root: Path | str,
    workers: int,
    continue_on_error: bool = False,
    poll_seconds: float = 0.2,
) -> ParallelBatchResult:
    """Run isolated acquisitions with bounded process concurrency.

    Each child process receives its own input and output root.  By default an
    observed failure stops scheduling later runs but lets already-running work
    finish, preserving every completed bundle for audit.
    """

    if workers < 1:
        raise ValueError("workers must be at least 1.")
    if not runs:
        raise ValueError("At least one acquisition run is required.")
    if poll_seconds <= 0:
        raise ValueError("poll_seconds must be positive.")

    output_path = Path(output_root).expanduser().resolve()
    project_path = Path(project_root).expanduser().resolve()
    validate_new_output_root(output_path, runs)
    logs_root = output_path / ORCHESTRATION_DIRECTORY_NAME / "logs"
    logs_root.mkdir(parents=True, exist_ok=False)
    manifest_path = output_path / ORCHESTRATION_DIRECTORY_NAME / "parallel_run_manifest.json"
    summary_path = output_path / "START_HERE_PARALLEL_RUN.md"

    records = [PerAcquisitionRunRecord(run=run) for run in runs]
    records_by_output = {record.run.output_root: record for record in records}
    pending: deque[PerAcquisitionRun] = deque(runs)
    active: dict[Path, tuple[subprocess.Popen[str], TextIO]] = {}
    scheduling_halted = False
    _write_run_artifacts(
        manifest_path=manifest_path,
        summary_path=summary_path,
        workers=workers,
        records=records,
        status="running",
    )

    try:
        while pending or active:
            while pending and len(active) < workers and not scheduling_halted:
                run = pending.popleft()
                record = records_by_output[run.output_root]
                record.started_at_utc = _utc_now()
                record._started_monotonic = time.monotonic()
                record.status = "running"
                log_handle = run.log_path.open("w", encoding="utf-8", newline="\n")
                log_handle.write("Command:\n")
                log_handle.write(" ".join(_quote_command_part(part) for part in run.command))
                log_handle.write("\n\nOutput:\n")
                log_handle.flush()
                try:
                    process = subprocess.Popen(
                        list(run.command),
                        cwd=str(project_path),
                        stdout=log_handle,
                        stderr=subprocess.STDOUT,
                        text=True,
                    )
                except OSError as exc:
                    log_handle.close()
                    record.status = "failed_to_start"
                    record.error = str(exc)
                    record.completed_at_utc = _utc_now()
                    record.duration_seconds = 0.0
                    if not continue_on_error:
                        scheduling_halted = True
                else:
                    active[run.output_root] = (process, log_handle)
                _write_run_artifacts(
                    manifest_path=manifest_path,
                    summary_path=summary_path,
                    workers=workers,
                    records=records,
                    status="running",
                )

            if active:
                time.sleep(poll_seconds)
                for output_key, (process, log_handle) in list(active.items()):
                    return_code = process.poll()
                    if return_code is None:
                        continue
                    log_handle.close()
                    del active[output_key]
                    record = records_by_output[output_key]
                    record.return_code = return_code
                    record.completed_at_utc = _utc_now()
                    record.duration_seconds = _duration_since(record._started_monotonic)
                    record._started_monotonic = None
                    record.status = "succeeded" if return_code == 0 else "failed"
                    if return_code != 0 and not continue_on_error:
                        scheduling_halted = True
                    _write_run_artifacts(
                        manifest_path=manifest_path,
                        summary_path=summary_path,
                        workers=workers,
                        records=records,
                        status="running",
                    )

            if scheduling_halted and pending:
                while pending:
                    run = pending.popleft()
                    record = records_by_output[run.output_root]
                    record.status = "not_started_due_to_earlier_failure"
                _write_run_artifacts(
                    manifest_path=manifest_path,
                    summary_path=summary_path,
                    workers=workers,
                    records=records,
                    status="running",
                )
    except KeyboardInterrupt:
        for process, log_handle in active.values():
            process.terminate()
            log_handle.close()
        for record in records:
            if record.status == "running":
                record.status = "interrupted"
                record.completed_at_utc = _utc_now()
                record.duration_seconds = _duration_since(record._started_monotonic)
                record._started_monotonic = None
        _write_run_artifacts(
            manifest_path=manifest_path,
            summary_path=summary_path,
            workers=workers,
            records=records,
            status="interrupted",
        )
        raise

    final_status = "completed" if all(record.status == "succeeded" for record in records) else "completed_with_errors"
    _write_run_artifacts(
        manifest_path=manifest_path,
        summary_path=summary_path,
        workers=workers,
        records=records,
        status=final_status,
    )
    return ParallelBatchResult(
        output_root=output_path,
        manifest_path=manifest_path,
        summary_path=summary_path,
        records=tuple(records),
    )


def _single_acquisition_command(
    *,
    python_path: Path,
    runner_path: Path,
    input_root: Path,
    output_root: Path,
    target_channel_id: str,
    dapi_channel_id: str,
    model_name: str,
    gpu: bool,
    background_value: float,
    flow_threshold: float,
    cellprob_threshold: float,
    diameter: float | None,
    max_images_per_acquisition: int | None,
    render_figures: bool,
    max_overlay_images: int | None,
) -> list[str]:
    command = [
        str(python_path),
        str(runner_path),
        "--input",
        str(input_root),
        "--output",
        str(output_root),
        "--model",
        model_name,
        "--gpu" if gpu else "--cpu",
        "--background-value",
        str(background_value),
        "--flow-threshold",
        str(flow_threshold),
        "--cellprob-threshold",
        str(cellprob_threshold),
        "--target-channel",
        target_channel_id,
        "--dapi-channel",
        dapi_channel_id,
        "--render-figures" if render_figures else "--skip-figures",
    ]
    if diameter is not None:
        command.extend(["--diameter", str(diameter)])
    if max_images_per_acquisition is not None:
        command.extend(["--max-images-per-acquisition", str(max_images_per_acquisition)])
    if max_overlay_images is not None:
        command.extend(["--max-overlay-images", str(max_overlay_images)])
    return command


def _write_run_artifacts(
    *,
    manifest_path: Path,
    summary_path: Path,
    workers: int,
    records: list[PerAcquisitionRunRecord],
    status: str,
) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema_version": 1,
        "status": status,
        "workers": workers,
        "updated_at_utc": _utc_now(),
        "runs": [_record_to_json(record) for record in records],
    }
    _atomic_write_text(manifest_path, json.dumps(manifest, indent=2) + "\n")
    _atomic_write_text(summary_path, _summary_markdown(status=status, workers=workers, records=records))


def _record_to_json(record: PerAcquisitionRunRecord) -> dict[str, Any]:
    return {
        "display_name": record.run.display_name,
        "input_root": str(record.run.input_root),
        "output_root": str(record.run.output_root),
        "output_key": record.run.output_key,
        "image_count": record.run.image_count,
        "command": list(record.run.command),
        "log_path": str(record.run.log_path),
        "expected_workbook_path": str(record.run.output_root / STANDARD_WORKBOOK_RELATIVE_PATH),
        "status": record.status,
        "return_code": record.return_code,
        "started_at_utc": record.started_at_utc,
        "completed_at_utc": record.completed_at_utc,
        "duration_seconds": record.duration_seconds,
        "error": record.error,
    }


def _summary_markdown(
    *, status: str,
    workers: int,
    records: list[PerAcquisitionRunRecord],
) -> str:
    lines = [
        "# Parallel per-acquisition Cellpose run",
        "",
        f"- Status: `{status}`",
        f"- Concurrent workers: `{workers}`",
        "- Each acquisition has an isolated standard Cellpose output bundle.",
        "",
        "| Acquisition | Image pairs | Status | Standard workbook | Child log |",
        "| --- | ---: | --- | --- | --- |",
    ]
    for record in records:
        workbook = record.run.output_root / STANDARD_WORKBOOK_RELATIVE_PATH
        lines.append(
            "| {name} | {count} | {status} | `{workbook}` | `{log}` |".format(
                name=record.run.display_name.replace("|", "\\|"),
                count=record.run.image_count,
                status=record.status,
                workbook=workbook,
                log=record.run.log_path,
            )
        )
    return "\n".join(lines) + "\n"


def _atomic_write_text(path: Path, content: str) -> None:
    temporary_path = path.with_suffix(path.suffix + ".tmp")
    temporary_path.write_text(content, encoding="utf-8", newline="\n")
    temporary_path.replace(path)


def _duration_since(started_monotonic: float | None) -> float | None:
    if started_monotonic is None:
        return None
    return round(time.monotonic() - started_monotonic, 3)


def _quote_command_part(part: str) -> str:
    return f'"{part}"' if " " in part else part


def _utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
