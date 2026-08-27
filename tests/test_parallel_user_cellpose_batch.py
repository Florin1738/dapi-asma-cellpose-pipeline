from __future__ import annotations

import json
from pathlib import Path

import pytest

from dapi_norm.parallel_user_cellpose_batch import (
    ORCHESTRATION_DIRECTORY_NAME,
    PerAcquisitionRun,
    build_per_acquisition_runs,
    execute_parallel_runs,
    validate_new_output_root,
)


def test_build_per_acquisition_runs_preserves_isolated_output_bundles(tmp_path: Path):
    input_root = tmp_path / "inputs"
    _touch_pair(input_root / "run a" / "XY01", "RunA_XY01")
    _touch_pair(input_root / "run b" / "XY02", "RunB_XY02")
    output_root = tmp_path / "outputs"

    runs = build_per_acquisition_runs(
        input_root=input_root,
        output_root=output_root,
        python_executable=tmp_path / "python.exe",
        runner_script=tmp_path / "run_user_cellpose_batch.py",
        target_channel_id="CH2",
        dapi_channel_id="CH4",
        gpu=True,
        max_images_per_acquisition=10,
    )

    assert [run.output_key for run in runs] == ["Plate_1/run_a", "Plate_1/run_b"]
    assert [run.output_root for run in runs] == [
        (output_root / "Plate_1" / "run_a").resolve(),
        (output_root / "Plate_1" / "run_b").resolve(),
    ]
    assert all("--gpu" in run.command for run in runs)
    assert all("--target-channel" in run.command and "CH2" in run.command for run in runs)
    assert all("--dapi-channel" in run.command and "CH4" in run.command for run in runs)
    assert all("--max-images-per-acquisition" in run.command for run in runs)
    assert runs[0].log_path == (
        output_root.resolve() / ORCHESTRATION_DIRECTORY_NAME / "logs" / "001__Plate_1__run_a.log"
    )


def test_validate_new_output_root_refuses_mixed_results(tmp_path: Path):
    output_root = tmp_path / "outputs"
    output_root.mkdir()
    (output_root / "prior_result.txt").write_text("do not mix", encoding="utf-8")
    run = _run_spec(output_root, "one")

    with pytest.raises(FileExistsError, match="new or empty"):
        validate_new_output_root(output_root, [run])


def test_build_per_acquisition_runs_selects_only_exact_requested_acquisitions(tmp_path: Path):
    input_root = tmp_path / "inputs"
    _touch_pair(input_root / "run a" / "XY01", "RunA_XY01")
    _touch_pair(input_root / "run b" / "XY02", "RunB_XY02")

    runs = build_per_acquisition_runs(
        input_root=input_root,
        output_root=tmp_path / "outputs",
        python_executable=tmp_path / "python.exe",
        runner_script=tmp_path / "run_user_cellpose_batch.py",
        acquisition_selectors=["run b"],
    )

    assert [run.display_name for run in runs] == ["run b"]


def test_build_per_acquisition_runs_rejects_unknown_acquisition_selector(tmp_path: Path):
    input_root = tmp_path / "inputs"
    _touch_pair(input_root / "run a" / "XY01", "RunA_XY01")

    with pytest.raises(ValueError, match="not found"):
        build_per_acquisition_runs(
            input_root=input_root,
            output_root=tmp_path / "outputs",
            python_executable=tmp_path / "python.exe",
            runner_script=tmp_path / "run_user_cellpose_batch.py",
            acquisition_selectors=["run missing"],
        )


def test_execute_parallel_runs_starts_two_children_before_polling(tmp_path: Path, monkeypatch):
    output_root = tmp_path / "outputs"
    runs = [_run_spec(output_root, "one"), _run_spec(output_root, "two")]
    created: list[_SuccessfulFakeProcess] = []
    active_counts_when_sleeping: list[int] = []

    class FakePopen(_SuccessfulFakeProcess):
        def __init__(self, *args, **kwargs):
            super().__init__()
            created.append(self)

    monkeypatch.setattr(
        "dapi_norm.parallel_user_cellpose_batch.subprocess.Popen", FakePopen
    )
    monkeypatch.setattr(
        "dapi_norm.parallel_user_cellpose_batch.time.sleep",
        lambda _: active_counts_when_sleeping.append(len(created)),
    )

    result = execute_parallel_runs(
        runs=runs,
        output_root=output_root,
        project_root=tmp_path,
        workers=2,
        poll_seconds=0.01,
    )

    assert active_counts_when_sleeping[0] == 2
    assert [record.status for record in result.records] == ["succeeded", "succeeded"]
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["status"] == "completed"
    assert manifest["workers"] == 2
    assert [entry["status"] for entry in manifest["runs"]] == ["succeeded", "succeeded"]
    assert all(
        Path(entry["expected_workbook_path"]).name
        == "aSMA Quantification Workbook - input.xlsx"
        for entry in manifest["runs"]
    )
    assert result.summary_path.read_text(encoding="utf-8").count("succeeded") == 2
    assert "aSMA Quantification Workbook - input.xlsx" in result.summary_path.read_text(
        encoding="utf-8"
    )
    assert all(run.log_path.exists() for run in runs)


def test_execute_parallel_runs_stops_scheduling_after_failure_by_default(
    tmp_path: Path, monkeypatch
):
    output_root = tmp_path / "outputs"
    runs = [_run_spec(output_root, "fails"), _run_spec(output_root, "not_started")]

    class FakePopen(_FailedFakeProcess):
        pass

    monkeypatch.setattr(
        "dapi_norm.parallel_user_cellpose_batch.subprocess.Popen", FakePopen
    )
    monkeypatch.setattr("dapi_norm.parallel_user_cellpose_batch.time.sleep", lambda _: None)

    result = execute_parallel_runs(
        runs=runs,
        output_root=output_root,
        project_root=tmp_path,
        workers=1,
        poll_seconds=0.01,
    )

    assert [record.status for record in result.records] == [
        "failed",
        "not_started_due_to_earlier_failure",
    ]
    assert not result.succeeded


class _SuccessfulFakeProcess:
    def __init__(self, *args, **kwargs):
        pass

    def poll(self) -> int:
        return 0


class _FailedFakeProcess:
    def __init__(self, *args, **kwargs):
        pass

    def poll(self) -> int:
        return 3


def _run_spec(output_root: Path, name: str) -> PerAcquisitionRun:
    return PerAcquisitionRun(
        input_root=Path(f"C:/{name}/input"),
        output_root=output_root / "Plate_1" / name,
        output_key=f"Plate_1/{name}",
        display_name=name,
        image_count=1,
        command=("python", "-c", "pass"),
        log_path=output_root / ORCHESTRATION_DIRECTORY_NAME / "logs" / f"{name}.log",
    )


def _touch_pair(xy_dir: Path, stem: str) -> None:
    xy_dir.mkdir(parents=True)
    (xy_dir / f"{stem}_CH2.tif").touch()
    (xy_dir / f"{stem}_CH4.tif").touch()
