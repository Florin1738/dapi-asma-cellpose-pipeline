from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_macos_setup_launcher_calls_project_setup_helper():
    text = (REPO_ROOT / "Setup Cellpose DAPI aSMA Pipeline.command").read_text(
        encoding="utf-8"
    )

    assert "scripts/setup_cellpose_pipeline_macos.sh" in text
    assert "/bin/bash" in text
    assert "Run Cellpose DAPI aSMA Pipeline.command" in text


def test_macos_setup_helper_creates_local_environment_and_delegates_finalization():
    text = (REPO_ROOT / "scripts" / "setup_cellpose_pipeline_macos.sh").read_text(
        encoding="utf-8"
    )

    # OS-specific bootstrap stays in the shell script.
    assert '"$UV_BIN" venv --python "$PYTHON_VERSION" .venv' in text
    assert ".venv/bin/python" in text
    assert ".[cellpose]" in text
    # The substantive, cross-platform steps are delegated to the shared module.
    assert "scripts/pipeline_env.py" in text
    assert "finish-setup" in text


def test_windows_setup_launcher_calls_project_setup_helper():
    text = (REPO_ROOT / "Setup Cellpose DAPI aSMA Pipeline Windows.cmd").read_text(
        encoding="utf-8"
    )

    assert "scripts\\setup_cellpose_pipeline_windows.ps1" in text
    assert "powershell.exe -NoProfile -ExecutionPolicy Bypass" in text
    assert '-ProjectDir "%PROJECT_DIR%"' in text


def test_windows_setup_helper_creates_local_environment_and_delegates_finalization():
    text = (REPO_ROOT / "scripts" / "setup_cellpose_pipeline_windows.ps1").read_text(
        encoding="utf-8"
    )

    # OS-specific bootstrap stays in the PowerShell script.
    assert "& $Uv venv --python $PythonVersion .venv" in text
    assert ".venv\\Scripts\\python.exe" in text
    assert ".[cellpose]" in text
    # The substantive, cross-platform steps are delegated to the shared module.
    assert "scripts\\pipeline_env.py" in text
    assert "finish-setup" in text
