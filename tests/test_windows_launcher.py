from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_windows_cmd_invokes_gui_helper_in_sta_mode():
    text = (REPO_ROOT / "Run Cellpose DAPI aSMA Pipeline Windows.cmd").read_text(
        encoding="utf-8"
    )

    assert "scripts\\run_cellpose_gui_windows.ps1" in text
    assert "powershell.exe -NoProfile -STA -ExecutionPolicy Bypass" in text
    assert '-ProjectDir "%PROJECT_DIR%"' in text


def test_windows_powershell_launcher_uses_prepared_windows_environment():
    text = (REPO_ROOT / "scripts" / "run_user_cellpose_batch_windows.ps1").read_text(
        encoding="utf-8"
    )

    assert '".venv\\Scripts\\python.exe"' in text
    assert '".models\\cellpose\\cpsam_v2"' in text
    assert '"scripts\\run_user_cellpose_batch.py"' in text
    assert "import cellpose; import dapi_norm.user_cellpose_batch" in text
    assert "$env:CELLPOSE_LOCAL_MODELS_PATH" in text
    assert "--input $InputFolder --output $RunOutput" in text
    assert "START_HERE_RUN_SUMMARY.html" in text
    assert "System.Windows.Forms.FolderBrowserDialog" in text
