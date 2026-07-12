from __future__ import annotations

import ast
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_gui_module_parses_and_targets_the_shared_runner():
    gui = REPO_ROOT / "scripts" / "cellpose_gui.py"
    text = gui.read_text(encoding="utf-8")
    # Valid Python.
    ast.parse(text)
    # Delegates to the exact same runner as the command-line path.
    assert "run_user_cellpose_batch.py" in text
    # Exposes the folder inputs and the advanced knobs behind a toggle.
    assert "Advanced options" in text
    assert "--input" in text and "--output" in text
    assert "--gpu" in text and "--cpu" in text
    assert "--flow-threshold" in text and "--cellprob-threshold" in text
    # Same DAPI-anchored endpoint output contract.
    assert "START_HERE_RUN_SUMMARY.html" in text
    assert "CELLPOSE_LOCAL_MODELS_PATH" in text


def test_macos_run_launcher_self_heals_and_opens_gui():
    text = (REPO_ROOT / "Run Cellpose DAPI aSMA Pipeline.command").read_text(
        encoding="utf-8"
    )
    # Self-heal: runs setup when the environment is missing.
    assert "scripts/setup_cellpose_pipeline_macos.sh" in text
    assert "env_ready" in text
    # Launches the GUI, not the old folder-picker flow.
    assert "scripts/cellpose_gui.py" in text
    assert "--project-dir" in text
    assert 'CELLPOSE_LOCAL_MODELS_PATH="$PROJECT_DIR/.models/cellpose"' in text


def test_windows_gui_launcher_self_heals_and_opens_gui():
    text = (REPO_ROOT / "scripts" / "run_cellpose_gui_windows.ps1").read_text(
        encoding="utf-8"
    )
    # Self-heal: runs setup when the environment is missing.
    assert "setup_cellpose_pipeline_windows.ps1" in text
    assert "Test-EnvReady" in text
    # Launches the GUI, preferring the no-console interpreter.
    assert "scripts\\cellpose_gui.py" in text
    assert "pythonw.exe" in text
    assert "--project-dir" in text
    # Readiness delegated to the shared cross-platform module.
    assert "scripts\\pipeline_env.py" in text
    assert "$env:CELLPOSE_LOCAL_MODELS_PATH" in text
