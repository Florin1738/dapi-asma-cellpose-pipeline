from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_gui():
    spec = importlib.util.spec_from_file_location(
        "cellpose_gui", REPO_ROOT / "scripts" / "cellpose_gui.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gui = _load_gui()


@pytest.fixture()
def app():
    tk = pytest.importorskip("tkinter")
    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("No display available for Tk.")
    root.withdraw()
    # Silence dialogs and record calls.
    calls: list[tuple[str, tuple]] = []
    gui.messagebox.showerror = lambda *a, **k: calls.append(("error", a))
    gui.messagebox.showwarning = lambda *a, **k: calls.append(("warning", a))
    gui.messagebox.showinfo = lambda *a, **k: calls.append(("info", a))
    instance = gui.PipelineGUI(root, REPO_ROOT)
    instance._toggle_advanced()  # build advanced widgets
    instance._dialog_calls = calls  # type: ignore[attr-defined]
    yield instance
    root.destroy()


DATA = REPO_ROOT / "data" / "aSMA_DAPI_plates"


@pytest.mark.skipif(not DATA.is_dir(), reason="Bundled plate data not present.")
def test_valid_inputs_build_a_forwarding_command(app):
    app.input_var.set(str(DATA))
    app.output_var.set(str(REPO_ROOT))
    cmd = app._build_command()
    assert cmd is not None
    # Pure flag-forwarding: the scientific contract is unchanged.
    assert "--input" in cmd and "--output" in cmd
    assert "--gpu" in cmd and "--model" in cmd
    assert "--target-channel" in cmd and "--dapi-channel" in cmd
    assert cmd[cmd.index("--target-channel") + 1] == "CH2"
    assert cmd[cmd.index("--dapi-channel") + 1] == "CH4"
    assert "--flow-threshold" in cmd and "--cellprob-threshold" in cmd


@pytest.mark.skipif(not DATA.is_dir(), reason="Bundled plate data not present.")
def test_channel_mapping_fields_forward_selected_channels(app):
    app.input_var.set(str(DATA))
    app.output_var.set(str(REPO_ROOT))
    app.target_channel_var.set("CH1")
    app.dapi_channel_var.set("CH4")

    cmd = app._build_command()

    assert cmd is not None
    assert cmd[cmd.index("--target-channel") + 1] == "CH1"
    assert cmd[cmd.index("--dapi-channel") + 1] == "CH4"


@pytest.mark.skipif(not DATA.is_dir(), reason="Bundled plate data not present.")
@pytest.mark.parametrize("field,value", [("flow_var", "nan"), ("bg_var", "inf"),
                                         ("cellprob_var", "abc")])
def test_bad_numeric_values_are_rejected(app, field, value):
    app.input_var.set(str(DATA))
    app.output_var.set(str(REPO_ROOT))
    getattr(app, field).set(value)
    assert app._build_command() is None


@pytest.mark.skipif(not DATA.is_dir(), reason="Bundled plate data not present.")
def test_max_images_zero_is_rejected(app):
    app.input_var.set(str(DATA))
    app.output_var.set(str(REPO_ROOT))
    app.maximg_var.set("0")
    assert app._build_command() is None


def test_missing_input_folder_is_rejected(app):
    app.input_var.set("/definitely/not/a/real/folder/xyz")
    app.output_var.set(str(REPO_ROOT))
    assert app._build_command() is None


def test_zero_row_run_warns_instead_of_claiming_success(app):
    app.cancelled = False
    app.rows_processed = 0
    app.run_output = REPO_ROOT / "nonexistent_output"
    app._finish_run(0)
    assert any(kind == "warning" for kind, _ in app._dialog_calls)
