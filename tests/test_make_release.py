from __future__ import annotations

import importlib.util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location(
        "make_release", REPO_ROOT / "scripts" / "make_release.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mr = _load()


def _namelist(zip_path: Path) -> list[str]:
    import zipfile

    return zipfile.ZipFile(zip_path).namelist()


def test_mac_zip_excludes_private_and_heavy_content(tmp_path):
    out = tmp_path / "mac.zip"
    mr.build_zip(out, "rel-mac", mr.MAC_LAUNCHERS, with_model=False, with_wheelhouse=False)
    names = _namelist(out)

    forbidden = ("/data/", "/output/", "/.venv/", "/.models/", "/reports/", "/manual_validation/")
    leaks = [
        n
        for n in names
        if any(p in n for p in forbidden) or n.endswith((".tif", ".tiff", ".npy"))
    ]
    assert leaks == [], f"private/heavy content leaked into release: {leaks}"


def test_mac_zip_includes_what_a_user_needs(tmp_path):
    out = tmp_path / "mac.zip"
    mr.build_zip(out, "rel-mac", mr.MAC_LAUNCHERS, with_model=False, with_wheelhouse=False)
    names = _namelist(out)

    assert any(n.endswith("READ ME FIRST.txt") for n in names)
    assert any(n.endswith("Run Cellpose DAPI aSMA Pipeline.command") for n in names)
    assert any(n.endswith("scripts/cellpose_gui.py") for n in names)
    assert any(n.endswith("scripts/pipeline_env.py") for n in names)
    # No Windows launchers in the Mac zip.
    assert not any(n.endswith(".cmd") for n in names)


def test_windows_zip_has_windows_launchers_only(tmp_path):
    out = tmp_path / "win.zip"
    mr.build_zip(out, "rel-win", mr.WINDOWS_LAUNCHERS, with_model=False, with_wheelhouse=False)
    names = _namelist(out)

    assert any(n.endswith("Run Cellpose DAPI aSMA Pipeline Windows.cmd") for n in names)
    assert not any(n.endswith(".command") for n in names)
