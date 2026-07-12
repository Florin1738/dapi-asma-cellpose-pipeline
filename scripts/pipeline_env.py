#!/usr/bin/env python
"""Cross-platform environment logic shared by the Mac and Windows launchers.

The OS-specific setup scripts only do what genuinely differs between platforms:
locate/install ``uv``, create the ``.venv``, and ``pip install`` the project.
Everything after that — downloading and checksum-verifying the Cellpose model,
verifying imports, the discovery dry run, and the "is the environment ready?"
check the run launchers use — lives here, in one place, so a test on macOS
exercises the exact same code path Windows runs.

Run under the project's ``.venv`` interpreter:

    python scripts/pipeline_env.py check        --project-dir PROJECT
    python scripts/pipeline_env.py finish-setup --project-dir PROJECT
    python scripts/pipeline_env.py verify-model  --project-dir PROJECT
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import time
import urllib.request
from pathlib import Path

MODEL_URL = "https://huggingface.co/mouseland/cellpose-sam/resolve/main/cpsam_v2"
MODEL_SHA256 = "0f1cc3f7ecdd8a037a57c6c48d9d8921391be4cbce3fa9f13c3e3a2e1253c667"


def model_dir(project_dir: Path) -> Path:
    return project_dir / ".models" / "cellpose"


def model_path(project_dir: Path) -> Path:
    return model_dir(project_dir) / "cpsam_v2"


def venv_python(project_dir: Path) -> Path:
    if sys.platform.startswith("win"):
        return project_dir / ".venv" / "Scripts" / "python.exe"
    return project_dir / ".venv" / "bin" / "python"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_model(project_dir: Path) -> tuple[bool, str]:
    """Return (ok, message). ok is True only if the model exists and matches."""
    path = model_path(project_dir)
    if not path.is_file():
        return False, f"Cellpose model not found at {path}."
    actual = _sha256(path)
    if actual != MODEL_SHA256:
        return False, (
            f"Cellpose model checksum mismatch at {path}. "
            f"Expected {MODEL_SHA256} but found {actual}."
        )
    return True, "Model present and checksum verified."


def download_model(project_dir: Path, attempts: int = 8) -> None:
    """Download cpsam_v2 if missing, then verify the checksum. Idempotent."""
    ok, _ = verify_model(project_dir)
    if ok:
        return
    md = model_dir(project_dir)
    md.mkdir(parents=True, exist_ok=True)
    target = model_path(project_dir)
    tmp = target.with_suffix(".download")
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(MODEL_URL, timeout=60) as response, tmp.open("wb") as out:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            # Verify the downloaded bytes BEFORE promoting to the real path, so a
            # "200 OK but wrong content" response (e.g. an error/rate-limit page)
            # is retried instead of accepted.
            actual = _sha256(tmp)
            if actual != MODEL_SHA256:
                raise RuntimeError(
                    f"Downloaded model checksum mismatch (got {actual}). "
                    "The download was incomplete or corrupted."
                )
            tmp.replace(target)
            return
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if tmp.exists():
                tmp.unlink()
            if attempt >= attempts:
                raise RuntimeError(
                    f"Could not download a valid Cellpose model after {attempts} attempts. "
                    f"Check your internet connection and try setup again. Last error: {exc}"
                ) from last_error
            time.sleep(5)


def verify_imports() -> tuple[bool, str]:
    try:
        import tkinter  # noqa: F401,PLC0415  (the desktop app needs it)
    except Exception as exc:  # noqa: BLE001
        return False, (
            "The graphical app toolkit (tkinter/Tk) is not available in this Python "
            f"build, so the app window cannot open: {exc}"
        )
    try:
        import cellpose  # noqa: PLC0415

        import dapi_norm.user_cellpose_batch  # noqa: F401,PLC0415
    except Exception as exc:  # noqa: BLE001
        return False, f"Cellpose or the project package could not be imported: {exc}"
    version = getattr(cellpose, "__version__", None)
    if not version:
        try:
            from importlib.metadata import version as _pkg_version  # noqa: PLC0415

            version = _pkg_version("cellpose")
        except Exception:  # noqa: BLE001
            version = "unknown"
    return True, f"cellpose {version}; project package import ok."


def env_ready(project_dir: Path) -> tuple[bool, str]:
    """Everything a run launcher needs before opening the app."""
    if not venv_python(project_dir).is_file():
        return False, "The analysis environment (.venv) is missing. Run setup first."
    ok, message = verify_model(project_dir)
    if not ok:
        return False, message
    ok, message = verify_imports()
    if not ok:
        return False, message
    return True, "Environment ready."


def _dry_run_discovery(project_dir: Path) -> None:
    data_root = project_dir / "data" / "aSMA_DAPI_plates"
    if not data_root.is_dir():
        return
    from dapi_norm.user_cellpose_batch import discover_acquisitions  # noqa: PLC0415

    acquisitions = discover_acquisitions(data_root)
    print(f"Discovered {len(acquisitions)} acquisition folder(s):")
    for acquisition in acquisitions:
        print(
            f"  - {acquisition.plate_name} / {acquisition.display_name}: "
            f"{acquisition.image_count} image pairs"
        )


def finish_setup(project_dir: Path) -> int:
    """Post-venv setup steps. Identical on macOS and Windows."""
    print("==> Ensuring Cellpose cpsam_v2 model is cached and verified")
    download_model(project_dir)
    ok, message = verify_model(project_dir)
    print(f"    {message}")
    if not ok:
        print(f"ERROR: {message}", file=sys.stderr)
        return 1

    print("==> Verifying imports")
    ok, message = verify_imports()
    print(f"    {message}")
    if not ok:
        print(f"ERROR: {message}", file=sys.stderr)
        return 1

    print("==> Discovery dry run on bundled plate data (if present)")
    _dry_run_discovery(project_dir)

    print("==> Setup verification complete")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["check", "finish-setup", "verify-model"])
    parser.add_argument("--project-dir", required=True)
    args = parser.parse_args(argv)
    project_dir = Path(args.project_dir).resolve()

    if args.command == "check":
        ok, message = env_ready(project_dir)
        print(message)
        return 0 if ok else 1
    if args.command == "verify-model":
        ok, message = verify_model(project_dir)
        print(message)
        return 0 if ok else 1
    return finish_setup(project_dir)


if __name__ == "__main__":
    raise SystemExit(main())
