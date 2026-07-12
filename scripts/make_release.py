#!/usr/bin/env python
"""Build the Mac and Windows release zips a nontechnical user downloads.

Each zip unzips to a single folder containing exactly what that operating
system needs: the shared source/scripts/docs, the OS-specific double-click
launchers, and READ ME FIRST.txt. Private microscopy data, generated outputs,
the virtual environment, and the model cache are always excluded.

Usage:

    python scripts/make_release.py                 # online installer zips (small)
    python scripts/make_release.py --with-model    # also bundle the Cellpose model
    python scripts/make_release.py --with-wheelhouse  # also bundle wheelhouse/ if present
    python scripts/make_release.py --version v1.0   # tag the zip names

Offline zips (--with-model and/or --with-wheelhouse) let the first launch run
without internet, at the cost of a much larger download.
"""

from __future__ import annotations

import argparse
import fnmatch
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

# Shared content included in both OS zips.
SHARED_TOP = [
    "src",
    "scripts",
    "docs",
    "constraints",
    "configs",
    "examples",
    "tests",
    "pyproject.toml",
    "README.md",
    "AGENTS.md",
    "PROJECT_CONTENTS.md",
    "READ ME FIRST.txt",
]

MAC_LAUNCHERS = [
    "Run Cellpose DAPI aSMA Pipeline.command",
    "Setup Cellpose DAPI aSMA Pipeline.command",
]
WINDOWS_LAUNCHERS = [
    "Run Cellpose DAPI aSMA Pipeline Windows.cmd",
    "Setup Cellpose DAPI aSMA Pipeline Windows.cmd",
]

# Never include these anywhere (heavy, machine-specific, or lab-private).
EXCLUDE_DIR_NAMES = {
    "__pycache__",
    ".git",
    ".venv",
    ".models",
    ".pytest_cache",
    ".ruff_cache",
    "node_modules",
    "data",
    "output",
    "reports",
    "manual_validation",
}
# Never include files matching these patterns — a hard guard against leaking
# microscopy data or generated artifacts even if they sit inside an included dir.
EXCLUDE_FILE_GLOBS = [
    "*.tif",
    "*.tiff",
    "*.npy",
    "*.png",
    "*.jpg",
    "*.jpeg",
    ".DS_Store",
    "*.pyc",
]


def _is_excluded_file(path: Path) -> bool:
    return any(fnmatch.fnmatch(path.name, pattern) for pattern in EXCLUDE_FILE_GLOBS)


def _iter_files(top: Path):
    if top.is_file():
        if not _is_excluded_file(top):
            yield top
        return
    for child in sorted(top.rglob("*")):
        if child.is_dir():
            if child.name in EXCLUDE_DIR_NAMES:
                # Skip the whole subtree by not descending; rglob still yields
                # descendants, so guard each file below too.
                continue
            continue
        if any(part in EXCLUDE_DIR_NAMES for part in child.relative_to(REPO_ROOT).parts):
            continue
        if _is_excluded_file(child):
            continue
        yield child


def build_zip(
    out_path: Path,
    release_name: str,
    launchers: list[str],
    with_model: bool,
    with_wheelhouse: bool,
) -> None:
    added: set[str] = set()
    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as zf:

        def add(source: Path) -> None:
            rel = source.relative_to(REPO_ROOT)
            arc = f"{release_name}/{rel.as_posix()}"
            if arc in added:
                return
            zf.write(source, arc)
            added.add(arc)

        for name in SHARED_TOP + launchers:
            top = REPO_ROOT / name
            if not top.exists():
                print(f"  ! missing, skipped: {name}")
                continue
            for f in _iter_files(top):
                add(f)

        if with_model:
            model = REPO_ROOT / ".models" / "cellpose" / "cpsam_v2"
            if model.is_file():
                add(model)
                print("  + bundled model cache (.models/cellpose/cpsam_v2)")
            else:
                print("  ! --with-model requested but .models/cellpose/cpsam_v2 not found")
        if with_wheelhouse:
            wh = REPO_ROOT / "wheelhouse"
            if wh.is_dir():
                for f in _iter_files(wh):
                    add(f)
                print("  + bundled wheelhouse/")
            else:
                print("  ! --with-wheelhouse requested but wheelhouse/ not found")

    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"  -> {out_path.name}  ({len(added)} files, {size_mb:.1f} MB)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default="", help="Optional version tag in the zip name.")
    parser.add_argument("--with-model", action="store_true", help="Bundle the Cellpose model.")
    parser.add_argument("--with-wheelhouse", action="store_true", help="Bundle wheelhouse/.")
    parser.add_argument("--out", default="dist", help="Output directory (default: dist).")
    args = parser.parse_args()

    suffix = f"-{args.version}" if args.version else ""
    offline = args.with_model or args.with_wheelhouse
    tier = "-offline" if offline else ""

    dist = REPO_ROOT / args.out
    dist.mkdir(parents=True, exist_ok=True)

    targets = [
        (f"cellpose-dapi-asma-pipeline-mac{tier}{suffix}", MAC_LAUNCHERS),
        (f"cellpose-dapi-asma-pipeline-windows{tier}{suffix}", WINDOWS_LAUNCHERS),
    ]
    for release_name, launchers in targets:
        print(f"Building {release_name} ...")
        build_zip(
            dist / f"{release_name}.zip",
            release_name,
            launchers,
            args.with_model,
            args.with_wheelhouse,
        )

    print(f"\nDone. Zips are in {dist}")
    if not offline:
        print("These are ONLINE installers: the first launch downloads the model + deps.")
        print("For an offline zip, add --with-model (and --with-wheelhouse if you have one).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
