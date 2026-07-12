#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
PREFIX="${CELLPROFILER_ARM64_PREFIX:-$REPO_ROOT/.cellprofiler-arm64}"
MICROMAMBA="${MICROMAMBA_EXE:-$HOME/.local/bin/micromamba}"

export JAVA_HOME="$PREFIX/lib/jvm"
export PATH="$JAVA_HOME/bin:$PATH"
export CELLPROFILER_ARM64_PREFIX_RESOLVED="$PREFIX"

if [ "${1:-}" = "--runtime-info" ]; then
  exec "$MICROMAMBA" run -p "$PREFIX" python - <<'PY'
import os
import platform
import sys

import cellprofiler
import cellprofiler_core

print(f"cellprofiler_arm64_prefix={os.environ.get('CELLPROFILER_ARM64_PREFIX_RESOLVED', '')}")
print(f"python_executable={sys.executable}")
print(f"platform_machine={platform.machine()}")
print(f"python_version={sys.version.split()[0]}")
print(f"cellprofiler_version={cellprofiler.__version__}")
print(f"cellprofiler_core_version={cellprofiler_core.__version__}")
PY
fi

exec "$MICROMAMBA" run -p "$PREFIX" python -m cellprofiler "$@"
