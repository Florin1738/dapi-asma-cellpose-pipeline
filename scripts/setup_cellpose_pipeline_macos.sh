#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PYTHON_VERSION="3.12"
# Model URL and checksum live in scripts/pipeline_env.py (single source of truth).

cd "$PROJECT_DIR"

log() {
  printf '\n==> %s\n' "$1"
}

fail() {
  printf '\nERROR: %s\n' "$1" >&2
  exit 1
}

find_uv() {
  if [[ -x "$PROJECT_DIR/.tools/uv/uv" ]]; then
    printf '%s\n' "$PROJECT_DIR/.tools/uv/uv"
    return 0
  fi
  if command -v uv >/dev/null 2>&1; then
    command -v uv
    return 0
  fi
  if [[ -x "$HOME/.local/bin/uv" ]]; then
    printf '%s\n' "$HOME/.local/bin/uv"
    return 0
  fi
  if [[ -x "$HOME/.cargo/bin/uv" ]]; then
    printf '%s\n' "$HOME/.cargo/bin/uv"
    return 0
  fi
  return 1
}

install_uv() {
  command -v curl >/dev/null 2>&1 || fail "curl is required to install uv and download the model."
  log "Installing uv for this user account"
  curl -LsSf https://astral.sh/uv/install.sh | sh
}

UV_BIN="$(find_uv || true)"
if [[ -z "$UV_BIN" ]]; then
  install_uv
  UV_BIN="$(find_uv || true)"
fi
[[ -n "$UV_BIN" ]] || fail "uv was not found after installation."

log "Using uv: $UV_BIN"
"$UV_BIN" --version

log "Creating project-local Python $PYTHON_VERSION environment"
"$UV_BIN" python install "$PYTHON_VERSION"
"$UV_BIN" venv --python "$PYTHON_VERSION" .venv
PYTHON_BIN="$PROJECT_DIR/.venv/bin/python"
[[ -x "$PYTHON_BIN" ]] || fail "Expected Python at .venv/bin/python after setup."

log "Installing project and Cellpose into .venv"
WHEELHOUSE="$PROJECT_DIR/wheelhouse"
if find "$WHEELHOUSE" -maxdepth 1 -name '*.whl' -print -quit >/dev/null 2>&1; then
  "$UV_BIN" pip install --python "$PYTHON_BIN" --no-index --find-links "$WHEELHOUSE" -e '.[cellpose]'
elif [[ -f "$PROJECT_DIR/constraints/cellpose-mac-2026-06-22.txt" ]]; then
  "$UV_BIN" pip install --python "$PYTHON_BIN" -c "$PROJECT_DIR/constraints/cellpose-mac-2026-06-22.txt" -e '.[cellpose]'
else
  "$UV_BIN" pip install --python "$PYTHON_BIN" -e '.[cellpose]'
fi

log "Finalizing setup: model download, checksum, import and discovery checks"
# Shared cross-platform logic in scripts/pipeline_env.py. macOS and Windows run
# this exact same code path, so testing setup on one OS exercises the
# substantive steps for the other.
CELLPOSE_LOCAL_MODELS_PATH="$PROJECT_DIR/.models/cellpose" "$PYTHON_BIN" \
  "$PROJECT_DIR/scripts/pipeline_env.py" finish-setup --project-dir "$PROJECT_DIR" \
  || fail "Setup finalization failed. See the messages above. If it mentions the internet, check your connection and run setup again."

log "Setup complete"
printf 'Double-click "Run Cellpose DAPI aSMA Pipeline.command" to process data.\n'
