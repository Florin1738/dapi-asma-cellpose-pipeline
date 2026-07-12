#!/bin/zsh
set -u

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON_BIN="$PROJECT_DIR/.venv/bin/python"
SETUP_SCRIPT="$PROJECT_DIR/scripts/setup_cellpose_pipeline_macos.sh"
GUI_SCRIPT="$PROJECT_DIR/scripts/cellpose_gui.py"

fail_with_dialog() {
  local message="$1"
  echo ""
  echo "$message"
  osascript -e "display dialog \"$message\" buttons {\"OK\"} default button \"OK\" with icon stop" >/dev/null 2>&1 || true
  echo ""
  echo "Press any key to close this window."
  read -k 1
  exit 1
}

env_ready() {
  [[ -x "$PYTHON_BIN" ]] || return 1
  # Shared cross-platform readiness check (venv + model checksum + imports).
  CELLPOSE_LOCAL_MODELS_PATH="$PROJECT_DIR/.models/cellpose" \
    "$PYTHON_BIN" "$PROJECT_DIR/scripts/pipeline_env.py" check --project-dir "$PROJECT_DIR" >/dev/null 2>&1
}

# Self-heal: install the environment automatically if it is missing or broken.
if ! env_ready; then
  if [[ -d "$PROJECT_DIR/.venv" ]]; then
    HEAL_MSG="The analysis environment needs to be repaired or reinstalled. This can happen after an interrupted setup or an update. Reinstalling now (a few minutes, needs internet). The app opens automatically when it finishes."
  else
    HEAL_MSG="First-time setup will now install the analysis environment (a few minutes, needs internet). The app opens automatically when it finishes."
  fi
  echo "$HEAL_MSG"
  echo ""
  osascript -e "display dialog \"$HEAL_MSG\" buttons {\"OK\"} default button \"OK\" with icon note" >/dev/null 2>&1 || true

  [[ -f "$SETUP_SCRIPT" ]] || fail_with_dialog "The setup helper was not found under scripts/setup_cellpose_pipeline_macos.sh."
  /bin/bash "$SETUP_SCRIPT" "$PROJECT_DIR"
  if [[ "$?" -ne 0 ]]; then
    fail_with_dialog "Setup failed. The Terminal window contains the error details. If it mentions the internet, check your connection and try again."
  fi
fi

if ! env_ready; then
  fail_with_dialog "The analysis environment is still not ready after setup. The Terminal window contains details."
fi

[[ -f "$GUI_SCRIPT" ]] || fail_with_dialog "The application file scripts/cellpose_gui.py was not found."

export CELLPOSE_LOCAL_MODELS_PATH="$PROJECT_DIR/.models/cellpose"

echo "Opening the Cellpose DAPI / aSMA app…"
"$PYTHON_BIN" "$GUI_SCRIPT" --project-dir "$PROJECT_DIR"
STATUS=$?

if [[ "$STATUS" -ne 0 ]]; then
  fail_with_dialog "The app could not start. The Terminal window contains the error details."
fi
