#!/bin/zsh
set -u

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
SETUP_SCRIPT="$PROJECT_DIR/scripts/setup_cellpose_pipeline_macos.sh"

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

if [[ ! -f "$SETUP_SCRIPT" ]]; then
  fail_with_dialog "The Mac setup helper was not found under scripts/setup_cellpose_pipeline_macos.sh."
fi

echo "Project: $PROJECT_DIR"
echo ""
echo "Setting up the project-local Cellpose environment. This can take a while."
echo ""

/bin/bash "$SETUP_SCRIPT" "$PROJECT_DIR"
STATUS=$?

if [[ "$STATUS" -ne 0 ]]; then
  fail_with_dialog "Setup failed. The Terminal window contains the error details."
fi

osascript -e "display dialog \"Setup finished. You can now double-click Run Cellpose DAPI aSMA Pipeline.command.\" buttons {\"OK\"} default button \"OK\" with icon note" >/dev/null 2>&1 || true
echo ""
echo "Setup finished. You can now double-click:"
echo "Run Cellpose DAPI aSMA Pipeline.command"
echo ""
echo "Press any key to close this window."
read -k 1
