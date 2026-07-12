param(
    [Parameter(Mandatory = $true)]
    [string]$ProjectDir
)

$ErrorActionPreference = "Stop"

function Show-ErrorAndExit {
    param([string]$Message)
    Write-Host ""
    Write-Host $Message -ForegroundColor Red
    try {
        Add-Type -AssemblyName System.Windows.Forms
        [System.Windows.Forms.MessageBox]::Show(
            $Message, "Cellpose DAPI / aSMA Pipeline",
            [System.Windows.Forms.MessageBoxButtons]::OK,
            [System.Windows.Forms.MessageBoxIcon]::Error
        ) | Out-Null
    } catch {
    }
    Write-Host ""
    Read-Host "Press Enter to close this window"
    exit 1
}

$ProjectDir = (Resolve-Path $ProjectDir).Path
$PythonBin  = Join-Path $ProjectDir ".venv\Scripts\python.exe"
$PythonwBin = Join-Path $ProjectDir ".venv\Scripts\pythonw.exe"
$SetupScript = Join-Path $ProjectDir "scripts\setup_cellpose_pipeline_windows.ps1"
$GuiScript  = Join-Path $ProjectDir "scripts\cellpose_gui.py"

function Test-EnvReady {
    if (-not (Test-Path $PythonBin -PathType Leaf)) { return $false }
    # Shared cross-platform readiness check (venv + model checksum + imports).
    $env:CELLPOSE_LOCAL_MODELS_PATH = Join-Path $ProjectDir ".models\cellpose"
    & $PythonBin (Join-Path $ProjectDir "scripts\pipeline_env.py") check --project-dir $ProjectDir *> $null
    return ($LASTEXITCODE -eq 0)
}

# Self-heal: install the environment automatically if it is missing or broken.
if (-not (Test-EnvReady)) {
    if (Test-Path (Join-Path $ProjectDir ".venv") -PathType Container) {
        $HealMsg = "The analysis environment needs to be repaired or reinstalled. This can happen after an interrupted setup or an update. Reinstalling now (a few minutes, needs internet). The app opens automatically when it finishes."
    } else {
        $HealMsg = "First-time setup will now install the analysis environment (a few minutes, needs internet). The app opens automatically when it finishes."
    }
    Write-Host $HealMsg
    Write-Host ""
    try {
        Add-Type -AssemblyName System.Windows.Forms
        [System.Windows.Forms.MessageBox]::Show(
            $HealMsg,
            "Cellpose DAPI / aSMA Pipeline",
            [System.Windows.Forms.MessageBoxButtons]::OK,
            [System.Windows.Forms.MessageBoxIcon]::Information
        ) | Out-Null
    } catch {
    }

    if (-not (Test-Path $SetupScript -PathType Leaf)) {
        Show-ErrorAndExit "The setup helper was not found under scripts\setup_cellpose_pipeline_windows.ps1."
    }
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File $SetupScript -ProjectDir $ProjectDir
    if ($LASTEXITCODE -ne 0) {
        Show-ErrorAndExit "Setup failed. This window contains the error details. If it mentions the internet, check your connection and try again."
    }
}

if (-not (Test-EnvReady)) {
    Show-ErrorAndExit "The analysis environment is still not ready after setup. This window contains details."
}

if (-not (Test-Path $GuiScript -PathType Leaf)) {
    Show-ErrorAndExit "The application file scripts\cellpose_gui.py was not found."
}

$env:CELLPOSE_LOCAL_MODELS_PATH = Join-Path $ProjectDir ".models\cellpose"

# Launch the GUI with pythonw so no extra console window lingers.
$Launcher = $PythonwBin
if (-not (Test-Path $Launcher -PathType Leaf)) { $Launcher = $PythonBin }

Write-Host "Opening the Cellpose DAPI / aSMA app…"
& $Launcher $GuiScript --project-dir $ProjectDir
$Status = $LASTEXITCODE
if ($Status -ne 0) {
    Show-ErrorAndExit "The app could not start. This window contains the error details."
}
