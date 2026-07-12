param(
    [Parameter(Mandatory = $true)]
    [string]$ProjectDir
)

$ErrorActionPreference = "Stop"
$PythonVersion = "3.12"
# Model URL and checksum live in scripts\pipeline_env.py (single source of truth).

function Write-Step {
    param([string]$Message)
    Write-Host ""
    Write-Host "==> $Message" -ForegroundColor Cyan
}

function Fail {
    param([string]$Message)
    Write-Host ""
    Write-Host "ERROR: $Message" -ForegroundColor Red
    Read-Host "Press Enter to close this window"
    exit 1
}

function Find-Uv {
    $candidates = @(
        (Join-Path $ProjectDir ".tools\uv\uv.exe"),
        (Join-Path $env:USERPROFILE ".local\bin\uv.exe"),
        (Join-Path $env:USERPROFILE ".cargo\bin\uv.exe")
    )
    foreach ($candidate in $candidates) {
        if (Test-Path $candidate -PathType Leaf) {
            return $candidate
        }
    }
    $command = Get-Command uv -ErrorAction SilentlyContinue
    if ($command) {
        return $command.Source
    }
    return $null
}

function Install-Uv {
    Write-Step "Installing uv for this user account"
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
    Invoke-RestMethod https://astral.sh/uv/install.ps1 | Invoke-Expression
}

$ProjectDir = (Resolve-Path $ProjectDir).Path
Set-Location $ProjectDir

$Uv = Find-Uv
if (-not $Uv) {
    Install-Uv
    $Uv = Find-Uv
}
if (-not $Uv) {
    Fail "uv was not found after installation."
}

Write-Step "Using uv: $Uv"
& $Uv --version

Write-Step "Creating project-local Python $PythonVersion environment"
& $Uv python install $PythonVersion
& $Uv venv --python $PythonVersion .venv
$PythonBin = Join-Path $ProjectDir ".venv\Scripts\python.exe"
if (-not (Test-Path $PythonBin -PathType Leaf)) {
    Fail "Expected Python at .venv\Scripts\python.exe after setup."
}

Write-Step "Installing project and Cellpose into .venv"
$Wheelhouse = Join-Path $ProjectDir "wheelhouse"
if ((Test-Path $Wheelhouse -PathType Container) -and (Get-ChildItem $Wheelhouse -Filter *.whl -ErrorAction SilentlyContinue | Select-Object -First 1)) {
    & $Uv pip install --python $PythonBin --no-index --find-links $Wheelhouse -e ".[cellpose]"
} else {
    & $Uv pip install --python $PythonBin -e ".[cellpose]"
}

Write-Step "Finalizing setup: model download, checksum, import and discovery checks"
# Shared cross-platform logic in scripts\pipeline_env.py. macOS and Windows run
# this exact same code path, so testing setup on one OS exercises the
# substantive steps for the other.
$env:CELLPOSE_LOCAL_MODELS_PATH = Join-Path $ProjectDir ".models\cellpose"
& $PythonBin scripts\pipeline_env.py finish-setup --project-dir $ProjectDir
if ($LASTEXITCODE -ne 0) {
    Fail "Setup finalization failed. See the messages above. If it mentions the internet, check your connection and run setup again."
}

Write-Step "Setup complete"
Write-Host "Double-click Run Cellpose DAPI aSMA Pipeline Windows.cmd to process data."
Read-Host "Press Enter to close this window"
