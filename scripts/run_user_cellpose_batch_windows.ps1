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
            $Message,
            "Cellpose batch pipeline",
            [System.Windows.Forms.MessageBoxButtons]::OK,
            [System.Windows.Forms.MessageBoxIcon]::Error
        ) | Out-Null
    } catch {
    }
    Write-Host ""
    Read-Host "Press Enter to close this window"
    exit 1
}

function Select-Folder {
    param([string]$Description)
    Add-Type -AssemblyName System.Windows.Forms
    $dialog = New-Object System.Windows.Forms.FolderBrowserDialog
    $dialog.Description = $Description
    $dialog.ShowNewFolderButton = $true
    $result = $dialog.ShowDialog()
    if ($result -ne [System.Windows.Forms.DialogResult]::OK) {
        exit 0
    }
    return $dialog.SelectedPath
}

$ProjectDir = (Resolve-Path $ProjectDir).Path
$PythonBin = Join-Path $ProjectDir ".venv\Scripts\python.exe"
$ModelPath = Join-Path $ProjectDir ".models\cellpose\cpsam_v2"
$Runner = Join-Path $ProjectDir "scripts\run_user_cellpose_batch.py"

if (-not (Test-Path $PythonBin -PathType Leaf)) {
    Show-ErrorAndExit "The Windows project Python environment was not found at .venv\Scripts\python.exe. Use the prepared Windows project folder or run Windows setup first."
}

if (-not (Test-Path $ModelPath -PathType Leaf)) {
    Show-ErrorAndExit "The Cellpose cpsam_v2 model was not found in .models\cellpose. Use the prepared Windows project folder or copy the prepared model cache."
}

if (-not (Test-Path $Runner -PathType Leaf)) {
    Show-ErrorAndExit "The batch runner script was not found under scripts\run_user_cellpose_batch.py."
}

& $PythonBin -c "import cellpose; import dapi_norm.user_cellpose_batch" *> $null
if ($LASTEXITCODE -ne 0) {
    Show-ErrorAndExit "The Windows Python environment exists, but Cellpose or this project package could not be imported. Use the prepared Windows project folder or run Windows setup first."
}

$InputFolder = Select-Folder "Select the folder that contains the plate or acquisition image data"
$OutputParent = Select-Folder "Select where the results folder should be created"
$RunStamp = Get-Date -Format "yyyyMMdd_HHmmss"
$RunOutput = Join-Path $OutputParent "Cellpose_DAPI_aSMA_$RunStamp"

$env:CELLPOSE_LOCAL_MODELS_PATH = Join-Path $ProjectDir ".models\cellpose"

Write-Host "Project: $ProjectDir"
Write-Host "Input:   $InputFolder"
Write-Host "Output:  $RunOutput"
Write-Host ""
Write-Host "Running Cellpose batch pipeline. This can take a while for full plates."
Write-Host ""

& $PythonBin $Runner --input $InputFolder --output $RunOutput
$Status = $LASTEXITCODE
if ($Status -ne 0) {
    Show-ErrorAndExit "The Cellpose batch pipeline failed. This window contains the error details."
}

$FinalDir = Join-Path $RunOutput "final"
$SummaryHtml = Join-Path $FinalDir "START_HERE_RUN_SUMMARY.html"
if (Test-Path $FinalDir) {
    Start-Process -FilePath explorer.exe -ArgumentList @("`"$FinalDir`"")
}
if (Test-Path $SummaryHtml) {
    Start-Process $SummaryHtml
}

try {
    Add-Type -AssemblyName System.Windows.Forms
    [System.Windows.Forms.MessageBox]::Show(
        "Cellpose batch pipeline finished. Results were written to $FinalDir",
        "Cellpose batch pipeline",
        [System.Windows.Forms.MessageBoxButtons]::OK,
        [System.Windows.Forms.MessageBoxIcon]::Information
    ) | Out-Null
} catch {
}

Write-Host ""
Write-Host "Done. Results: $FinalDir"
Read-Host "Press Enter to close this window"
