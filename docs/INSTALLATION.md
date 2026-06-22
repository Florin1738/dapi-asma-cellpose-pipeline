# Installation and Runtime Setup

This document describes the current Cellpose-first setup. The Mac setup has been run successfully for the local sample data; the Windows NVIDIA path is the matching setup plan for larger GPU-backed batches.

## Common Policy

- Use a project-local environment.
- Do not install into system Python.
- Use Python 3.12 for the Cellpose-first path.
- Keep Cellpose model weights outside git.
- Verify CPU/GPU availability before running real data.
- Install StarDist, CellProfiler, Docker, or QuPath only after Cellpose v1 works.
- Use Cellpose 4 with `cpsam_v2` for the first runnable workflow.

## Mac Apple Silicon Setup

This current Mac has Homebrew and `uv`, no conda/mamba, no Docker, and global Python 3.14.5. Python 3.14 is too new for a conservative bioimage stack, so use `uv` to provision Python 3.12.

```bash
cd /Users/florinselaru/Documents/DAPI_intensity_quantification_tata
uv python install 3.12
uv venv --python 3.12 .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev,cellpose]'
```

For the closest reproduction of the successful 2026-06-22 Mac pilot run, use the dated constraints file:

```bash
python -m pip install -c constraints/cellpose-mac-2026-06-22.txt -e '.[dev,cellpose]'
```

Optional: keep Cellpose model downloads inside the project folder rather than the home directory.

```bash
mkdir -p .models/cellpose
export CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose"
```

The current Mac model cache is:

```text
.models/cellpose/cpsam_v2
size: 1.1G
sha256: 0f1cc3f7ecdd8a037a57c6c48d9d8921391be4cbce3fa9f13c3e3a2e1253c667
```

If the automatic Cellpose model download is unreliable, this command downloads the model directly into the project cache:

```bash
mkdir -p .models/cellpose
curl -L --fail --retry 8 --retry-all-errors --retry-delay 5 \
  --connect-timeout 30 \
  -o .models/cellpose/cpsam_v2.download \
  https://huggingface.co/mouseland/cellpose-sam/resolve/main/cpsam_v2
mv .models/cellpose/cpsam_v2.download .models/cellpose/cpsam_v2
```

Verify Python and PyTorch:

```bash
python - <<'PY'
import platform
import torch

print("python ok")
print("machine:", platform.machine())
print("torch:", torch.__version__)
print("mps_available:", torch.backends.mps.is_available())
PY
```

Verify Cellpose imports:

```bash
python - <<'PY'
import cellpose
print("cellpose:", getattr(cellpose, "__version__", "unknown"))
PY
```

If MPS is unavailable, the first local version can still run on CPU for small synthetic and pilot images. Large batch runs should move to the Windows NVIDIA machine.

Run the current count-only pipeline:

```bash
CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose" \
  .venv/bin/python scripts/run_cellpose_counts.py \
  --input ApYYM20AGGSMA_02 \
  --output output/cellpose_counts \
  --channel CH4 \
  --model cpsam_v2 \
  --gpu
```

## Windows NVIDIA Setup

Use this path if you want to run larger batches on a Windows PC with an NVIDIA GPU.

1. Install or update the NVIDIA driver from NVIDIA.
2. Install Python 3.12 from python.org or another trusted Python manager.
3. Create a project virtual environment.

PowerShell:

```powershell
cd C:\path\to\DAPI_intensity_quantification_tata
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
```

Install PyTorch using the official selector at https://pytorch.org/get-started/locally/. Select:

```text
OS: Windows
Package: Pip
Language: Python
Compute Platform: CUDA version supported by the installed NVIDIA driver
```

The current PyTorch selector should generate the exact command for the PC. The command will look like this, but the CUDA wheel tag may differ:

```powershell
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128
```

Use the live selector output for the actual PC. Newer GPUs or drivers may use a newer CUDA wheel, and older drivers may require a different supported CUDA wheel.

Then install the project:

```powershell
python -m pip install -e ".[dev,cellpose]"
```

Do not blindly apply `constraints/cellpose-mac-2026-06-22.txt` on Windows before installing the CUDA-enabled PyTorch wheel chosen by the official selector. That file records the successful Mac pilot stack and is mainly useful for comparing package versions.

Verify CUDA:

```powershell
python - <<'PY'
import torch
print("torch:", torch.__version__)
print("cuda_available:", torch.cuda.is_available())
print("cuda_device_count:", torch.cuda.device_count())
if torch.cuda.is_available():
    print("cuda_device:", torch.cuda.get_device_name(0))
PY
```

Verify Cellpose:

```powershell
python - <<'PY'
import cellpose
print("cellpose:", getattr(cellpose, "__version__", "unknown"))
PY
```

Set a project-local model cache and run the same count command from PowerShell:

```powershell
New-Item -ItemType Directory -Force .models\cellpose | Out-Null
$env:CELLPOSE_LOCAL_MODELS_PATH = "$PWD\.models\cellpose"
python scripts\run_cellpose_counts.py `
  --input ApYYM20AGGSMA_02 `
  --output output\cellpose_counts `
  --channel CH4 `
  --model cpsam_v2 `
  --gpu
```

## Cellpose Model Download Behavior

Cellpose downloads pretrained model weights on first use. For reproducibility:

- Record the Cellpose package version in the run config, currently `output/cellpose_counts/logs/config_resolved.yaml`.
- Record the model name and resolved model path.
- Keep model weights out of git.
- Prefer one fixed model name for v1, currently `cpsam_v2`.

Cellpose 4.2.1.1 exposes the current built-in model names `cpsam_v2`, `cpdino`, `cpdino-vitb`, and `cpsam` in this environment. The old `nuclei` model name from earlier Cellpose workflows is not the active v4 path here.

## Source Links Checked

- Cellpose installation docs: https://cellpose.readthedocs.io/en/latest/installation.html
- Cellpose model docs: https://cellpose.readthedocs.io/en/latest/models.html
- PyTorch local install selector: https://pytorch.org/get-started/locally/

## Optional GUI Tools

### Cellpose GUI

Install only if interactive review/training is needed:

```bash
python -m pip install 'cellpose[gui]'
python -m cellpose
```

### QuPath

Install separately from the Python environment if whole-slide images, tissue annotations, or interactive ROI work become central. On Apple Silicon, use the arm64 package for speed unless file-format compatibility requires the Intel build.

### CellProfiler

Keep CellProfiler separate from this CLI project. It can be used later as a reference GUI pipeline or plugin-based comparison workflow.
