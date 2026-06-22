# Installation Plan

This document describes how I would set up the project. No heavy install has been run yet.

## Common Policy

- Use a project-local environment.
- Do not install into system Python.
- Use Python 3.12 for the Cellpose-first path.
- Keep Cellpose model weights outside git.
- Verify CPU/GPU availability before running real data.
- Install StarDist, CellProfiler, Docker, or QuPath only after Cellpose v1 works.

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

Optional: keep Cellpose model downloads inside the project folder rather than the home directory.

```bash
mkdir -p .models/cellpose
export CELLPOSE_LOCAL_MODELS_PATH="$PWD/.models/cellpose"
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

The Cellpose docs give this as the shape of the command for a CUDA-enabled PyTorch install:

```powershell
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu126
```

Use the current command from the PyTorch selector for the actual PC. Newer GPUs or drivers may use a newer CUDA wheel.

Then install the project:

```powershell
python -m pip install -e ".[dev,cellpose]"
```

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

## Cellpose Model Download Behavior

Cellpose downloads pretrained model weights on first use. For reproducibility:

- Record the Cellpose package version in `output/logs/config_resolved.yaml`.
- Record the model name and resolved model path.
- Keep model weights out of git.
- Prefer one fixed model name for v1, initially `nuclei`.

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

