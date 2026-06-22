# Software Research

Research snapshot date: 2026-06-22.

## Recommendation

Use Cellpose as the v1 segmentation engine. It is the best fit for a folder-based Python CLI that can run on this Apple Silicon Mac now and later on a Windows workstation with an NVIDIA GPU.

StarDist is worth testing later as a comparison backend, especially for compact nuclei, but it depends on TensorFlow and has more fragile setup paths across platforms. CellProfiler and QuPath remain useful GUI/reference tools but should not be required for v1.

## Cellpose

Cellpose is the first implementation target.

Relevant facts from official sources:

- The Cellpose GitHub install guide recommends Python 3.12, with Python 3.9 and 3.11 also expected to work.
- Cellpose can be installed into a conda environment or a normal Python virtual environment with `python -m pip install cellpose`.
- The optional GUI install is `python -m pip install 'cellpose[gui]'`.
- On Apple Silicon, Cellpose can use PyTorch MPS acceleration when run with `--use_gpu` and a compatible PyTorch install.
- Pretrained Cellpose model weights download automatically on first use. The default model location is under `$HOME/.cellpose/models`, unless `CELLPOSE_LOCAL_MODELS_PATH` is set.
- For Windows/Linux NVIDIA GPU use, Cellpose points users to the PyTorch CUDA install path and gives an example of installing CUDA-enabled `torch` and `torchvision`.

Sources:

- [Cellpose GitHub install section](https://github.com/MouseLand/cellpose)
- [Cellpose installation docs](https://cellpose.readthedocs.io/en/latest/installation.html)
- [Cellpose model docs](https://cellpose.readthedocs.io/en/latest/models.html)

## PyTorch Platform Layer

Cellpose depends on PyTorch for acceleration.

Apple Silicon:

- Apple documents PyTorch MPS as the Metal backend for GPU acceleration on Macs.
- Current Apple guidance lists Apple Silicon, macOS 14+, Python 3.10+, and Xcode command-line tools as requirements for the latest stable PyTorch at the time of this research.
- Verification should check `torch.backends.mps.is_available()`.

Windows NVIDIA:

- PyTorch officially recommends Windows systems with NVIDIA GPUs for CUDA acceleration.
- The official install selector should be used on the target Windows PC to choose the correct CUDA wheel for the installed driver/GPU.
- Verification should check `torch.cuda.is_available()`.

Sources:

- [PyTorch local install guide](https://pytorch.org/get-started/locally/)
- [Apple PyTorch on Metal guide](https://developer.apple.com/metal/pytorch/)

## StarDist

StarDist remains a good scientific comparison backend but not the v1 blocker.

Relevant facts:

- StarDist provides pretrained-style instance segmentation workflows for star-convex objects such as nuclei.
- StarDist includes matching utilities for segmentation metrics such as true positives, false positives, false negatives, precision, recall, and F1.
- StarDist has Apple Silicon notes, but the recommended path uses an arm64 conda environment plus TensorFlow-specific setup.
- On Windows, StarDist may need Microsoft C++ build tools if wheels are unavailable or compilation is required.

Source:

- [StarDist GitHub README](https://github.com/stardist/stardist)

## CellProfiler

CellProfiler is useful as a GUI/reference workflow but not the default for this project.

Relevant facts:

- CellProfiler plugins can expose RunCellpose and RunStarDist modules.
- The plugin docs say CellProfiler detects plugins from a configured plugin folder when dependencies are present.
- RunCellpose can use a Docker image so the plugin can run without installing Cellpose directly into the CellProfiler Python environment.
- The plugin docs caution that having Cellpose and StarDist in the same Python environment can be troublesome and recommend separate environments when needed.

Local note: Homebrew has a CellProfiler cask on this Mac, but it is Intel/Rosetta and marked deprecated by Homebrew because of Gatekeeper issues. I would not use that as the primary path.

Sources:

- [CellProfiler plugin usage docs](https://plugins.cellprofiler.org/using_plugins.html)
- [RunCellpose plugin docs](https://plugins.cellprofiler.org/runcellpose.html)
- [Supported CellProfiler plugins](https://plugins.cellprofiler.org/supported_plugins.html)

## QuPath

QuPath is useful for whole-slide/tissue workflows, ROI review, and visual inspection, especially if large tissue sections become central.

Relevant facts:

- QuPath provides both Intel and Apple Silicon macOS builds.
- The Apple Silicon build is faster on recent Macs but has Bio-Formats caveats for some `.czi` and `.ndpi/.ndpis` files.
- QuPath docs recommend OME-TIFF as a well-supported open format when possible.

Sources:

- [QuPath installation docs](https://qupath.readthedocs.io/en/stable/docs/intro/installation.html)
- [QuPath image format docs](https://qupath.readthedocs.io/en/latest/docs/intro/formats.html)

## Decision

The v1 project should install and validate this stack first:

```text
Python 3.12
project-local virtual environment
scientific base packages
Cellpose
PyTorch acceleration appropriate to the machine
```

The first runnable version should not require StarDist, CellProfiler, Docker, QuPath, or proprietary image readers.

