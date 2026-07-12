#!/usr/bin/env python
"""Single-window desktop GUI for the Cellpose DAPI/aSMA pipeline.

Designed for non-technical users. Pick an input image folder and an output
folder, press Run, watch progress. Scientific parameters are hidden behind a
collapsible "Advanced options" panel and default to the validated values.

Launched by the double-click run launchers with:

    python scripts/cellpose_gui.py --project-dir /path/to/project

It runs the existing ``scripts/run_user_cellpose_batch.py`` as a subprocess so
the analysis contract and outputs are identical to the command-line path.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import math
import os
import queue
import re
import subprocess
import sys
import threading
from pathlib import Path

import tkinter as tk
from tkinter import filedialog, messagebox, ttk


# --------------------------------------------------------------------------- #
# Backend invocation
# --------------------------------------------------------------------------- #

def _python_executable(project_dir: Path) -> Path:
    if sys.platform.startswith("win"):
        return project_dir / ".venv" / "Scripts" / "python.exe"
    return project_dir / ".venv" / "bin" / "python"


def _open_in_file_manager(path: Path) -> None:
    try:
        if sys.platform.startswith("win"):
            os.startfile(str(path))  # noqa: S606  (intended: open in Explorer)
        elif sys.platform == "darwin":
            subprocess.run(["open", str(path)], check=False)
        else:
            subprocess.run(["xdg-open", str(path)], check=False)
    except Exception:
        pass


# --------------------------------------------------------------------------- #
# Hover tooltip
# --------------------------------------------------------------------------- #

class Tooltip:
    """A small description box that appears when the pointer rests on a widget."""

    def __init__(self, widget: tk.Widget, text: str, delay_ms: int = 450) -> None:
        self.widget = widget
        self.text = text
        self.delay_ms = delay_ms
        self._after_id: str | None = None
        self._tip: tk.Toplevel | None = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, _event: object = None) -> None:
        self._cancel()
        self._after_id = self.widget.after(self.delay_ms, self._show)

    def _show(self) -> None:
        if self._tip is not None:
            return
        try:
            x = self.widget.winfo_rootx() + 16
            y = self.widget.winfo_rooty() + self.widget.winfo_height() + 8
        except tk.TclError:
            return
        self._tip = tk.Toplevel(self.widget)
        self._tip.wm_overrideredirect(True)
        self._tip.wm_geometry(f"+{x}+{y}")
        label = tk.Label(
            self._tip,
            text=self.text,
            justify="left",
            wraplength=340,
            background="#ffffe0",
            foreground="#222222",
            relief="solid",
            borderwidth=1,
            padx=8,
            pady=6,
            font=("Helvetica", 11),
        )
        label.pack()

    def _hide(self, _event: object = None) -> None:
        self._cancel()
        if self._tip is not None:
            self._tip.destroy()
            self._tip = None

    def _cancel(self) -> None:
        if self._after_id is not None:
            try:
                self.widget.after_cancel(self._after_id)
            except tk.TclError:
                pass
            self._after_id = None


def attach_tooltip(widget: tk.Widget, text: str) -> None:
    Tooltip(widget, text)


# --------------------------------------------------------------------------- #
# GUI
# --------------------------------------------------------------------------- #

PAD = 12

# Plain-language descriptions shown on hover. Written for a bench scientist,
# not a programmer: what the control does and what a safe default is.
TIPS = {
    "input": (
        "The folder that holds your microscopy images. Pick the plate or "
        "acquisition folder — the app looks inside its XY sub-folders for "
        "TIFF pairs using the selected target/aSMA and DAPI channels."
    ),
    "output": (
        "Where results are saved. Each run creates a new time-stamped folder "
        "here, so an earlier run is never overwritten."
    ),
    "run": "Start the analysis. Progress appears below; you can leave it running.",
    "stop": "Stop the current run. Partial results may be left in the output folder.",
    "target_channel": (
        "The channel measured as alpha smooth muscle actin intensity. For the original dataset this "
        "was CH2; change it only when a preview or acquisition notes show a different aSMA channel."
    ),
    "dapi_channel": (
        "The channel used to segment and count DAPI-positive nuclei. For the current datasets this "
        "is usually CH4, which appears as blue punctate nuclei in the overlay."
    ),
    "model": (
        "The Cellpose model used to find nuclei and regions. Leave this as "
        "cpsam_v2 unless a maintainer tells you to change it."
    ),
    "background_value": (
        "Subtracts a flat brightness level from the aSMA channel before measuring "
        "signal, to correct for camera/background haze. 0 = no correction (default "
        "and recommended). Only change this if a lab member specifically asked you to."
    ),
    "flow_threshold": (
        "Controls how strict the software is about cell shape. The default (0.4) "
        "works for most images. Raise it only if oddly-shaped or merged blobs are "
        "being called single cells; lower it only if real cells are being missed. "
        "When in doubt, leave as is."
    ),
    "cellprob_threshold": (
        "How sure the software must be before it counts a pixel as inside a cell. "
        "The default (0.0) works for most images. If it finds too many small false "
        "objects, raise it a little (e.g. 0.2); if it misses real objects, lower it "
        "a little. When in doubt, leave as is."
    ),
    "diameter": (
        "Rough object size in pixels. Leave blank to let the software estimate it "
        "automatically — recommended for most users."
    ),
    "max_images": (
        "Limit how many images per folder are processed. Leave blank to process all "
        "of them; set a small number (e.g. 2) for a quick test run first."
    ),
    "gpu": (
        "Speeds things up on computers with a compatible graphics card. If the "
        "progress log shows red text mentioning 'CUDA' or 'GPU' and the run stops, "
        "come back here, switch this off, and run again — it will just be slower."
    ),
    "figures": (
        "Also render QC overlay images so you can visually check the "
        "segmentation. Turn this off to finish a little faster."
    ),
}
DEFAULTS = {
    "model": "cpsam_v2",
    "gpu": True,
    "background_value": "0.0",
    "flow_threshold": "0.4",
    "cellprob_threshold": "0.0",
    "diameter": "",
    "max_images": "",
    "target_channel": "CH2",
    "dapi_channel": "CH4",
    "render_figures": True,
}


class PipelineGUI:
    def __init__(self, root: tk.Tk, project_dir: Path) -> None:
        self.root = root
        self.project_dir = project_dir
        self.python_bin = _python_executable(project_dir)
        self.runner = project_dir / "scripts" / "run_user_cellpose_batch.py"

        self.log_queue: "queue.Queue[str | tuple[str, object]]" = queue.Queue()
        self.worker: threading.Thread | None = None
        self.proc: subprocess.Popen[str] | None = None
        self.cancelled = False
        self.rows_processed: int | None = None
        self.run_output: Path | None = None
        self.advanced_visible = False

        self._build_style()
        self._build_widgets()
        self.root.after(100, self._drain_log_queue)

    # -- layout ------------------------------------------------------------- #

    def _build_style(self) -> None:
        self.root.title("Cellpose DAPI / aSMA Pipeline")
        self.root.minsize(720, 560)
        try:
            self.root.tk.call("tk", "scaling", 1.15)
        except tk.TclError:
            pass
        style = ttk.Style(self.root)
        if "aqua" not in style.theme_names() or sys.platform != "darwin":
            try:
                style.theme_use("clam")
            except tk.TclError:
                pass
        style.configure("Title.TLabel", font=("Helvetica", 18, "bold"))
        style.configure("Sub.TLabel", foreground="#555555")
        style.configure("Run.TButton", font=("Helvetica", 13, "bold"), padding=10)

    def _build_widgets(self) -> None:
        outer = ttk.Frame(self.root, padding=PAD)
        outer.pack(fill="both", expand=True)
        outer.columnconfigure(0, weight=1)

        # Header
        ttk.Label(outer, text="Cellpose DAPI / aSMA Pipeline", style="Title.TLabel").grid(
            row=0, column=0, sticky="w"
        )
        ttk.Label(
            outer,
            text="Pick your image folder and where to save results, then press Run.",
            style="Sub.TLabel",
        ).grid(row=1, column=0, sticky="w", pady=(2, PAD))

        # Folder pickers
        folders = ttk.LabelFrame(outer, text="Folders", padding=PAD)
        folders.grid(row=2, column=0, sticky="ew")
        folders.columnconfigure(1, weight=1)

        self.input_var = tk.StringVar()
        self.output_var = tk.StringVar()
        self.target_channel_var = tk.StringVar(value=DEFAULTS["target_channel"])
        self.dapi_channel_var = tk.StringVar(value=DEFAULTS["dapi_channel"])

        in_label = ttk.Label(folders, text="Image folder")
        in_label.grid(row=0, column=0, sticky="w", padx=(0, PAD))
        self.input_entry = ttk.Entry(folders, textvariable=self.input_var)
        self.input_entry.grid(row=0, column=1, sticky="ew")
        in_btn = ttk.Button(folders, text="Browse…", command=self._pick_input)
        in_btn.grid(row=0, column=2, padx=(PAD, 0))
        for widget in (in_label, self.input_entry, in_btn):
            attach_tooltip(widget, TIPS["input"])

        out_label = ttk.Label(folders, text="Results go to")
        out_label.grid(row=1, column=0, sticky="w", padx=(0, PAD), pady=(PAD, 0))
        self.output_entry = ttk.Entry(folders, textvariable=self.output_var)
        self.output_entry.grid(row=1, column=1, sticky="ew", pady=(PAD, 0))
        out_btn = ttk.Button(folders, text="Browse…", command=self._pick_output)
        out_btn.grid(row=1, column=2, padx=(PAD, 0), pady=(PAD, 0))
        for widget in (out_label, self.output_entry, out_btn):
            attach_tooltip(widget, TIPS["output"])

        target_label = ttk.Label(folders, text="aSMA target channel")
        target_label.grid(row=2, column=0, sticky="w", padx=(0, PAD), pady=(PAD, 0))
        target_box = ttk.Combobox(
            folders,
            textvariable=self.target_channel_var,
            values=("CH1", "CH2", "CH3", "CH4", "CH5"),
            width=8,
        )
        target_box.grid(row=2, column=1, sticky="w", pady=(PAD, 0))
        dapi_label = ttk.Label(folders, text="DAPI nuclei channel")
        dapi_label.grid(row=3, column=0, sticky="w", padx=(0, PAD), pady=(6, 0))
        dapi_box = ttk.Combobox(
            folders,
            textvariable=self.dapi_channel_var,
            values=("CH1", "CH2", "CH3", "CH4", "CH5"),
            width=8,
        )
        dapi_box.grid(row=3, column=1, sticky="w", pady=(6, 0))
        for widget in (target_label, target_box):
            attach_tooltip(widget, TIPS["target_channel"])
        for widget in (dapi_label, dapi_box):
            attach_tooltip(widget, TIPS["dapi_channel"])
        target_box.bind("<<ComboboxSelected>>", self._channel_mapping_changed, add="+")
        dapi_box.bind("<<ComboboxSelected>>", self._channel_mapping_changed, add="+")
        target_box.bind("<FocusOut>", self._channel_mapping_changed, add="+")
        dapi_box.bind("<FocusOut>", self._channel_mapping_changed, add="+")

        self.discovery_var = tk.StringVar(value="")
        ttk.Label(folders, textvariable=self.discovery_var, style="Sub.TLabel").grid(
            row=4, column=1, sticky="w", pady=(PAD, 0)
        )

        # Advanced toggle
        self.adv_button = ttk.Button(
            outer, text="▸ Advanced options", command=self._toggle_advanced, width=22
        )
        self.adv_button.grid(row=3, column=0, sticky="w", pady=(PAD, 0))

        self.adv_frame = ttk.LabelFrame(outer, text="Advanced options", padding=PAD)
        self.adv_frame.columnconfigure(1, weight=1)
        self._build_advanced(self.adv_frame)
        # not gridded until toggled

        # Run button + progress
        run_row = ttk.Frame(outer)
        run_row.grid(row=5, column=0, sticky="ew", pady=(PAD, 6))
        run_row.columnconfigure(0, weight=1)
        self.run_button = ttk.Button(
            run_row, text="Run analysis", style="Run.TButton", command=self._start_run
        )
        self.run_button.grid(row=0, column=0, sticky="ew")
        attach_tooltip(self.run_button, TIPS["run"])
        self.stop_button = ttk.Button(run_row, text="Stop", command=self._stop_run, width=8)
        self.stop_button.grid(row=0, column=1, padx=(8, 0))
        self.stop_button.state(["disabled"])
        attach_tooltip(self.stop_button, TIPS["stop"])

        self.progress = ttk.Progressbar(outer, mode="indeterminate")
        self.progress.grid(row=6, column=0, sticky="ew")

        self.status_var = tk.StringVar(value="Ready.")
        ttk.Label(outer, textvariable=self.status_var, style="Sub.TLabel").grid(
            row=7, column=0, sticky="w", pady=(4, 6)
        )

        # Log
        logframe = ttk.LabelFrame(outer, text="Progress log", padding=6)
        logframe.grid(row=8, column=0, sticky="nsew")
        outer.rowconfigure(8, weight=1)
        logframe.rowconfigure(0, weight=1)
        logframe.columnconfigure(0, weight=1)
        self.log = tk.Text(logframe, height=10, wrap="word", state="disabled",
                           background="#1e1e1e", foreground="#e0e0e0", relief="flat")
        self.log.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(logframe, command=self.log.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.log["yscrollcommand"] = scroll.set

        self._sync_run_enabled()

    def _build_advanced(self, frame: ttk.Frame) -> None:
        self.model_var = tk.StringVar(value=DEFAULTS["model"])
        self.gpu_var = tk.BooleanVar(value=DEFAULTS["gpu"])
        self.bg_var = tk.StringVar(value=DEFAULTS["background_value"])
        self.flow_var = tk.StringVar(value=DEFAULTS["flow_threshold"])
        self.cellprob_var = tk.StringVar(value=DEFAULTS["cellprob_threshold"])
        self.diameter_var = tk.StringVar(value=DEFAULTS["diameter"])
        self.maximg_var = tk.StringVar(value=DEFAULTS["max_images"])
        self.figures_var = tk.BooleanVar(value=DEFAULTS["render_figures"])

        ttk.Label(
            frame,
            text="Most users can leave these as-is. Hover any option (ⓘ) for a description.",
            style="Sub.TLabel",
        ).grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, PAD))

        rows = [
            ("Model name or path", self.model_var, "model"),
            ("CH2 background value", self.bg_var, "background_value"),
            ("Flow threshold", self.flow_var, "flow_threshold"),
            ("Cell probability threshold", self.cellprob_var, "cellprob_threshold"),
            ("Diameter (blank = auto)", self.diameter_var, "diameter"),
            ("Max images per acquisition", self.maximg_var, "max_images"),
        ]
        r = 1
        for label_text, var, key in rows:
            label = ttk.Label(frame, text=label_text)
            label.grid(row=r, column=0, sticky="w", pady=4, padx=(0, PAD))
            entry = ttk.Entry(frame, textvariable=var, width=24)
            entry.grid(row=r, column=1, sticky="w", pady=4)
            info = ttk.Label(frame, text="ⓘ", foreground="#3573b9", cursor="question_arrow")
            info.grid(row=r, column=2, sticky="w", padx=(PAD, 0))
            for widget in (label, entry, info):
                attach_tooltip(widget, TIPS[key])
            r += 1

        gpu_cb = ttk.Checkbutton(
            frame, text="Use GPU acceleration when available  ⓘ", variable=self.gpu_var
        )
        gpu_cb.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
        attach_tooltip(gpu_cb, TIPS["gpu"])
        r += 1
        fig_cb = ttk.Checkbutton(
            frame, text="Render QC overlay figures  ⓘ", variable=self.figures_var
        )
        fig_cb.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
        attach_tooltip(fig_cb, TIPS["figures"])
        r += 1
        ttk.Button(frame, text="Reset to defaults", command=self._reset_advanced).grid(
            row=r, column=0, sticky="w", pady=(PAD, 0)
        )

    # -- behaviour ---------------------------------------------------------- #

    def _toggle_advanced(self) -> None:
        self.advanced_visible = not self.advanced_visible
        if self.advanced_visible:
            self.adv_frame.grid(row=4, column=0, sticky="ew", pady=(6, 0))
            self.adv_button.config(text="▾ Advanced options")
        else:
            self.adv_frame.grid_forget()
            self.adv_button.config(text="▸ Advanced options")

    def _reset_advanced(self) -> None:
        self.model_var.set(DEFAULTS["model"])
        self.gpu_var.set(DEFAULTS["gpu"])
        self.bg_var.set(DEFAULTS["background_value"])
        self.flow_var.set(DEFAULTS["flow_threshold"])
        self.cellprob_var.set(DEFAULTS["cellprob_threshold"])
        self.diameter_var.set(DEFAULTS["diameter"])
        self.maximg_var.set(DEFAULTS["max_images"])
        self.figures_var.set(DEFAULTS["render_figures"])

    def _channel_mapping_changed(self, _event: object = None) -> None:
        input_path = self.input_var.get().strip()
        if input_path and Path(input_path).is_dir():
            self._run_discovery(Path(input_path))

    def _pick_input(self) -> None:
        path = filedialog.askdirectory(title="Select the folder that contains the image data")
        if path:
            self.input_var.set(path)
            self._sync_run_enabled()
            self._run_discovery(Path(path))

    def _pick_output(self) -> None:
        path = filedialog.askdirectory(title="Select where the results folder should be created")
        if path:
            self.output_var.set(path)
            self._sync_run_enabled()

    def _sync_run_enabled(self) -> None:
        ready = bool(self.input_var.get()) and bool(self.output_var.get())
        self.run_button.state(["!disabled"] if (ready and self.worker is None) else ["disabled"])

    # -- discovery preview -------------------------------------------------- #

    def _run_discovery(self, input_root: Path) -> None:
        self.discovery_var.set("Scanning folder…")
        target_channel = self.target_channel_var.get().strip().upper() or DEFAULTS["target_channel"]
        dapi_channel = self.dapi_channel_var.get().strip().upper() or DEFAULTS["dapi_channel"]

        def work() -> None:
            try:
                sys.path.insert(0, str(self.project_dir / "src"))
                from dapi_norm.user_cellpose_batch import discover_acquisitions

                acqs = discover_acquisitions(
                    input_root,
                    target_channel_id=target_channel,
                    dapi_channel_id=dapi_channel,
                )
                pairs = sum(a.image_count for a in acqs)
                if acqs:
                    msg = (
                        f"Found {len(acqs)} acquisition folder(s), {pairs} "
                        f"{target_channel}/{dapi_channel} image pair(s)."
                    )
                else:
                    msg = f"No {target_channel}/{dapi_channel} image pairs found in that folder."
            except Exception as exc:  # noqa: BLE001
                msg = f"Could not scan folder: {exc}"
            self.log_queue.put(("discovery", msg))

        threading.Thread(target=work, daemon=True).start()

    # -- run ---------------------------------------------------------------- #

    def _build_command(self) -> list[str] | None:
        input_path = self.input_var.get().strip()
        output_path = self.output_var.get().strip()
        if not Path(input_path).is_dir():
            messagebox.showerror(
                "Folder not found",
                f"The image folder does not exist:\n{input_path}\n\nPick it again with Browse.",
            )
            return None
        if not Path(output_path).is_dir():
            messagebox.showerror(
                "Folder not found",
                f"The results folder does not exist:\n{output_path}\n\nPick it again with Browse.",
            )
            return None
        target_channel = self._validated_channel(self.target_channel_var.get(), "aSMA target channel")
        if target_channel is None:
            return None
        dapi_channel = self._validated_channel(self.dapi_channel_var.get(), "DAPI nuclei channel")
        if dapi_channel is None:
            return None
        if target_channel == dapi_channel:
            messagebox.showerror(
                "Invalid channel mapping",
                "The aSMA target channel and DAPI nuclei channel must be different.",
            )
            return None

        stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        self.run_output = Path(output_path) / f"Cellpose_DAPI_aSMA_{stamp}"

        cmd = [
            str(self.python_bin),
            "-u",
            str(self.runner),
            "--input", input_path,
            "--output", str(self.run_output),
            "--model", self.model_var.get().strip() or DEFAULTS["model"],
            "--target-channel", target_channel,
            "--dapi-channel", dapi_channel,
            "--gpu" if self.gpu_var.get() else "--cpu",
            "--render-figures" if self.figures_var.get() else "--skip-figures",
        ]

        def add_float(flag: str, raw: str, name: str) -> bool:
            raw = raw.strip()
            if raw == "":
                return True
            try:
                value = float(raw)
            except ValueError:
                messagebox.showerror("Invalid value", f"{name} must be a number, got '{raw}'.")
                return False
            if not math.isfinite(value):
                messagebox.showerror(
                    "Invalid value",
                    f"{name} must be an ordinary number, not '{raw}'.",
                )
                return False
            cmd.extend([flag, raw])
            return True

        if not add_float("--background-value", self.bg_var.get(), "CH2 background value"):
            return None
        if not add_float("--flow-threshold", self.flow_var.get(), "Flow threshold"):
            return None
        if not add_float("--cellprob-threshold", self.cellprob_var.get(),
                         "Cell probability threshold"):
            return None
        if self.diameter_var.get().strip():
            if not add_float("--diameter", self.diameter_var.get(), "Diameter"):
                return None
        if self.maximg_var.get().strip():
            raw = self.maximg_var.get().strip()
            if not raw.isdigit() or int(raw) < 1:
                messagebox.showerror(
                    "Invalid value",
                    "Max images per acquisition must be a whole number of 1 or more "
                    f"(leave blank to process all), got '{raw}'.",
                )
                return None
            cmd.extend(["--max-images-per-acquisition", raw])
        return cmd

    def _validated_channel(self, raw_value: str, field_name: str) -> str | None:
        value = raw_value.strip().upper()
        if not re.fullmatch(r"CH\d+", value):
            messagebox.showerror(
                "Invalid channel",
                f"{field_name} must look like CH1, CH2, CH4, etc.; got '{raw_value}'.",
            )
            return None
        return value

    def _start_run(self) -> None:
        if self.worker is not None:
            return
        if not self.python_bin.exists():
            messagebox.showerror("Setup needed",
                                 "The analysis environment is missing. Run the Setup launcher first.")
            return
        cmd = self._build_command()
        if cmd is None:
            return

        self._log_clear()
        self._log_line(f"Input:   {self.input_var.get()}")
        self._log_line(f"Output:  {self.run_output}")
        self._log_line("Starting Cellpose pipeline. Full plates can take a while.\n")

        self.cancelled = False
        self.rows_processed = None
        self.run_button.state(["disabled"])
        self.stop_button.state(["!disabled"])
        self.progress.start(12)
        self.status_var.set("Running…")

        env = dict(os.environ)
        env["CELLPOSE_LOCAL_MODELS_PATH"] = str(self.project_dir / ".models" / "cellpose")

        def work() -> None:
            code = 1
            try:
                proc = subprocess.Popen(
                    cmd,
                    cwd=str(self.project_dir),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                    env=env,
                )
                self.proc = proc
                assert proc.stdout is not None
                for line in proc.stdout:
                    stripped = line.rstrip("\n")
                    if "Processed image rows:" in stripped:
                        try:
                            self.rows_processed = int(stripped.split(":")[-1].strip())
                        except ValueError:
                            pass
                    self.log_queue.put(stripped)
                code = proc.wait()
            except Exception as exc:  # noqa: BLE001
                self.log_queue.put(f"ERROR: {exc}")
            finally:
                self.proc = None
            self.log_queue.put(("done", code))

        self.worker = threading.Thread(target=work, daemon=True)
        self.worker.start()

    def _stop_run(self) -> None:
        proc = self.proc
        if proc is None:
            return
        self.cancelled = True
        self.status_var.set("Stopping…")
        self._log_line("\nStopping — asking the analysis to shut down…")
        try:
            proc.terminate()
        except Exception:  # noqa: BLE001
            pass

    def _finish_run(self, code: int) -> None:
        self.progress.stop()
        self.worker = None
        self.stop_button.state(["disabled"])
        self._sync_run_enabled()

        if self.cancelled:
            self.status_var.set("Stopped.")
            self._log_line("\nStopped by user. Any partial results are in the output folder.")
            messagebox.showinfo("Stopped", "The run was stopped. No complete results were produced.")
            return

        if code != 0 or self.run_output is None:
            self.status_var.set("Failed. See the progress log.")
            messagebox.showerror(
                "Analysis failed",
                "The pipeline did not finish. The progress log shows what went wrong — "
                "copy that text if you need to ask for help.",
            )
            return

        if self.rows_processed == 0:
            self.status_var.set("Finished, but no images were processed.")
            messagebox.showwarning(
                "Nothing to process",
                "The run finished but processed 0 images. The image folder probably does not "
                "contain the selected channel TIFF pairs inside XY sub-folders. Check the folder "
                "and try again.",
            )
            return

        final = self.run_output / "final"
        summary = final / "START_HERE_RUN_SUMMARY.html"
        self.status_var.set("Done. Results opened in a new window.")
        self._log_line(f"\nDone. Results: {final}")
        if final.exists():
            _open_in_file_manager(final)
        if summary.exists():
            _open_in_file_manager(summary)
        messagebox.showinfo(
            "Finished",
            f"Analysis finished.\n\nResults were written to:\n{final}",
        )

    # -- log plumbing ------------------------------------------------------- #

    def _drain_log_queue(self) -> None:
        try:
            while True:
                item = self.log_queue.get_nowait()
                if isinstance(item, tuple):
                    kind, payload = item
                    if kind == "done":
                        self._finish_run(int(payload))
                    elif kind == "discovery":
                        self.discovery_var.set(str(payload))
                else:
                    self._log_line(item)
        except queue.Empty:
            pass
        self.root.after(100, self._drain_log_queue)

    def _log_clear(self) -> None:
        self.log.config(state="normal")
        self.log.delete("1.0", "end")
        self.log.config(state="disabled")

    def _log_line(self, text: str) -> None:
        self.log.config(state="normal")
        self.log.insert("end", text + "\n")
        self.log.see("end")
        self.log.config(state="disabled")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    args = parser.parse_args()
    project_dir = Path(args.project_dir).resolve()

    root = tk.Tk()
    PipelineGUI(root, project_dir)
    # Bring the window to the front on launch. When started from a double-click
    # launcher (not an .app bundle), macOS otherwise leaves it behind other
    # windows with no Dock icon to click.
    root.update_idletasks()
    root.deiconify()
    root.lift()
    root.attributes("-topmost", True)
    root.after(800, lambda: root.attributes("-topmost", False))
    try:
        root.focus_force()
    except Exception:  # noqa: BLE001
        pass
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
