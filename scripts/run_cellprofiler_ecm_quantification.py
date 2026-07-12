from __future__ import annotations

import argparse
import csv
import os
from collections import Counter
from pathlib import Path
import shutil
import subprocess
import textwrap

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import tifffile
import yaml

from dapi_norm.ecm_cellprofiler import (
    ECM_ENDPOINT,
    EcmImageRecord,
    EcmMeasurement,
    discover_ecm_records,
    find_cellprofiler_mask,
    find_cellprofiler_overlay,
    load_dapi_counts,
    mask_boundary_rgba,
    measure_ecm_from_mask,
    normalize_for_display,
    robust_background_threshold_preview,
    select_representative_records,
    stage_ecm_images,
    write_csv,
    write_ecm_workbook,
)
from dapi_norm.image_arrays import read_primary_intensity_plane


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_ROOT = (
    REPO_ROOT
    / "data"
    / "new_download_2026_07_04_partial_test"
    / "CFB1 a-SMA ECM for 6 drugs"
)
DEFAULT_CELLPOSE_REPORT_ROOT = REPO_ROOT / "reports" / "new_download_2026_07_04_cellpose_CH2_CH4"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "reports" / "cellprofiler_ecm_ch1_2026_07_06_k2_arm64"
DEFAULT_CELLPROFILER_APP = Path.home() / "Applications" / "CellProfiler.app"
DEFAULT_CELLPROFILER_ARM64_PREFIX = REPO_ROOT / ".cellprofiler-arm64"
DEFAULT_CELLPROFILER_ARM64_WRAPPER = REPO_ROOT / "scripts" / "cellprofiler_arm64.sh"
DEFAULT_SELECTED_DEVIATIONS = 2.0


def main() -> None:
    args = parse_args()
    output_root = args.output_root.resolve()
    if args.overwrite and output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    logs_dir = output_root / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    records = discover_ecm_records(
        args.input_root.resolve(),
        ecm_channel_id=args.ecm_channel,
        dapi_channel_id=args.dapi_channel,
    )
    if args.max_images_per_acquisition is not None:
        records = limit_per_acquisition(records, args.max_images_per_acquisition)
    dapi_counts = load_dapi_counts(args.cellpose_report_root.resolve())

    representative_records = select_representative_records(
        records,
        representative_count=args.representative_count,
    )
    write_manifest(logs_dir / "resolved_ecm_input_manifest.csv", records)
    write_manifest(logs_dir / "representative_field_manifest.csv", representative_records)
    (logs_dir / "cellprofiler_version.txt").write_text(
        describe_cellprofiler(
            cellprofiler_app=args.cellprofiler_app,
            cellprofiler_executable=args.cellprofiler_executable,
        ),
        encoding="utf-8",
    )

    candidate_summaries = run_parameter_sweep(
        args=args,
        output_root=output_root,
        records=representative_records,
        dapi_counts=dapi_counts,
    )
    candidate_summary_path = output_root / "tables" / "ecm_parameter_sweep_summary.csv"
    write_csv(candidate_summary_path, candidate_summaries)
    candidate_aggregate_path = output_root / "tables" / "ecm_parameter_sweep_aggregate.csv"
    write_csv(candidate_aggregate_path, aggregate_candidate_rows(candidate_summaries))
    selected_override = None if args.auto_select_deviations else args.selected_deviations
    selected_deviations = resolve_selected_deviations(
        candidate_summaries,
        selected_deviations=selected_override,
    )
    selection_source = "automatic_heuristic" if args.auto_select_deviations else "human_qc_override"
    selected_path = logs_dir / "selected_ecm_threshold.yaml"
    selected_path.write_text(
        yaml.safe_dump(
            {
                "selected_method": "cellprofiler_robust_background",
                "selected_deviations": selected_deviations,
                "selection_source": selection_source,
                "selection_rule": (
                    "Automatic heuristic prefers candidates without near-full masks and with "
                    "low QC flag burden. Human QC override is used when visual review shows "
                    "the automatic candidate is too strict or too permissive."
                ),
                "candidate_summary_csv": str(candidate_summary_path.relative_to(REPO_ROOT)),
                "candidate_aggregate_csv": str(candidate_aggregate_path.relative_to(REPO_ROOT)),
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    if args.sweep_only:
        print(f"Processed ECM representative rows: {len(representative_records)}")
        print(f"Selected CellProfiler Robust Background deviations: {selected_deviations:g}")
        print(f"Parameter sweep CSV: {candidate_summary_path}")
        print(f"Parameter sweep aggregate CSV: {candidate_aggregate_path}")
        print(
            "Parameter sweep panel: "
            f"{output_root / 'qc' / 'parameter_selection' / 'ecm_parameter_sweep_panel.png'}"
        )
        return

    final_measurements = run_final_cellprofiler(
        args=args,
        output_root=output_root,
        records=records,
        dapi_counts=dapi_counts,
        selected_deviations=selected_deviations,
    )
    write_final_outputs(
        output_root=output_root,
        measurements=final_measurements,
        selected_deviations=selected_deviations,
        ecm_channel_id=args.ecm_channel,
        dapi_channel_id=args.dapi_channel,
    )
    render_representative_review_panel(
        output_root=output_root,
        records=representative_records,
        measurement_lookup=measurement_lookup(final_measurements),
        output_path=output_root / "qc" / "review" / "ecm_representative_review_panel.png",
    )
    render_acquisition_qc_pages(
        output_root=output_root,
        records=records,
        measurement_lookup=measurement_lookup(final_measurements),
    )
    write_run_summary(
        output_root=output_root,
        input_root=args.input_root.resolve(),
        records=records,
        measurements=final_measurements,
        selected_deviations=selected_deviations,
        candidate_summaries=candidate_summaries,
    )

    print(f"Processed ECM image rows: {len(final_measurements)}")
    print(f"Selected CellProfiler Robust Background deviations: {selected_deviations:g}")
    print(f"Summary CSV: {output_root / 'final' / 'tables' / 'cellprofiler_ecm_ch1_summary.csv'}")
    print(f"Workbook: {output_root / 'final' / 'workbooks' / 'cellprofiler_ecm_ch1_summary.xlsx'}")
    print(f"QC review panel: {output_root / 'qc' / 'review' / 'ecm_representative_review_panel.png'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run actual CellProfiler ECM CH1 quantification and QC panels."
    )
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--cellpose-report-root", type=Path, default=DEFAULT_CELLPOSE_REPORT_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--cellprofiler-app", type=Path, default=DEFAULT_CELLPROFILER_APP)
    parser.add_argument(
        "--cellprofiler-executable",
        type=Path,
        default=default_cellprofiler_executable(),
        help=(
            "Executable or wrapper to run CellProfiler. If provided, this overrides "
            "--cellprofiler-app and is useful for project-local native Apple Silicon installs. "
            "Defaults to scripts/cellprofiler_arm64.sh on non-Windows systems when the "
            "project-local .cellprofiler-arm64 environment exists."
        ),
    )
    parser.add_argument("--ecm-channel", default="CH1")
    parser.add_argument("--dapi-channel", default="CH4")
    parser.add_argument("--candidate-deviations", nargs="+", type=float, default=[2.0, 3.0, 5.0, 8.0])
    parser.add_argument(
        "--selected-deviations",
        type=float,
        default=DEFAULT_SELECTED_DEVIATIONS,
        help=(
            "Force the final CellProfiler Robust Background '# of deviations' value. "
            "The value must be present in --candidate-deviations so it is included in the QC sweep. "
            f"Default: {DEFAULT_SELECTED_DEVIATIONS:g}, the current human-QC-corrected ECM setting."
        ),
    )
    parser.add_argument(
        "--auto-select-deviations",
        action="store_true",
        help=(
            "Use the legacy automatic heuristic instead of the current human-QC-corrected "
            "--selected-deviations value."
        ),
    )
    parser.add_argument("--representative-count", type=int, default=24)
    parser.add_argument("--max-images-per-acquisition", type=int)
    parser.add_argument("--skip-cellprofiler", action="store_true")
    parser.add_argument(
        "--sweep-only",
        action="store_true",
        help="Run only the representative CellProfiler threshold sweep and QC panels.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def default_cellprofiler_executable() -> Path | None:
    if os.name == "nt":
        return None
    if DEFAULT_CELLPROFILER_ARM64_PREFIX.exists() and DEFAULT_CELLPROFILER_ARM64_WRAPPER.exists():
        return DEFAULT_CELLPROFILER_ARM64_WRAPPER
    return None


def limit_per_acquisition(records: list[EcmImageRecord], max_images: int) -> list[EcmImageRecord]:
    kept: list[EcmImageRecord] = []
    counts: Counter[str] = Counter()
    for record in records:
        if counts[record.acquisition_rel] >= max_images:
            continue
        kept.append(record)
        counts[record.acquisition_rel] += 1
    return kept


def write_manifest(path: Path, records: list[EcmImageRecord]) -> None:
    rows = [
        {
            "acquisition": record.acquisition_name,
            "acquisition_rel": record.acquisition_rel,
            "location": record.location,
            "source_id": record.source_id,
            "ecm_path": str(record.ecm_path),
            "dapi_path": str(record.dapi_path),
            "staged_name": record.staged_name,
        }
        for record in records
    ]
    write_csv(path, rows)


def run_parameter_sweep(
    *,
    args: argparse.Namespace,
    output_root: Path,
    records: list[EcmImageRecord],
    dapi_counts: dict[tuple[str, str], int],
) -> list[dict[str, object]]:
    all_rows: list[dict[str, object]] = []
    measurements_by_candidate: dict[float, list[EcmMeasurement]] = {}
    for deviations in args.candidate_deviations:
        candidate_root = output_root / "parameter_selection" / f"robust_bg_dev_{deviations:g}"
        run_root = run_cellprofiler_for_records(
            cellprofiler_app=args.cellprofiler_app,
            cellprofiler_executable=args.cellprofiler_executable,
            output_root=candidate_root,
            records=records,
            deviations=deviations,
            skip_cellprofiler=args.skip_cellprofiler,
        )
        measurements = measure_records(
            records=records,
            dapi_counts=dapi_counts,
            run_root=run_root,
            threshold_deviations=deviations,
            ecm_channel_id=args.ecm_channel,
            dapi_channel_id=args.dapi_channel,
            qc_panel_root=output_root,
        )
        measurements_by_candidate[deviations] = measurements
        for measurement in measurements:
            row = measurement.as_row()
            row["candidate_deviations"] = deviations
            all_rows.append(row)
    render_parameter_selection_panel(
        output_root=output_root,
        records=records,
        measurements_by_candidate=measurements_by_candidate,
        output_path=output_root / "qc" / "parameter_selection" / "ecm_parameter_sweep_panel.png",
    )
    render_parameter_selection_summary(
        measurements_by_candidate=measurements_by_candidate,
        output_path=output_root / "qc" / "parameter_selection" / "ecm_parameter_sweep_summary.png",
    )
    return all_rows


def run_final_cellprofiler(
    *,
    args: argparse.Namespace,
    output_root: Path,
    records: list[EcmImageRecord],
    dapi_counts: dict[tuple[str, str], int],
    selected_deviations: float,
) -> list[EcmMeasurement]:
    final_cp_root = run_cellprofiler_for_records(
        cellprofiler_app=args.cellprofiler_app,
        cellprofiler_executable=args.cellprofiler_executable,
        output_root=output_root / "cellprofiler_final",
        records=records,
        deviations=selected_deviations,
        skip_cellprofiler=args.skip_cellprofiler,
    )
    return measure_records(
        records=records,
        dapi_counts=dapi_counts,
        run_root=final_cp_root,
        threshold_deviations=selected_deviations,
        ecm_channel_id=args.ecm_channel,
        dapi_channel_id=args.dapi_channel,
        qc_panel_root=output_root,
    )


def run_cellprofiler_for_records(
    *,
    cellprofiler_app: Path,
    cellprofiler_executable: Path | None,
    output_root: Path,
    records: list[EcmImageRecord],
    deviations: float,
    skip_cellprofiler: bool,
) -> Path:
    staging_dir = output_root / "staging_grayscale"
    pipeline_dir = output_root / "pipelines"
    cp_run_dir = output_root / "cp_run"
    logs_dir = output_root / "logs"
    for path in [staging_dir, pipeline_dir, cp_run_dir, logs_dir]:
        path.mkdir(parents=True, exist_ok=True)
    stage_ecm_images(records, staging_dir, logs_dir / "resolved_input_manifest.csv")
    pipeline_path = pipeline_dir / f"ecm_ch1_robust_bg_dev_{deviations:g}.cppipe"
    pipeline_path.write_text(build_ecm_cellprofiler_pipeline(deviations=deviations), encoding="utf-8")
    if not skip_cellprofiler:
        if cp_run_dir.exists():
            shutil.rmtree(cp_run_dir)
        cp_run_dir.mkdir(parents=True, exist_ok=True)
        run_cellprofiler(
            cellprofiler_app=cellprofiler_app,
            cellprofiler_executable=cellprofiler_executable,
            pipeline_path=pipeline_path,
            input_dir=staging_dir,
            output_dir=cp_run_dir,
            log_path=logs_dir / "cellprofiler_run_log.txt",
        )
    copied_root = output_root / "cellprofiler_outputs"
    copy_cellprofiler_outputs(cp_run_dir=cp_run_dir, output_root=copied_root)
    return copied_root


def build_ecm_cellprofiler_pipeline(*, deviations: float) -> str:
    return textwrap.dedent(
        f"""\
        CellProfiler Pipeline: http://www.cellprofiler.org
        Version:5
        DateRevision:428
        GitHash:
        ModuleCount:10
        HasImagePlaneDetails:False

        Images:[module_num:1|svn_version:'Unknown'|variable_revision_number:2|show_window:False|notes:['Collect staged grayscale CH1 ECM TIFFs.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            :
            Filter images?:Images only
            Select the rule criteria:and (extension does isimage) (directory doesnot containregexp "[\\\\\\\\/]\\\\.")

        Metadata:[module_num:2|svn_version:'Unknown'|variable_revision_number:6|show_window:False|notes:['Extract acquisition and XY field from staged ECM file names.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Extract metadata?:Yes
            Metadata data type:Text
            Metadata types:{{}}
            Extraction method count:1
            Metadata extraction method:Extract from file/folder names
            Metadata source:File name
            Regular expression to extract from file name:^(?P<Acquisition>.+)__(?P<Field>XY[0-9]+)_ECM\\.tif$
            Regular expression to extract from folder name:(?P<Date>[0-9]{{4}}_[0-9]{{2}}_[0-9]{{2}})$
            Extract metadata from:All images
            Select the filtering criteria:and (file does contain "")
            Metadata file location:Elsewhere...|
            Match file and image metadata:[]
            Use case insensitive matching?:No
            Metadata file name:
            Does cached metadata exist?:No

        NamesAndTypes:[module_num:3|svn_version:'Unknown'|variable_revision_number:8|show_window:False|notes:['Assign CH1 ECM grayscale channel.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Assign a name to:Images matching rules
            Select the image type:Grayscale image
            Name to assign these images:ECM
            Match metadata:[]
            Image set matching method:Order
            Set intensity range from:Image metadata
            Assignments count:1
            Single images count:0
            Maximum intensity:65535.0
            Process as 3D?:No
            Relative pixel spacing in X:1.0
            Relative pixel spacing in Y:1.0
            Relative pixel spacing in Z:1.0
            Select the rule criteria:and (file does contain "_ECM.tif")
            Name to assign these images:ECM
            Name to assign these objects:ECMPositive
            Select the image type:Grayscale image
            Set intensity range from:Image metadata
            Maximum intensity:65535.0

        Groups:[module_num:4|svn_version:'Unknown'|variable_revision_number:2|show_window:False|notes:['No grouping.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Do you want to group your images?:No
            grouping metadata count:1
            Metadata category:None

        IdentifyPrimaryObjects:[module_num:5|svn_version:'Unknown'|variable_revision_number:14|show_window:True|notes:['Identify ECM-positive pixels using actual CellProfiler Robust Background thresholding.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the input image:ECM
            Name the primary objects to be identified:ECMPositive
            Typical diameter of objects, in pixel units (Min,Max):1,5000
            Discard objects outside the diameter range?:No
            Discard objects touching the border of the image?:No
            Method to distinguish clumped objects:Intensity
            Method to draw dividing lines between clumped objects:Intensity
            Size of smoothing filter:0
            Suppress local maxima that are closer than this minimum allowed distance:1.0
            Speed up by using lower-resolution image to find local maxima?:No
            Fill holes in identified objects?:Never
            Automatically calculate size of smoothing filter for declumping?:No
            Automatically calculate minimum allowed distance between local maxima?:No
            Handling of objects if excessive number of objects identified:Continue
            Maximum number of objects:100000
            Display accepted local maxima?:No
            Select maxima color:Blue
            Use advanced settings?:Yes
            Threshold setting version:11
            Threshold strategy:Global
            Thresholding method:Robust Background
            Threshold smoothing scale:0.0
            Threshold correction factor:1.0
            Lower and upper bounds on threshold:0.0,1.0
            Manual threshold:0.0
            Select the measurement to threshold with:None
            Two-class or three-class thresholding?:Two classes
            Assign pixels in the middle intensity class to the foreground or the background?:Foreground
            Size of adaptive window:50
            Lower outlier fraction:0.05
            Upper outlier fraction:0.05
            Averaging method:Mean
            Variance method:Standard deviation
            # of deviations:{deviations:g}
            Thresholding method:Robust Background

        ConvertObjectsToImage:[module_num:6|svn_version:'Unknown'|variable_revision_number:1|show_window:True|notes:['Convert ECM-positive CellProfiler objects to uint16 label image.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the input objects:ECMPositive
            Name the output image:ECMPositiveImage
            Select the color format:uint16
            Select the colormap:Default

        SaveImages:[module_num:7|svn_version:'Unknown'|variable_revision_number:16|show_window:True|notes:['Save CellProfiler ECM-positive labels.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the type of image to save:Image
            Select the image to save:ECMPositiveImage
            Select method for constructing file names:From image filename
            Select image name for file prefix:ECM
            Enter single file name:OrigBlue
            Number of digits:4
            Append a suffix to the image file name?:Yes
            Text to append to the image name:_CellProfilerECMLabels
            Saved file format:tiff
            Output file location:Default Output Folder|
            Image bit depth:16-bit integer
            Overwrite existing files without warning?:Yes
            When to save:Every cycle
            Record the file and path information to the saved image?:No
            Create subfolders in the output folder?:No
            Base image folder:Elsewhere...|
            How to save the series:T (Time)
            Save with lossless compression?:No

        OverlayOutlines:[module_num:8|svn_version:'Unknown'|variable_revision_number:4|show_window:True|notes:['Overlay ECM-positive outlines on CH1 ECM.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Display outlines on a blank image?:No
            Select image on which to display outlines:ECM
            Name the output image:ECMOverlay
            Outline display mode:Color
            Select method to determine brightness of outlines:Max of image
            How to outline:Thick
            Select outline color:Green
            Select objects to display:ECMPositive

        SaveImages:[module_num:9|svn_version:'Unknown'|variable_revision_number:16|show_window:True|notes:['Save CellProfiler ECM overlay.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the type of image to save:Image
            Select the image to save:ECMOverlay
            Select method for constructing file names:From image filename
            Select image name for file prefix:ECM
            Enter single file name:OrigBlue
            Number of digits:4
            Append a suffix to the image file name?:Yes
            Text to append to the image name:_CellProfilerECMOverlay
            Saved file format:png
            Output file location:Default Output Folder|
            Image bit depth:8-bit integer
            Overwrite existing files without warning?:Yes
            When to save:Every cycle
            Record the file and path information to the saved image?:No
            Create subfolders in the output folder?:No
            Base image folder:Elsewhere...|
            How to save the series:T (Time)
            Save with lossless compression?:Yes

        ExportToSpreadsheet:[module_num:10|svn_version:'Unknown'|variable_revision_number:13|show_window:True|notes:['Export actual CellProfiler ECM measurements.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the column delimiter:Comma (",")
            Add image metadata columns to your object data file?:Yes
            Add image file and folder names to your object data file?:Yes
            Select the measurements to export:No
            Calculate the per-image mean values for object measurements?:Yes
            Calculate the per-image median values for object measurements?:No
            Calculate the per-image standard deviation values for object measurements?:No
            Output file location:Default Output Folder|
            Create a GenePattern GCT file?:No
            Select source of sample row name:Metadata
            Select the image to use as the identifier:None
            Select the metadata to use as the identifier:Field
            Export all measurement types?:Yes
            Press button to select measurements:
            Representation of Nan/Inf:NaN
            Add a prefix to file names?:No
            Filename prefix:CellProfiler_
            Overwrite existing files without warning?:Yes
            Data to export:Do not use
            Combine these object measurements with those of the previous object?:No
            File name:DATA.csv
            Use the object name for the file name?:Yes
        """
    )


def run_cellprofiler(
    *,
    cellprofiler_app: Path,
    cellprofiler_executable: Path | None,
    pipeline_path: Path,
    input_dir: Path,
    output_dir: Path,
    log_path: Path,
) -> None:
    cp_executable = resolve_cellprofiler_executable(
        cellprofiler_app=cellprofiler_app,
        cellprofiler_executable=cellprofiler_executable,
    )
    env = os.environ.copy()
    if cellprofiler_executable is None:
        java_home = cellprofiler_app / "Contents" / "Resources" / "Home"
        if not java_home.exists():
            raise FileNotFoundError(java_home)
        env["JAVA_HOME"] = str(java_home)
    cmd = [
        str(cp_executable),
        "-c",
        "-r",
        "-p",
        str(pipeline_path),
        "-i",
        str(input_dir),
        "-o",
        str(output_dir),
    ]
    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(" ".join(cmd) + "\n\n")
        log_handle.flush()
        subprocess.run(cmd, check=True, env=env, stdout=log_handle, stderr=subprocess.STDOUT)


def resolve_cellprofiler_executable(
    *,
    cellprofiler_app: Path,
    cellprofiler_executable: Path | None,
) -> Path:
    if cellprofiler_executable is not None:
        executable = cellprofiler_executable.expanduser().resolve()
    else:
        executable = cellprofiler_app.expanduser().resolve() / "Contents" / "MacOS" / "cp"
    if not executable.exists():
        raise FileNotFoundError(executable)
    return executable


def copy_cellprofiler_outputs(*, cp_run_dir: Path, output_root: Path) -> None:
    masks_dir = output_root / "masks"
    exports_dir = output_root / "exports"
    qc_dir = output_root / "qc"
    for path in [masks_dir, exports_dir, qc_dir]:
        path.mkdir(parents=True, exist_ok=True)
    for path in cp_run_dir.iterdir():
        if not path.is_file():
            continue
        lower = path.name.lower()
        if lower.endswith(".csv"):
            shutil.copy2(path, exports_dir / path.name)
        elif "labels" in lower and lower.endswith((".tif", ".tiff")):
            shutil.copy2(path, masks_dir / normalize_cp_output_name(path.name))
        elif "overlay" in lower and lower.endswith(".png"):
            shutil.copy2(path, qc_dir / normalize_cp_output_name(path.name))


def normalize_cp_output_name(name: str) -> str:
    return name


def measure_records(
    *,
    records: list[EcmImageRecord],
    dapi_counts: dict[tuple[str, str], int],
    run_root: Path,
    threshold_deviations: float,
    ecm_channel_id: str,
    dapi_channel_id: str,
    qc_panel_root: Path,
) -> list[EcmMeasurement]:
    masks_dir = run_root / "masks"
    qc_dir = run_root / "qc"
    rows: list[EcmMeasurement] = []
    for record in records:
        image, _ = read_primary_intensity_plane(record.ecm_path)
        mask_path = find_cellprofiler_mask(masks_dir, record)
        mask = np.asarray(tifffile.imread(mask_path))
        overlay_path = find_cellprofiler_overlay(qc_dir, record)
        rows.append(
            measure_ecm_from_mask(
                record=record,
                image=image,
                mask=mask,
                dapi_count=dapi_count_for(record, dapi_counts),
                threshold_deviations=threshold_deviations,
                ecm_channel_id=ecm_channel_id,
                dapi_channel_id=dapi_channel_id,
                mask_path=mask_path,
                overlay_path=overlay_path,
                qc_panel_path=qc_panel_root / "qc" / "pages",
                root_for_paths=REPO_ROOT,
            )
        )
    return rows


def choose_threshold_deviations(candidate_rows: list[dict[str, object]]) -> float:
    by_candidate: dict[float, list[dict[str, object]]] = {}
    for row in candidate_rows:
        by_candidate.setdefault(float(row["candidate_deviations"]), []).append(row)
    summaries: list[tuple[tuple[int, int, float, float], float]] = []
    for deviations, rows in by_candidate.items():
        fractions = np.array([float(row["ecm_positive_area_fraction"]) for row in rows])
        near_full = int(np.count_nonzero(fractions > 0.60))
        near_empty = int(np.count_nonzero(fractions < 0.001))
        median_fraction = float(np.median(fractions))
        flag_count = sum(1 for row in rows if str(row["qc_flags"]))
        # Prefer reviewable masks, then lower flag burden, then a moderate area fraction.
        score = (near_full, flag_count, abs(median_fraction - 0.05), near_empty)
        summaries.append((score, deviations))
    summaries.sort(key=lambda item: item[0])
    return summaries[0][1]


def resolve_selected_deviations(
    candidate_rows: list[dict[str, object]],
    *,
    selected_deviations: float | None,
) -> float:
    candidate_values = sorted({float(row["candidate_deviations"]) for row in candidate_rows})
    if selected_deviations is None:
        return choose_threshold_deviations(candidate_rows)
    selected = float(selected_deviations)
    if selected not in candidate_values:
        values = ", ".join(f"{value:g}" for value in candidate_values)
        raise ValueError(
            f"--selected-deviations {selected:g} must be one of the swept candidate values: {values}"
        )
    return selected


def aggregate_candidate_rows(candidate_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_candidate: dict[float, list[dict[str, object]]] = {}
    for row in candidate_rows:
        by_candidate.setdefault(float(row["candidate_deviations"]), []).append(row)
    rows: list[dict[str, object]] = []
    for deviations, items in sorted(by_candidate.items()):
        area = np.array([float(row["ecm_positive_area_fraction"]) for row in items], dtype=float)
        endpoints = np.array(
            [
                float(row[ECM_ENDPOINT])
                for row in items
                if row.get(ECM_ENDPOINT) is not None and str(row.get(ECM_ENDPOINT, "")).strip()
            ],
            dtype=float,
        )
        rows.append(
            {
                "candidate_deviations": deviations,
                "representative_field_count": len(items),
                "median_ecm_positive_area_fraction": float(np.median(area)),
                "min_ecm_positive_area_fraction": float(np.min(area)),
                "max_ecm_positive_area_fraction": float(np.max(area)),
                "near_empty_mask_count": int(np.count_nonzero(area < 0.001)),
                "near_full_field_mask_count": int(np.count_nonzero(area > 0.60)),
                "review_status_count": sum(1 for row in items if row["qc_status"] == "review"),
                "dapi_count_missing_or_zero_context_count": sum(
                    1
                    for row in items
                    if row.get("dapi_positive_nucleus_count") in (None, "")
                    or int(float(row["dapi_positive_nucleus_count"])) == 0
                ),
                "median_ecm_integrated_background_corrected": (
                    float(np.median(endpoints)) if endpoints.size else ""
                ),
                "q1_ecm_integrated_background_corrected": (
                    float(np.quantile(endpoints, 0.25)) if endpoints.size else ""
                ),
                "q3_ecm_integrated_background_corrected": (
                    float(np.quantile(endpoints, 0.75)) if endpoints.size else ""
                ),
            }
        )
    return rows


def write_final_outputs(
    *,
    output_root: Path,
    measurements: list[EcmMeasurement],
    selected_deviations: float,
    ecm_channel_id: str,
    dapi_channel_id: str,
) -> None:
    final_dir = output_root / "final"
    tables_dir = final_dir / "tables"
    workbook_dir = final_dir / "workbooks"
    rows = [measurement.as_row() for measurement in measurements]
    write_csv(tables_dir / "cellprofiler_ecm_ch1_summary.csv", rows)
    write_ecm_workbook(workbook_dir / "cellprofiler_ecm_ch1_summary.xlsx", rows)
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["acquisition"]), []).append(row)
    for acquisition, acquisition_rows in grouped.items():
        safe = acquisition.replace("/", "_").replace(" ", "_")
        write_csv(tables_dir / f"{safe}_cellprofiler_ecm_ch1_summary.csv", acquisition_rows)
    (output_root / "logs" / "final_parameters.yaml").write_text(
        yaml.safe_dump(
            {
                "ecm_channel_id": ecm_channel_id,
                "dapi_channel_id": dapi_channel_id,
                "cellprofiler_threshold_method": "Robust Background",
                "cellprofiler_threshold_deviations": selected_deviations,
                "background_value_per_px": f"median of mask-negative {ecm_channel_id} pixels per image",
                "primary_endpoint": ECM_ENDPOINT,
                "normalization_denominator": "none",
                "dapi_count_role": "context_only_not_used_for_ecm_normalization",
                "manual_ground_truth_validation": "not_available",
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )


def render_parameter_selection_panel(
    *,
    output_root: Path,
    records: list[EcmImageRecord],
    measurements_by_candidate: dict[float, list[EcmMeasurement]],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    candidates = sorted(measurements_by_candidate)
    rows = min(len(records), 16)
    records = records[:rows]
    fig, axes = plt.subplots(rows, len(candidates) + 1, figsize=(3.8 * (len(candidates) + 1), 3.3 * rows), squeeze=False)
    for row_idx, record in enumerate(records):
        image, _ = read_primary_intensity_plane(record.ecm_path)
        ax = axes[row_idx, 0]
        ax.imshow(normalize_for_display(image), cmap="gray")
        ax.set_title("Raw CH1 ECM" if row_idx == 0 else "")
        ax.set_ylabel(f"{record.acquisition_name}\n{record.location}", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for col_idx, deviations in enumerate(candidates, start=1):
            measurement = {
                (m.acquisition, m.location): m for m in measurements_by_candidate[deviations]
            }[(record.acquisition_name, record.location)]
            mask_path = REPO_ROOT / measurement.cellprofiler_mask_path
            mask = np.asarray(tifffile.imread(mask_path))
            ax = axes[row_idx, col_idx]
            ax.imshow(normalize_for_display(image), cmap="gray")
            ax.imshow(mask_boundary_rgba(mask, color=(0.0, 1.0, 0.25)))
            if row_idx == 0:
                ax.set_title(f"Robust BG\n{k_label(deviations)}", fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.text(
                0.02,
                0.98,
                f"area {measurement.ecm_positive_area_fraction:.1%}\n"
                f"bg {measurement.ecm_background_value_per_px:.0f}\n"
                f"ECM corr {format_metric(measurement.ecm_positive_integrated_background_corrected)}",
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=7,
                color="white",
                bbox={"facecolor": "black", "alpha": 0.62, "pad": 2, "edgecolor": "none"},
            )
    fig.suptitle("CellProfiler ECM Robust Background threshold sweep", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def render_parameter_selection_summary(
    *,
    measurements_by_candidate: dict[float, list[EcmMeasurement]],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    candidates = sorted(measurements_by_candidate)
    area_data = [
        [m.ecm_positive_area_fraction for m in measurements_by_candidate[candidate]]
        for candidate in candidates
    ]
    endpoint_data = [
        [
            m.ecm_positive_integrated_background_corrected
            for m in measurements_by_candidate[candidate]
        ]
        for candidate in candidates
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].boxplot(area_data, tick_labels=[k_label(k) for k in candidates], showfliers=True)
    axes[0].set_ylabel("ECM positive area fraction")
    axes[0].set_title("Foreground mask size")
    axes[1].boxplot(endpoint_data, tick_labels=[k_label(k) for k in candidates], showfliers=True)
    axes[1].set_ylabel("Background-corrected ECM integrated intensity")
    axes[1].set_title("ECM endpoint distribution")
    fig.suptitle("ECM parameter sweep diagnostics")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def render_representative_review_panel(
    *,
    output_root: Path,
    records: list[EcmImageRecord],
    measurement_lookup: dict[tuple[str, str], EcmMeasurement],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    records = records[:16]
    cols = ["DAPI context", "Raw CH1 ECM", "Background-corrected ECM", "CellProfiler ECM mask + value"]
    fig, axes = plt.subplots(len(records), len(cols), figsize=(18, 3.1 * len(records)), squeeze=False)
    for row_idx, record in enumerate(records):
        ecm, _ = read_primary_intensity_plane(record.ecm_path)
        dapi, _ = read_primary_intensity_plane(record.dapi_path)
        measurement = measurement_lookup[(record.acquisition_name, record.location)]
        mask = np.asarray(tifffile.imread(REPO_ROOT / measurement.cellprofiler_mask_path))
        corrected = np.clip(ecm.astype(float) - measurement.ecm_background_value_per_px, 0, None)
        images = [dapi, ecm, corrected, ecm]
        cmaps = ["Blues", "gray", "Greens", "gray"]
        for col_idx, title in enumerate(cols):
            ax = axes[row_idx, col_idx]
            ax.imshow(normalize_for_display(images[col_idx]), cmap=cmaps[col_idx])
            ax.set_xticks([])
            ax.set_yticks([])
            if row_idx == 0:
                ax.set_title(title, fontsize=11)
            if col_idx == 0:
                ax.set_ylabel(f"{record.acquisition_name}\n{record.location}", fontsize=8)
            if col_idx == 3:
                ax.imshow(mask_boundary_rgba(mask, color=(0.0, 1.0, 0.25)))
                ax.text(
                    0.02,
                    0.98,
                    metric_text(measurement),
                    transform=ax.transAxes,
                    va="top",
                    ha="left",
                    fontsize=7,
                    color="white",
                    bbox={"facecolor": "black", "alpha": 0.65, "pad": 2, "edgecolor": "none"},
                )
    fig.suptitle("CellProfiler ECM quantification review: image, correction, mask, endpoint", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def render_acquisition_qc_pages(
    *,
    output_root: Path,
    records: list[EcmImageRecord],
    measurement_lookup: dict[tuple[str, str], EcmMeasurement],
) -> None:
    pages_root = output_root / "qc" / "pages"
    pages_root.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, list[EcmImageRecord]] = {}
    for record in records:
        grouped.setdefault(record.acquisition_name, []).append(record)
    for acquisition, acquisition_records in grouped.items():
        acquisition_records = sorted(acquisition_records, key=lambda row: row.location)
        for page_idx in range(0, len(acquisition_records), 12):
            chunk = acquisition_records[page_idx : page_idx + 12]
            fig, axes = plt.subplots(4, 3, figsize=(14, 16), squeeze=False)
            for ax in axes.flat:
                ax.axis("off")
            for ax, record in zip(axes.flat, chunk):
                ecm, _ = read_primary_intensity_plane(record.ecm_path)
                measurement = measurement_lookup[(record.acquisition_name, record.location)]
                mask = np.asarray(tifffile.imread(REPO_ROOT / measurement.cellprofiler_mask_path))
                ax.imshow(normalize_for_display(ecm), cmap="gray")
                ax.imshow(mask_boundary_rgba(mask, color=(0.0, 1.0, 0.25)))
                ax.set_title(record.location, fontsize=9)
                ax.text(
                    0.02,
                    0.98,
                    metric_text(measurement),
                    transform=ax.transAxes,
                    va="top",
                    ha="left",
                    fontsize=7,
                    color="white",
                    bbox={"facecolor": "black", "alpha": 0.65, "pad": 2, "edgecolor": "none"},
                )
                ax.axis("off")
            safe_acq = acquisition.replace(" ", "_").replace("/", "_")
            page_number = page_idx // 12 + 1
            output_path = pages_root / f"{safe_acq}_ecm_qc_page_{page_number:03d}.png"
            fig.suptitle(f"{acquisition}: CellProfiler ECM mask overlay and values, page {page_number}", fontsize=14)
            fig.tight_layout(rect=(0, 0, 1, 0.98))
            fig.savefig(output_path, dpi=170)
            plt.close(fig)


def metric_text(measurement: EcmMeasurement) -> str:
    return (
        f"area {measurement.ecm_positive_area_fraction:.1%}\n"
        f"bg/px {measurement.ecm_background_value_per_px:.0f}\n"
        f"ECM corr {measurement.ecm_positive_integrated_background_corrected:.2e}"
    )


def write_run_summary(
    *,
    output_root: Path,
    input_root: Path,
    records: list[EcmImageRecord],
    measurements: list[EcmMeasurement],
    selected_deviations: float,
    candidate_summaries: list[dict[str, object]],
) -> None:
    final_dir = output_root / "final"
    status_counts = Counter(m.qc_status for m in measurements)
    by_acquisition: dict[str, list[EcmMeasurement]] = {}
    for measurement in measurements:
        by_acquisition.setdefault(measurement.acquisition, []).append(measurement)
    flag_counts: Counter[str] = Counter()
    for measurement in measurements:
        for flag in measurement.qc_flags.split(";"):
            if flag:
                flag_counts[flag] += 1
    review_rows = [measurement for measurement in measurements if measurement.qc_status != "pass"]
    mask_count = len(list((output_root / "cellprofiler_final" / "cellprofiler_outputs" / "masks").glob("*.tif*")))
    overlay_count = len(list((output_root / "cellprofiler_final" / "cellprofiler_outputs" / "qc").glob("*.png")))
    qc_page_count = len(list((output_root / "qc" / "pages").glob("*.png")))
    candidate_aggregate_path = output_root / "tables" / "ecm_parameter_sweep_aggregate.csv"
    lines = [
        "# CellProfiler ECM CH1 Run Summary",
        "",
        f"Processed ECM image rows: {len(measurements)}",
        f"Selected CellProfiler Robust Background deviations: {selected_deviations:g}",
        f"Input root: `{input_root}`",
        "",
        "## Outputs",
        "",
        "- `final/tables/cellprofiler_ecm_ch1_summary.csv`",
        "- `final/workbooks/cellprofiler_ecm_ch1_summary.xlsx`",
        "- `tables/ecm_parameter_sweep_aggregate.csv`",
        "- `qc/parameter_selection/ecm_parameter_sweep_panel.png`",
        "- `qc/parameter_selection/ecm_parameter_sweep_summary.png`",
        "- `qc/review/ecm_representative_review_panel.png`",
        "- `qc/pages/`",
        "",
        "## Endpoint",
        "",
        f"`{ECM_ENDPOINT}` = sum of CH1 ECM signal inside the CellProfiler ECM-positive "
        "mask after subtracting the per-image mask-negative median CH1 background. "
        "No DAPI or nuclei normalization is applied for ECM.",
        "",
        "## Automated Integrity Checks",
        "",
        f"- ECM records discovered: {len(records)}",
        f"- final CellProfiler mask files: {mask_count}",
        f"- final CellProfiler overlay files: {overlay_count}",
        f"- acquisition QC pages: {qc_page_count}",
        f"- parameter aggregate table exists: {candidate_aggregate_path.exists()}",
        "",
        "## Per-Acquisition Endpoint Summary",
        "",
        *per_acquisition_summary_lines(by_acquisition),
        "",
        "## QC Status Counts",
        "",
        *[f"- {key}: {value}" for key, value in sorted(status_counts.items())],
        "",
        "## QC Flag Counts",
        "",
        *([f"- {key}: {value}" for key, value in sorted(flag_counts.items())] or ["- none"]),
        "",
        "## Review Fields",
        "",
        *(
            [
                f"- {measurement.acquisition} {measurement.location}: {measurement.qc_flags}"
                for measurement in review_rows
            ]
            or ["- none"]
        ),
        "",
        "## Interpretation",
        "",
        "CellProfiler generated the ECM-positive masks from staged CH1 images. The project "
        "postprocessing calculated the background-corrected ECM endpoint using the median "
        "of mask-negative CH1 pixels as the per-image background value. Manual ground-truth "
        "ECM masks were not available, so this run does not report precision, recall, F1, "
        "false positives, false negatives, or IoU.",
    ]
    markdown = "\n".join(lines) + "\n"
    (final_dir / "START_HERE_ECM_RUN_SUMMARY.md").write_text(markdown, encoding="utf-8")
    html = "<html><body>" + markdown_to_simple_html(markdown) + "</body></html>\n"
    (final_dir / "START_HERE_ECM_RUN_SUMMARY.html").write_text(html, encoding="utf-8")
    (output_root / "logs" / "run_record.yaml").write_text(
        yaml.safe_dump(
            {
                "input_root": str(input_root),
                "output_root": str(output_root),
                "records_processed": len(records),
                "selected_deviations": selected_deviations,
                "candidate_deviations": sorted({float(row["candidate_deviations"]) for row in candidate_summaries}),
                "qc_status_counts": dict(status_counts),
                "qc_flag_counts": dict(flag_counts),
                "final_mask_files": mask_count,
                "final_overlay_files": overlay_count,
                "qc_page_files": qc_page_count,
                "primary_endpoint": ECM_ENDPOINT,
                "normalization_denominator": "none",
                "dapi_count_role": "context_only_not_used_for_ecm_normalization",
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )


def per_acquisition_summary_lines(by_acquisition: dict[str, list[EcmMeasurement]]) -> list[str]:
    lines: list[str] = []
    for acquisition, measurements in sorted(by_acquisition.items()):
        endpoints = np.array(
            [
                measurement.ecm_positive_integrated_background_corrected
                for measurement in measurements
            ],
            dtype=float,
        )
        if endpoints.size:
            lines.append(
                f"- {acquisition}: rows={len(measurements)}, endpoint rows={endpoints.size}, "
                f"median={np.median(endpoints):.3g}, q1={np.quantile(endpoints, 0.25):.3g}, "
                f"q3={np.quantile(endpoints, 0.75):.3g}"
            )
        else:
            lines.append(f"- {acquisition}: rows={len(measurements)}, endpoint rows=0")
    return lines


def markdown_to_simple_html(markdown: str) -> str:
    html_lines: list[str] = []
    in_list = False
    for line in markdown.splitlines():
        if line.startswith("# "):
            if in_list:
                html_lines.append("</ul>")
                in_list = False
            html_lines.append(f"<h1>{line[2:]}</h1>")
        elif line.startswith("## "):
            if in_list:
                html_lines.append("</ul>")
                in_list = False
            html_lines.append(f"<h2>{line[3:]}</h2>")
        elif line.startswith("- "):
            if not in_list:
                html_lines.append("<ul>")
                in_list = True
            html_lines.append(f"<li>{line[2:]}</li>")
        elif line.strip():
            if in_list:
                html_lines.append("</ul>")
                in_list = False
            html_lines.append(f"<p>{line}</p>")
        else:
            if in_list:
                html_lines.append("</ul>")
                in_list = False
    if in_list:
        html_lines.append("</ul>")
    return "\n".join(html_lines)


def measurement_lookup(measurements: list[EcmMeasurement]) -> dict[tuple[str, str], EcmMeasurement]:
    return {(measurement.acquisition, measurement.location): measurement for measurement in measurements}


def dapi_count_for(record: EcmImageRecord, dapi_counts: dict[tuple[str, str], int]) -> int | None:
    return dapi_counts.get((record.acquisition_name, record.location)) or dapi_counts.get(
        (record.acquisition_name.split("/")[-1], record.location)
    )


def describe_cellprofiler(
    *,
    cellprofiler_app: Path,
    cellprofiler_executable: Path | None,
) -> str:
    cp_executable = resolve_cellprofiler_executable(
        cellprofiler_app=cellprofiler_app,
        cellprofiler_executable=cellprofiler_executable,
    )
    java_home = cellprofiler_app.expanduser().resolve() / "Contents" / "Resources" / "Home"
    lines = [
        f"cellprofiler_app={cellprofiler_app}",
        f"cellprofiler_executable={cp_executable}",
        f"cellprofiler_executable_source={'explicit' if cellprofiler_executable else 'app_bundle'}",
    ]
    if cellprofiler_executable is None:
        lines.append(f"bundled_java_home={java_home}")
    if cp_executable.exists():
        result = subprocess.run(["file", str(cp_executable)], check=False, capture_output=True, text=True)
        lines.append(result.stdout.strip())
        if cellprofiler_executable is not None:
            runtime_result = subprocess.run(
                [str(cp_executable), "--runtime-info"],
                check=False,
                capture_output=True,
                text=True,
            )
            if runtime_result.returncode == 0 and runtime_result.stdout.strip():
                lines.append(runtime_result.stdout.strip())
            else:
                version_result = subprocess.run(
                    [str(cp_executable), "--version"],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                version_text = "\n".join(
                    part.strip() for part in [version_result.stdout, version_result.stderr] if part.strip()
                )
                if version_text:
                    lines.append(version_text)
        elif cellprofiler_executable is None:
            version_result = subprocess.run(
                [str(cp_executable), "--version"],
                check=False,
                capture_output=True,
                text=True,
            )
            version_text = "\n".join(
                part.strip() for part in [version_result.stdout, version_result.stderr] if part.strip()
            )
            if version_text:
                lines.append(version_text)
    java = java_home / "bin" / "java"
    if cellprofiler_executable is None and java.exists():
        result = subprocess.run([str(java), "-version"], check=False, capture_output=True, text=True)
        lines.append(result.stderr.strip())
    return "\n".join(lines) + "\n"


def format_metric(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.2e}"


def k_label(value: float) -> str:
    return f"k={value:g}"


if __name__ == "__main__":
    main()
