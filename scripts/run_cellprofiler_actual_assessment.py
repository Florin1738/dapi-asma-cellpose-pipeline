from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess
import textwrap

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from skimage import measure, segmentation
import tifffile

from dapi_norm.image_arrays import read_primary_intensity_plane


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "output" / "cellprofiler_actual_assessment"
DEFAULT_SOURCE_ROOT = (
    REPO_ROOT
    / "data"
    / "aSMA_DAPI_plates"
    / "plate 1"
    / "ApYYM20AGGSMA_01"
)
DEFAULT_CELLPROFILER_APP = Path.home() / "Applications" / "CellProfiler.app"
DEFAULT_FIELDS = ["XY22", "XY23", "XY24", "XY40", "XY41"]
CELLPOSE_NUCLEI_COUNTS = (
    REPO_ROOT
    / "output"
    / "pi_simple_summary"
    / "cellpose_counts"
    / "Plate_1"
    / "ApYYM20AGGSMA_01"
    / "summaries"
    / "nucleus_counts.csv"
)
CELLPOSE_REGION_SUMMARY = (
    REPO_ROOT
    / "output"
    / "cellpose_cell_regions"
    / "full_plate_cpsam_v2"
    / "Plate_1"
    / "ApYYM20AGGSMA_01"
    / "summaries"
    / "cellpose_cell_region_image_metrics.csv"
)
CELLPOSE_REGION_MASK_DIR = (
    REPO_ROOT
    / "output"
    / "cellpose_cell_regions"
    / "full_plate_cpsam_v2"
    / "Plate_1"
    / "ApYYM20AGGSMA_01"
    / "masks"
)
CELLPOSE_NUCLEI_MASK_DIR = (
    REPO_ROOT
    / "output"
    / "pi_simple_summary"
    / "cellpose_counts"
    / "Plate_1"
    / "ApYYM20AGGSMA_01"
    / "masks"
)
METHOD_REVIEW_SUMMARY = (
    REPO_ROOT
    / "output"
    / "method_review"
    / "plate1_selected15_propagation_vs_cellpose"
    / "method_comparison_review_summary.csv"
)


def main() -> None:
    args = parse_args()
    fields = [field.upper() for field in args.fields]
    output_root = args.output_root.resolve()
    staging_dir = output_root / "staging_grayscale"
    pipeline_dir = output_root / "pipelines"
    cp_run_dir = output_root / "cp_run"
    logs_dir = output_root / "logs"
    masks_dir = output_root / "masks"
    qc_dir = output_root / "qc"
    tables_dir = output_root / "tables"
    exports_dir = output_root / "exports"

    for path in [
        staging_dir,
        pipeline_dir,
        cp_run_dir,
        logs_dir,
        masks_dir,
        qc_dir,
        tables_dir,
        exports_dir,
    ]:
        path.mkdir(parents=True, exist_ok=True)

    manifest_path = logs_dir / "resolved_input_manifest.csv"
    stage_grayscale_inputs(
        source_root=args.source_root.resolve(),
        staging_dir=staging_dir,
        fields=fields,
        manifest_path=manifest_path,
    )

    threshold_suffix = (
        "otsu" if args.primary_threshold_mode == "otsu" else f"manual{args.manual_primary_threshold:g}"
    )
    pipeline_path = pipeline_dir / f"plate1_dapi_asma_native_cellprofiler_{threshold_suffix}.cppipe"
    pipeline_path.write_text(
        build_pipeline(
            primary_threshold_mode=args.primary_threshold_mode,
            manual_primary_threshold=args.manual_primary_threshold,
        ),
        encoding="utf-8",
    )

    version_path = logs_dir / "cellprofiler_version.txt"
    version_path.write_text(describe_cellprofiler(args.cellprofiler_app), encoding="utf-8")

    if not args.skip_cellprofiler:
        reset_directory(cp_run_dir)
        cp_run_dir.mkdir(parents=True, exist_ok=True)
        run_cellprofiler(
            cellprofiler_app=args.cellprofiler_app,
            pipeline_path=pipeline_path,
            input_dir=staging_dir,
            output_dir=cp_run_dir,
            log_path=logs_dir / "run_log.txt",
        )

    copy_cellprofiler_exports(cp_run_dir=cp_run_dir, output_root=output_root)
    metrics = build_metrics(fields=fields, output_root=output_root, pipeline_path=pipeline_path)
    metrics_path = tables_dir / "cellprofiler_image_metrics.csv"
    write_csv(metrics_path, metrics)
    render_comparison_panel(
        fields=fields,
        metrics=metrics,
        output_root=output_root,
        output_path=qc_dir / "cellprofiler_actual_fullfield_comparison.png",
        crop=False,
    )
    render_comparison_panel(
        fields=fields,
        metrics=metrics,
        output_root=output_root,
        output_path=qc_dir / "cellprofiler_actual_crop_comparison.png",
        crop=True,
    )
    write_readme(output_root=output_root, fields=fields)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run an actual CellProfiler app pipeline on representative DAPI/aSMA fields."
    )
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--cellprofiler-app", type=Path, default=DEFAULT_CELLPROFILER_APP)
    parser.add_argument("--fields", nargs="+", default=DEFAULT_FIELDS)
    parser.add_argument("--skip-cellprofiler", action="store_true")
    parser.add_argument(
        "--primary-threshold-mode",
        choices=["otsu", "manual"],
        default="otsu",
        help="DAPI IdentifyPrimaryObjects threshold mode for actual CellProfiler nuclei.",
    )
    parser.add_argument(
        "--manual-primary-threshold",
        type=float,
        default=0.7,
        help="Manual DAPI threshold in CellProfiler-normalized 0-1 units when using manual mode.",
    )
    return parser.parse_args()


def stage_grayscale_inputs(
    *, source_root: Path, staging_dir: Path, fields: list[str], manifest_path: Path
) -> None:
    rows: list[dict[str, str]] = []
    for field in fields:
        field_dir = source_root / field
        ch2_path = field_dir / f"ApYYM20AGGSMA_{field}_CH2.tif"
        ch4_path = field_dir / f"ApYYM20AGGSMA_{field}_CH4.tif"
        if not ch2_path.exists():
            raise FileNotFoundError(ch2_path)
        if not ch4_path.exists():
            raise FileNotFoundError(ch4_path)
        asma, _ = read_primary_intensity_plane(ch2_path)
        dapi, _ = read_primary_intensity_plane(ch4_path)
        if asma.shape != dapi.shape:
            raise ValueError(f"Shape mismatch for {field}: CH2={asma.shape}, CH4={dapi.shape}")

        dapi_out = staging_dir / f"{field}_DAPI.tif"
        asma_out = staging_dir / f"{field}_ASMA.tif"
        tifffile.imwrite(dapi_out, dapi.astype(np.uint16), photometric="minisblack")
        tifffile.imwrite(asma_out, asma.astype(np.uint16), photometric="minisblack")
        rows.append(
            {
                "image_id": field,
                "source_ch2_asma_path": str(ch2_path.relative_to(REPO_ROOT)),
                "source_ch4_dapi_path": str(ch4_path.relative_to(REPO_ROOT)),
                "staged_asma_path": str(asma_out.relative_to(REPO_ROOT)),
                "staged_dapi_path": str(dapi_out.relative_to(REPO_ROOT)),
                "staging_note": "RGB pseudocolor TIFF collapsed to active single channel before CellProfiler",
            }
        )
    write_csv(manifest_path, rows)


def build_pipeline(*, primary_threshold_mode: str, manual_primary_threshold: float) -> str:
    pipeline = textwrap.dedent(
        """\
        CellProfiler Pipeline: http://www.cellprofiler.org
        Version:5
        DateRevision:428
        GitHash:
        ModuleCount:20
        HasImagePlaneDetails:False

        Images:[module_num:1|svn_version:'Unknown'|variable_revision_number:2|show_window:False|notes:['Collect staged grayscale DAPI and ASMA TIFFs for actual CellProfiler segmentation.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            :
            Filter images?:Images only
            Select the rule criteria:and (extension does isimage) (directory doesnot containregexp "[\\\\\\\\/]\\\\.")

        Metadata:[module_num:2|svn_version:'Unknown'|variable_revision_number:6|show_window:False|notes:['Extract image ID and staged channel from file names.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Extract metadata?:Yes
            Metadata data type:Text
            Metadata types:{}
            Extraction method count:1
            Metadata extraction method:Extract from file/folder names
            Metadata source:File name
            Regular expression to extract from file name:^(?P<Field>XY[0-9]{2})_(?P<Channel>DAPI|ASMA)\\.tif$
            Regular expression to extract from folder name:(?P<Date>[0-9]{4}_[0-9]{2}_[0-9]{2})$
            Extract metadata from:All images
            Select the filtering criteria:and (file does contain "")
            Metadata file location:Elsewhere...|
            Match file and image metadata:[]
            Use case insensitive matching?:No
            Metadata file name:
            Does cached metadata exist?:No

        NamesAndTypes:[module_num:3|svn_version:'Unknown'|variable_revision_number:8|show_window:False|notes:['Assign DAPI and ASMA grayscale channels.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Assign a name to:Images matching rules
            Select the image type:Grayscale image
            Name to assign these images:DAPI
            Match metadata:[]
            Image set matching method:Order
            Set intensity range from:Image metadata
            Assignments count:2
            Single images count:0
            Maximum intensity:65535.0
            Process as 3D?:No
            Relative pixel spacing in X:1.0
            Relative pixel spacing in Y:1.0
            Relative pixel spacing in Z:1.0
            Select the rule criteria:and (file does contain "_DAPI.tif")
            Name to assign these images:DAPI
            Name to assign these objects:Cell
            Select the image type:Grayscale image
            Set intensity range from:Image metadata
            Maximum intensity:65535.0
            Select the rule criteria:and (file does contain "_ASMA.tif")
            Name to assign these images:ASMA
            Name to assign these objects:Cell
            Select the image type:Grayscale image
            Set intensity range from:Image metadata
            Maximum intensity:65535.0

        Groups:[module_num:4|svn_version:'Unknown'|variable_revision_number:2|show_window:False|notes:['No grouping for the small representative run.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Do you want to group your images?:No
            grouping metadata count:1
            Metadata category:None

        IdentifyPrimaryObjects:[module_num:5|svn_version:'Unknown'|variable_revision_number:14|show_window:True|notes:['Identify DAPI nuclei using native CellProfiler.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the input image:DAPI
            Name the primary objects to be identified:Nuclei
            Typical diameter of objects, in pixel units (Min,Max):8,80
            Discard objects outside the diameter range?:Yes
            Discard objects touching the border of the image?:Yes
            Method to distinguish clumped objects:Intensity
            Method to draw dividing lines between clumped objects:Intensity
            Size of smoothing filter:10
            Suppress local maxima that are closer than this minimum allowed distance:7.0
            Speed up by using lower-resolution image to find local maxima?:Yes
            Fill holes in identified objects?:After declumping only
            Automatically calculate size of smoothing filter for declumping?:Yes
            Automatically calculate minimum allowed distance between local maxima?:Yes
            Handling of objects if excessive number of objects identified:Continue
            Maximum number of objects:1000
            Display accepted local maxima?:No
            Select maxima color:Blue
            Use advanced settings?:Yes
            Threshold setting version:11
            Threshold strategy:Global
            Thresholding method:Otsu
            Threshold smoothing scale:1.3488
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
            # of deviations:2.0
            Thresholding method:Otsu

        IdentifySecondaryObjects:[module_num:6|svn_version:'Unknown'|variable_revision_number:10|show_window:True|notes:['ASMA-guided CellProfiler secondary objects via Propagation. Exploratory because ASMA is also the measurement channel.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the input objects:Nuclei
            Name the objects to be identified:CellsPropagation
            Select the method to identify the secondary objects:Propagation
            Select the input image:ASMA
            Number of pixels by which to expand the primary objects:20
            Regularization factor:0.05
            Discard secondary objects touching the border of the image?:No
            Discard the associated primary objects?:No
            Name the new primary objects:FilteredNucleiPropagation
            Fill holes in identified objects?:Yes
            Threshold setting version:11
            Threshold strategy:Global
            Thresholding method:Minimum Cross-Entropy
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
            # of deviations:2.0
            Thresholding method:Otsu

        IdentifySecondaryObjects:[module_num:7|svn_version:'Unknown'|variable_revision_number:10|show_window:True|notes:['Fixed-distance CellProfiler secondary objects from DAPI nuclei via Distance - N.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the input objects:Nuclei
            Name the objects to be identified:CellsDistanceN
            Select the method to identify the secondary objects:Distance - N
            Select the input image:DAPI
            Number of pixels by which to expand the primary objects:20
            Regularization factor:0.05
            Discard secondary objects touching the border of the image?:No
            Discard the associated primary objects?:No
            Name the new primary objects:FilteredNucleiDistanceN
            Fill holes in identified objects?:Yes
            Threshold setting version:11
            Threshold strategy:Global
            Thresholding method:Otsu
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
            # of deviations:2.0
            Thresholding method:Otsu

        MeasureObjectIntensity:[module_num:8|svn_version:'Unknown'|variable_revision_number:4|show_window:True|notes:['Measure DAPI and ASMA intensity inside CellProfiler objects.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select images to measure:DAPI, ASMA
            Select objects to measure:Nuclei, CellsPropagation, CellsDistanceN

        MeasureObjectSizeShape:[module_num:9|svn_version:'Unknown'|variable_revision_number:3|show_window:True|notes:['Measure CellProfiler object size and shape.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select object sets to measure:Nuclei, CellsPropagation, CellsDistanceN
            Calculate the Zernike features?:No
            Calculate the advanced features?:No

        ConvertObjectsToImage:[module_num:10|svn_version:'Unknown'|variable_revision_number:1|show_window:True|notes:['Convert CellProfiler nuclei objects to uint16 label image.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the input objects:Nuclei
            Name the output image:NucleiImage
            Select the color format:uint16
            Select the colormap:Default

        SaveImages:[module_num:11|svn_version:'Unknown'|variable_revision_number:16|show_window:True|notes:['Save CellProfiler nuclei labels.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the type of image to save:Image
            Select the image to save:NucleiImage
            Select method for constructing file names:From image filename
            Select image name for file prefix:DAPI
            Enter single file name:OrigBlue
            Number of digits:4
            Append a suffix to the image file name?:Yes
            Text to append to the image name:_CellProfilerNucleiLabels
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

        ConvertObjectsToImage:[module_num:12|svn_version:'Unknown'|variable_revision_number:1|show_window:True|notes:['Convert Propagation cells to uint16 label image.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the input objects:CellsPropagation
            Name the output image:CellsPropagationImage
            Select the color format:uint16
            Select the colormap:Default

        SaveImages:[module_num:13|svn_version:'Unknown'|variable_revision_number:16|show_window:True|notes:['Save CellProfiler Propagation labels.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the type of image to save:Image
            Select the image to save:CellsPropagationImage
            Select method for constructing file names:From image filename
            Select image name for file prefix:ASMA
            Enter single file name:OrigBlue
            Number of digits:4
            Append a suffix to the image file name?:Yes
            Text to append to the image name:_CellProfilerPropagationLabels
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

        ConvertObjectsToImage:[module_num:14|svn_version:'Unknown'|variable_revision_number:1|show_window:True|notes:['Convert Distance - N cells to uint16 label image.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the input objects:CellsDistanceN
            Name the output image:CellsDistanceNImage
            Select the color format:uint16
            Select the colormap:Default

        SaveImages:[module_num:15|svn_version:'Unknown'|variable_revision_number:16|show_window:True|notes:['Save CellProfiler Distance - N labels.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the type of image to save:Image
            Select the image to save:CellsDistanceNImage
            Select method for constructing file names:From image filename
            Select image name for file prefix:DAPI
            Enter single file name:OrigBlue
            Number of digits:4
            Append a suffix to the image file name?:Yes
            Text to append to the image name:_CellProfilerDistanceNLabels
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

        OverlayOutlines:[module_num:16|svn_version:'Unknown'|variable_revision_number:4|show_window:True|notes:['Overlay Propagation cell and nucleus outlines on ASMA.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Display outlines on a blank image?:No
            Select image on which to display outlines:ASMA
            Name the output image:PropagationOverlay
            Outline display mode:Color
            Select method to determine brightness of outlines:Max of image
            How to outline:Thick
            Select outline color:Blue
            Select objects to display:CellsPropagation
            Select outline color:yellow
            Select objects to display:Nuclei

        SaveImages:[module_num:17|svn_version:'Unknown'|variable_revision_number:16|show_window:True|notes:['Save CellProfiler Propagation overlay.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the type of image to save:Image
            Select the image to save:PropagationOverlay
            Select method for constructing file names:From image filename
            Select image name for file prefix:ASMA
            Enter single file name:OrigBlue
            Number of digits:4
            Append a suffix to the image file name?:Yes
            Text to append to the image name:_CellProfilerPropagationOverlay
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

        OverlayOutlines:[module_num:18|svn_version:'Unknown'|variable_revision_number:4|show_window:True|notes:['Overlay Distance - N cell and nucleus outlines on ASMA.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Display outlines on a blank image?:No
            Select image on which to display outlines:ASMA
            Name the output image:DistanceNOverlay
            Outline display mode:Color
            Select method to determine brightness of outlines:Max of image
            How to outline:Thick
            Select outline color:Green
            Select objects to display:CellsDistanceN
            Select outline color:yellow
            Select objects to display:Nuclei

        SaveImages:[module_num:19|svn_version:'Unknown'|variable_revision_number:16|show_window:True|notes:['Save CellProfiler Distance - N overlay.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
            Select the type of image to save:Image
            Select the image to save:DistanceNOverlay
            Select method for constructing file names:From image filename
            Select image name for file prefix:ASMA
            Enter single file name:OrigBlue
            Number of digits:4
            Append a suffix to the image file name?:Yes
            Text to append to the image name:_CellProfilerDistanceNOverlay
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

        ExportToSpreadsheet:[module_num:20|svn_version:'Unknown'|variable_revision_number:13|show_window:True|notes:['Export actual CellProfiler measurements.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
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
    if primary_threshold_mode == "manual":
        pipeline = pipeline.replace(
            "    Thresholding method:Otsu\n"
            "    Threshold smoothing scale:1.3488\n",
            "    Thresholding method:Manual\n"
            "    Threshold smoothing scale:1.3488\n",
            1,
        )
        pipeline = pipeline.replace(
            "    Manual threshold:0.0\n",
            f"    Manual threshold:{manual_primary_threshold:g}\n",
            1,
        )
    return pipeline


def describe_cellprofiler(cellprofiler_app: Path) -> str:
    info_plist = cellprofiler_app / "Contents" / "Info.plist"
    cp_executable = cellprofiler_app / "Contents" / "MacOS" / "cp"
    java_home = cellprofiler_app / "Contents" / "Resources" / "Home"
    lines = [
        f"cellprofiler_app={cellprofiler_app}",
        f"cellprofiler_executable={cp_executable}",
        f"bundled_java_home={java_home}",
    ]
    if info_plist.exists():
        result = subprocess.run(
            ["plutil", "-extract", "CFBundleShortVersionString", "raw", str(info_plist)],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            lines.append(f"CFBundleShortVersionString={result.stdout.strip()}")
    if cp_executable.exists():
        result = subprocess.run(["file", str(cp_executable)], check=False, capture_output=True, text=True)
        lines.append(result.stdout.strip())
    java = java_home / "bin" / "java"
    if java.exists():
        result = subprocess.run([str(java), "-version"], check=False, capture_output=True, text=True)
        lines.append(result.stderr.strip())
    return "\n".join(lines) + "\n"


def run_cellprofiler(
    *, cellprofiler_app: Path, pipeline_path: Path, input_dir: Path, output_dir: Path, log_path: Path
) -> None:
    cp_executable = cellprofiler_app / "Contents" / "MacOS" / "cp"
    java_home = cellprofiler_app / "Contents" / "Resources" / "Home"
    if not cp_executable.exists():
        raise FileNotFoundError(cp_executable)
    if not java_home.exists():
        raise FileNotFoundError(java_home)
    env = os.environ.copy()
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


def copy_cellprofiler_exports(*, cp_run_dir: Path, output_root: Path) -> None:
    masks_dir = output_root / "masks"
    exports_dir = output_root / "exports"
    qc_dir = output_root / "qc"
    for path in cp_run_dir.iterdir():
        if not path.is_file():
            continue
        lower_name = path.name.lower()
        if lower_name.endswith(".csv"):
            shutil.copy2(path, exports_dir / path.name)
        elif "labels" in lower_name and lower_name.endswith((".tif", ".tiff")):
            shutil.copy2(path, masks_dir / normalize_cp_output_name(path.name))
        elif "overlay" in lower_name and lower_name.endswith(".png"):
            shutil.copy2(path, qc_dir / normalize_cp_output_name(path.name))


def normalize_cp_output_name(name: str) -> str:
    return name.replace("_DAPI_", "_").replace("_ASMA_", "_")


def build_metrics(*, fields: list[str], output_root: Path, pipeline_path: Path) -> list[dict[str, object]]:
    cellpose_counts = read_lookup(CELLPOSE_NUCLEI_COUNTS, key="image_id")
    cellpose_summary = read_lookup(CELLPOSE_REGION_SUMMARY, key="image_id")
    rows: list[dict[str, object]] = []
    for field in fields:
        asma = read_image(output_root / "staging_grayscale" / f"{field}_ASMA.tif")
        nuclei = read_label(find_one(output_root / "masks", f"{field}_CellProfilerNucleiLabels*.tif*"))
        propagation = read_label(
            find_one(output_root / "masks", f"{field}_CellProfilerPropagationLabels*.tif*")
        )
        distance_n = read_label(
            find_one(output_root / "masks", f"{field}_CellProfilerDistanceNLabels*.tif*")
        )
        cellpose = read_label(CELLPOSE_REGION_MASK_DIR / f"{field}_cellpose_ch2_ch4_cpsam_v2_labels.tif")
        cellpose_nuclei = read_label(CELLPOSE_NUCLEI_MASK_DIR / f"{field}_CH4_cpsam_v2_labels.tif")
        anchored_cellpose = anchor_cellpose_by_dapi_centroid(
            cellpose_labels=cellpose, nuclei_labels=cellpose_nuclei
        )
        denominator = int(cellpose_counts[field]["nucleus_count"])
        base = {
            "image_id": field,
            "source_id": field,
            "dapi_positive_nucleus_count": denominator,
            "normalization_denominator_count": denominator,
            "cellprofiler_pipeline_path": str(
                pipeline_path.relative_to(REPO_ROOT)
            ),
            "cellprofiler_version": "4.2.8 app bundle",
            "cellprofiler_primary_object_count": count_labels(nuclei),
            "cellpose_object_count": int(cellpose_summary[field]["dapi_anchored_cellpose_object_count"]),
            "dapi_anchored_cellpose_masked_area_px": int(
                cellpose_summary[field]["dapi_anchored_candidate_region_area_px"]
            ),
            "dapi_anchored_cellpose_ch2_integrated_raw": float(
                cellpose_summary[field]["dapi_anchored_target_integrated_raw"]
            ),
            "dapi_anchored_cellpose_ch2_integrated_raw_per_DAPI_positive_nucleus": float(
                cellpose_summary[field][
                    "dapi_anchored_target_integrated_intensity_per_DAPI_positive_nucleus"
                ]
            ),
            "cellpose_mask_path": str(
                (CELLPOSE_REGION_MASK_DIR / f"{field}_cellpose_ch2_ch4_cpsam_v2_labels.tif").relative_to(REPO_ROOT)
            ),
        }
        rows.append(
            method_row(
                base=base,
                method="cellprofiler_propagation_asma_guided",
                labels=propagation,
                asma=asma,
                denominator=denominator,
                cellpose_mask=anchored_cellpose,
                mask_path=output_root / "masks" / find_one(output_root / "masks", f"{field}_CellProfilerPropagationLabels*.tif*").name,
            )
        )
        rows.append(
            method_row(
                base=base,
                method="cellprofiler_distance_n_20px",
                labels=distance_n,
                asma=asma,
                denominator=denominator,
                cellpose_mask=anchored_cellpose,
                mask_path=output_root / "masks" / find_one(output_root / "masks", f"{field}_CellProfilerDistanceNLabels*.tif*").name,
            )
        )
    return rows


def method_row(
    *,
    base: dict[str, object],
    method: str,
    labels: np.ndarray,
    asma: np.ndarray,
    denominator: int,
    cellpose_mask: np.ndarray,
    mask_path: Path,
) -> dict[str, object]:
    region = labels > 0
    area = int(np.count_nonzero(region))
    integrated_raw = float(np.sum(asma[region]))
    union = np.logical_or(region, cellpose_mask)
    both = np.logical_and(region, cellpose_mask)
    cp_only = np.logical_and(region, ~cellpose_mask)
    cellpose_only = np.logical_and(~region, cellpose_mask)
    return {
        **base,
        "cellprofiler_method": method,
        "cellprofiler_secondary_object_count": count_labels(labels),
        "cellprofiler_secondary_region_area_px": area,
        "cellprofiler_secondary_region_fraction": area / labels.size,
        "cellprofiler_ch2_integrated_raw": integrated_raw,
        "cellprofiler_target_integrated_intensity_per_DAPI_positive_nucleus": (
            integrated_raw / denominator if denominator else float("nan")
        ),
        "cellprofiler_mask_path": str(mask_path.relative_to(REPO_ROOT)),
        "both_region_area_px": int(np.count_nonzero(both)),
        "cellprofiler_only_area_px": int(np.count_nonzero(cp_only)),
        "cellpose_only_area_px": int(np.count_nonzero(cellpose_only)),
        "union_area_px": int(np.count_nonzero(union)),
        "cellprofiler_vs_cellpose_region_jaccard": (
            float(np.count_nonzero(both) / np.count_nonzero(union))
            if np.count_nonzero(union)
            else float("nan")
        ),
        "cellprofiler_to_cellpose_endpoint_ratio": (
            integrated_raw / float(base["dapi_anchored_cellpose_ch2_integrated_raw"])
            if float(base["dapi_anchored_cellpose_ch2_integrated_raw"]) > 0
            else float("nan")
        ),
        "comparison_status": "exploratory_not_manual_ground_truth",
    }


def render_comparison_panel(
    *,
    fields: list[str],
    metrics: list[dict[str, object]],
    output_root: Path,
    output_path: Path,
    crop: bool,
) -> None:
    metrics_by_field = {
        row["image_id"]: row for row in metrics if row["cellprofiler_method"] == "cellprofiler_propagation_asma_guided"
    }
    cols = [
        "DAPI + CP nuclei",
        "ASMA raw",
        "CP Propagation",
        "CP Distance-N",
        "Cellpose anchored",
        "Propagation vs Cellpose",
    ]
    fig, axes = plt.subplots(len(fields), len(cols), figsize=(22, 3.8 * len(fields)), squeeze=False)
    for row_idx, field in enumerate(fields):
        dapi = read_image(output_root / "staging_grayscale" / f"{field}_DAPI.tif")
        asma = read_image(output_root / "staging_grayscale" / f"{field}_ASMA.tif")
        nuclei = read_label(find_one(output_root / "masks", f"{field}_CellProfilerNucleiLabels*.tif*"))
        propagation = read_label(
            find_one(output_root / "masks", f"{field}_CellProfilerPropagationLabels*.tif*")
        )
        distance_n = read_label(
            find_one(output_root / "masks", f"{field}_CellProfilerDistanceNLabels*.tif*")
        )
        cellpose = read_label(CELLPOSE_REGION_MASK_DIR / f"{field}_cellpose_ch2_ch4_cpsam_v2_labels.tif")
        cellpose_nuclei = read_label(CELLPOSE_NUCLEI_MASK_DIR / f"{field}_CH4_cpsam_v2_labels.tif")
        anchored_cellpose = anchor_cellpose_by_dapi_centroid(
            cellpose_labels=cellpose, nuclei_labels=cellpose_nuclei
        )
        bounds = crop_bounds_for_field(field) if crop else None
        images = [dapi, asma, asma, asma, asma, asma]
        masks = [nuclei, None, propagation, distance_n, anchored_cellpose, None]
        for col_idx, title in enumerate(cols):
            ax = axes[row_idx, col_idx]
            base = crop_array(images[col_idx], bounds)
            ax.imshow(normalize_for_display(base), cmap="gray")
            ax.set_xticks([])
            ax.set_yticks([])
            if row_idx == 0:
                ax.set_title(title, fontsize=11)
            if col_idx == 0:
                ax.set_ylabel(field, rotation=0, labelpad=28, fontsize=10, weight="bold")
            if masks[col_idx] is not None:
                overlay_boundaries(ax, crop_array(masks[col_idx], bounds), color="cyan")
            if col_idx == 5:
                disagreement = disagreement_rgb(
                    crop_array(propagation > 0, bounds),
                    crop_array(anchored_cellpose, bounds),
                )
                ax.imshow(disagreement, alpha=0.65)
            if col_idx == 2:
                add_metrics_text(ax, metrics_by_field[field])
    fig.suptitle(
        "Actual CellProfiler native segmentation vs Cellpose anchored masks"
        + (" (matched crops)" if crop else " (full field)"),
        fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def add_metrics_text(ax: plt.Axes, row: dict[str, object]) -> None:
    text = (
        f"CP obj {row['cellprofiler_secondary_object_count']} | "
        f"area {int(row['cellprofiler_secondary_region_area_px'])}\n"
        f"per DAPI {float(row['cellprofiler_target_integrated_intensity_per_DAPI_positive_nucleus']):.2e}\n"
        f"Jaccard {float(row['cellprofiler_vs_cellpose_region_jaccard']):.2f}"
    )
    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.55, "pad": 3, "edgecolor": "none"},
    )


def crop_bounds_for_field(field: str) -> tuple[slice, slice]:
    rows = read_lookup(METHOD_REVIEW_SUMMARY, key="image_id")
    crop_box = rows[field]["crop_box"]
    y0, x0, y1, x1 = [int(part) for part in crop_box.split(",")]
    return (slice(y0, y1), slice(x0, x1))


def crop_array(arr: np.ndarray, bounds: tuple[slice, slice] | None) -> np.ndarray:
    if bounds is None:
        return arr
    return arr[bounds]


def overlay_boundaries(ax: plt.Axes, labels: np.ndarray, *, color: str) -> None:
    boundaries = segmentation.find_boundaries(labels > 0, mode="outer")
    rgba = np.zeros((*boundaries.shape, 4), dtype=float)
    if color == "cyan":
        rgba[..., 1] = 1.0
        rgba[..., 2] = 1.0
    else:
        rgba[..., 0] = 1.0
    rgba[..., 3] = boundaries.astype(float) * 0.9
    ax.imshow(rgba)


def disagreement_rgb(cp_mask: np.ndarray, cellpose_mask: np.ndarray) -> np.ndarray:
    both = np.logical_and(cp_mask, cellpose_mask)
    cp_only = np.logical_and(cp_mask, ~cellpose_mask)
    cellpose_only = np.logical_and(~cp_mask, cellpose_mask)
    rgb = np.zeros((*cp_mask.shape, 3), dtype=float)
    rgb[both] = [1.0, 1.0, 1.0]
    rgb[cp_only] = [0.1, 0.85, 1.0]
    rgb[cellpose_only] = [1.0, 0.15, 0.15]
    return rgb


def normalize_for_display(image: np.ndarray) -> np.ndarray:
    arr = image.astype(float)
    lo, hi = np.percentile(arr, [1, 99.5])
    if hi <= lo:
        return np.zeros_like(arr)
    return np.clip((arr - lo) / (hi - lo), 0, 1)


def anchor_cellpose_by_dapi_centroid(*, cellpose_labels: np.ndarray, nuclei_labels: np.ndarray) -> np.ndarray:
    keep: set[int] = set()
    for prop in measure.regionprops(nuclei_labels.astype(np.int32)):
        y, x = prop.centroid
        label = int(cellpose_labels[int(round(y)), int(round(x))])
        if label > 0:
            keep.add(label)
    if not keep:
        return np.zeros(cellpose_labels.shape, dtype=bool)
    return np.isin(cellpose_labels, list(keep))


def read_lookup(path: Path, *, key: str) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle)}


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_image(path: Path) -> np.ndarray:
    return np.asarray(tifffile.imread(path))


def read_label(path: Path) -> np.ndarray:
    return np.asarray(tifffile.imread(path)).astype(np.int32)


def find_one(root: Path, pattern: str) -> Path:
    matches = sorted(root.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected exactly one match for {pattern} in {root}, found {len(matches)}")
    return matches[0]


def count_labels(labels: np.ndarray) -> int:
    values = np.unique(labels)
    return int(np.count_nonzero(values))


def reset_directory(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)


def write_readme(*, output_root: Path, fields: list[str]) -> None:
    readme = output_root / "README.md"
    readme.write_text(
        textwrap.dedent(
            f"""\
            # Actual CellProfiler Assessment

            This folder was generated by `scripts/run_cellprofiler_actual_assessment.py`.

            Fields: {", ".join(fields)}

            Actual CellProfiler artifacts:

            - `pipelines/plate1_dapi_asma_native_cellprofiler_*.cppipe`
            - `logs/run_log.txt`
            - `logs/cellprofiler_version.txt`
            - `logs/resolved_input_manifest.csv`
            - `exports/*.csv`
            - `masks/*CellProfiler*Labels*.tif`
            - `qc/*CellProfiler*Overlay*.png`
            - `qc/cellprofiler_actual_fullfield_comparison.png`
            - `qc/cellprofiler_actual_crop_comparison.png`
            - `tables/cellprofiler_image_metrics.csv`

            Interpretation status: exploratory. The Propagation workflow uses ASMA as both
            boundary cue and measurement channel, so it is not a validated unbiased whole-cell
            segmentation. The Distance-N workflow is a fixed expansion from DAPI nuclei and is
            not true whole-cell segmentation either.
            """
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
