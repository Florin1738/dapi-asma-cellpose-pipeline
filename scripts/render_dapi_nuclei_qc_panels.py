#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

from dapi_norm.cellpose_endpoint_figures import render_dapi_nuclei_qc_pages


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render comprehensive CH4/DAPI nuclei detection QC panels from the full-plate summary CSV."
    )
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--panel-page-size", type=int, default=12)
    args = parser.parse_args()

    outputs = render_dapi_nuclei_qc_pages(
        summary_csv=args.summary,
        output_dir=args.output,
        page_size=args.panel_page_size,
    )
    print(f"fields={outputs['field_count']}")
    print(f"plates={outputs['plate_count']}")
    print(f"pages={len(outputs['pages'])}")
    print(f"index={outputs['index']}")


if __name__ == "__main__":
    main()

