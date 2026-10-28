# Image, well and treatment identity

## Required microscope provenance

The user-facing single and parallel Cellpose runners validate well identity
before model loading or output creation, including in dry-run mode. Each
raw TIFF field folder must contain exactly one microscope marker named like
`_A01` (A-H, 01-12). Markers can be empty: their names carry the well identity.
Preserve them when downloading or extracting data. Missing, conflicting or
duplicate-alias markers stop the run. XY order is never a well map, and root
shortcuts are not executed or used as a fallback.

The canonical implementation is `src/dapi_norm/well_mapping.py`. Join each
measurement through its actual input image path to the field's marker.
Keep repeated fields from the same well as distinct rows. `field_in_well`
is the numeric-XY ordering within that well; `fields_in_well` counts all raw
TIFF fields bearing that well marker, including fields excluded from analysis.
These fields are provenance, not biological replicates or well averages.

Unrepresented wells and unpaired fields must be reported rather than fabricated.
A valid run need not contain exactly 96 fields or represent every well.

## Experimenter-supplied treatment map

The microscope tells us the well, not the compound. Put `plate_map.csv` inside
each acquisition folder, or pass `--plate-map <CSV>` to either runner. A shared
file is allowed only when acquisition leaf names are unique. Required columns:

```csv
assay_folder,well_id,compound,batch,dose,condition,tgf_beta,serum,cells_seeded
Example assay,A01,vehicle,n/a,0 uM,control,no,1%,5000
```

Use the exact original acquisition folder name and uppercase zero-padded well
ID. Join by assay plus well, never row number or scan order. Provide one row
per well, covering every represented raw well. Duplicate wells, invalid IDs,
missing coverage and blank metadata in supplied maps are errors. Use explicit
`n/a` where appropriate. The pipeline records the filename and file SHA-256.
Keep real plate maps outside Git; the synthetic template above is not a map.

Without a map, image quantification is allowed, treatment fields remain blank,
and `treatment_mapping_status` is `not_provided`. Never infer treatment names
from XY order. Treatment comparisons require a verified experimenter map and
separate control-image review; pipeline completion is not approval of biology.
When a map is supplied and validated, status is `provided`, which means schema
and join validation, not independent verification of what was pipetted.

## Output schema

The polished workbook preserves intensity in B, nucleus count in C and ratio
formula in D. Column A is now explicitly `Image field ID` rather than LOCATION.
Columns E-I are Well, Field in well, Fields in well, Treatment mapping and QC
status. Channel Mapping describes the stain channels. The separate Well Mapping
sheet has one row per measured image, including source ID and these CSV fields:

```text
assay_folder, field_id, well_id, field_in_well, fields_in_well,
well_marker, well_marker_sha256, treatment_mapping_status,
plate_map_source, plate_map_sha256, compound, batch, dose, condition,
tgf_beta, serum, cells_seeded
```

Both endpoint CSVs carry these fields. `final/tables/image_to_well.csv` records
all raw TIFF fields, not just successfully paired or capped analysis fields;
`mapping_validation.json` records represented wells, repeat counts, absent
wells and treatment status. Therefore mapping table and endpoint row counts
can legitimately differ. Hashes bind file bytes, while marker filenames bind
the well identity; an empty marker's hash alone is not a unique identifier.

This is a provenance change. It does not change segmentation, scientific
thresholds, normalization, or the meaning of QC flags. Generic historical
PI-summary tools retain their legacy export schema; use the user-facing
Cellpose runners for the mandatory mapped assay workflow.
