#!/usr/bin/env node
import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const args = parseArgs(process.argv.slice(2));
if (!args.summary || !args.output) {
  console.error("usage: build_cellpose_full_plate_workbook.mjs --summary <csv> --output <xlsx>");
  process.exit(2);
}

const summaryPath = path.resolve(args.summary);
const outputPath = path.resolve(args.output);
const csvText = await fs.readFile(summaryPath, "utf8");
const records = parseCsv(csvText);
if (records.length === 0) {
  throw new Error(`No rows found in ${summaryPath}`);
}

const workbook = Workbook.create();
workbook.comments.setSelf({ displayName: "User" });
const readme = workbook.worksheets.add("README");
const allFields = workbook.worksheets.add("All Fields");
const plates = unique(records.map((row) => row.plate)).sort(plateSort);
const plateSheets = new Map(plates.map((plate) => [plate, workbook.worksheets.add(plate)]));
const createdSheets = [readme, allFields, ...plateSheets.values()];

writeReadme(readme, records, summaryPath);
writeDataSheet(allFields, records, { includePlate: true, tableName: "AllFieldsTable" });
for (const plate of plates) {
  const rows = records.filter((row) => row.plate === plate).sort(fieldSort);
  writeDataSheet(plateSheets.get(plate), rows, {
    includePlate: false,
    tableName: sanitizeTableName(`${plate}Table`),
  });
}

for (const sheet of createdSheets) {
  sheet.showGridLines = false;
}

await fs.mkdir(path.dirname(outputPath), { recursive: true });
const exported = await SpreadsheetFile.exportXlsx(workbook);
await exported.save(outputPath);

const preview = await workbook.render({
  sheetName: plates[0] || "All Fields",
  autoCrop: "all",
  scale: 1,
  format: "png",
});
await fs.writeFile(
  outputPath.replace(/\.xlsx$/i, "_preview.png"),
  new Uint8Array(await preview.arrayBuffer()),
);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "formula error scan",
});
console.log(errors.ndjson);
console.log(`saved=${outputPath}`);

function writeReadme(sheet, rows, sourceCsv) {
  const plateCounts = unique(rows.map((row) => row.plate))
    .sort(plateSort)
    .map((plate) => `${plate}: ${rows.filter((row) => row.plate === plate).length} fields`)
    .join("; ");
  sheet.getRange("A1:H1").merge();
  sheet.getRange("A1").values = [["Cellpose CH2+CH4 aSMA/DAPI Full-Plate Workbook"]];
  sheet.getRange("A1").format = {
    fill: "#17324D",
    font: { bold: true, color: "#FFFFFF", size: 16 },
  };
  sheet.getRange("A3:B11").values = [
    ["Source CSV", sourceCsv],
    ["Generated", new Date().toISOString()],
    ["Fields", rows.length],
    ["Plates", plateCounts],
    ["CH2", "Alpha smooth muscle actin (aSMA) signal channel"],
    ["CH4", "DAPI nuclei channel"],
    [
      "Primary denominator",
      "Count of DAPI-positive nuclei, not DAPI fluorescence intensity",
    ],
    [
      "Cellpose aSMA intensity",
      "Raw CH2 integrated intensity summed inside Cellpose CH2+CH4 candidate regions",
    ],
    [
      "DAPI-anchored variant",
      "Cellpose objects are retained only if at least one DAPI nucleus centroid falls inside the object",
    ],
  ];
  sheet.getRange("A3:A11").format = {
    fill: "#E8EEF5",
    font: { bold: true },
  };
  sheet.getRange("A3:B11").format.borders = {
    preset: "inside",
    style: "thin",
    color: "#D4DEE8",
  };
  sheet.getRange("A13:H13").merge();
  sheet.getRange("A13").values = [[
    "Interpretation note: this workbook computes auditable image-analysis metrics. It does not establish biological significance, and the Cellpose CH2+CH4 regions should be treated as candidate aSMA-associated regions unless manually validated. Review QC_STATUS and QC_FLAGS before using any row.",
  ]];
  sheet.getRange("A13").format = {
    fill: "#FFF7D6",
    font: { bold: true, color: "#5F4100" },
    wrapText: true,
  };
  sheet.getRange("A:A").format.columnWidth = 28;
  sheet.getRange("B:B").format.columnWidth = 110;
}

function writeDataSheet(sheet, rows, { includePlate, tableName }) {
  const headers = [
    ...(includePlate ? ["PLATE"] : []),
    "SOURCE_ID",
    "LOCATION",
    "aSMA intensity - Cellpose raw CH2",
    "Nuclei Count",
    "Ratio - Cellpose intensity per DAPI-positive nucleus",
    "DAPI-anchored aSMA intensity",
    "DAPI-anchored Ratio",
    "Whole-field aSMA intensity",
    "Whole-field Ratio",
    "Cellpose positive area px",
    "Cellpose area per DAPI-positive nucleus",
    "DAPI-anchored positive area px",
    "DAPI-anchored area per DAPI-positive nucleus",
    "No-DAPI Cellpose objects excluded",
    "Cellpose object count",
    "QC_STATUS",
    "QC_FLAGS",
    "SOURCE_WARNINGS",
    "Cellpose mask path",
    "DAPI nuclei mask path",
    "Source QC panel path",
    "Source excluded-signal check path",
    "CH2 path",
    "CH4 path",
  ];
  const valueRows = rows.map((row) => [
    ...(includePlate ? [row.plate] : []),
    row.source_id,
    row.location,
    number(row.cellpose_masked_ch2_integrated_raw),
    number(row.dapi_positive_nucleus_count),
    null,
    number(row.dapi_anchored_cellpose_ch2_integrated_raw),
    null,
    number(row.whole_field_ch2_integrated_raw),
    null,
    number(row.cellpose_masked_area_px),
    null,
    number(row.dapi_anchored_cellpose_masked_area_px),
    null,
    number(row.no_dapi_cellpose_object_count_excluded_in_anchored_variant),
    number(row.cellpose_object_count),
    row.qc_status,
    row.qc_flags,
    row.source_warnings,
    row.cellpose_mask_path,
    row.dapi_nuclei_mask_path,
    row.source_qc_panel_path,
    row.source_excluded_signal_check_path,
    row.ch2_path,
    row.ch4_path,
  ]);
  const matrix = [headers, ...valueRows];
  const lastCol = colName(headers.length);
  const lastRow = matrix.length;
  sheet.getRange(`A1:${lastCol}${lastRow}`).values = matrix;

  const offset = includePlate ? 1 : 0;
  const firstDataRow = 2;
  const lastDataRow = Math.max(firstDataRow, lastRow);
  const cellposeIntensityCol = colName(3 + offset);
  const nucleiCol = colName(4 + offset);
  const cellposeRatioCol = colName(5 + offset);
  const anchoredIntensityCol = colName(6 + offset);
  const anchoredRatioCol = colName(7 + offset);
  const wholeIntensityCol = colName(8 + offset);
  const wholeRatioCol = colName(9 + offset);
  const areaCol = colName(10 + offset);
  const areaRatioCol = colName(11 + offset);
  const anchoredAreaCol = colName(12 + offset);
  const anchoredAreaRatioCol = colName(13 + offset);

  if (rows.length > 0) {
    sheet.getRange(`${cellposeRatioCol}${firstDataRow}:${cellposeRatioCol}${lastDataRow}`).formulas =
      rows.map((_row, index) => [
        `=IFERROR(${cellposeIntensityCol}${firstDataRow + index}/${nucleiCol}${firstDataRow + index},"")`,
      ]);
    sheet.getRange(`${anchoredRatioCol}${firstDataRow}:${anchoredRatioCol}${lastDataRow}`).formulas =
      rows.map((_row, index) => [
        `=IFERROR(${anchoredIntensityCol}${firstDataRow + index}/${nucleiCol}${firstDataRow + index},"")`,
      ]);
    sheet.getRange(`${wholeRatioCol}${firstDataRow}:${wholeRatioCol}${lastDataRow}`).formulas = rows.map(
      (_row, index) => [
        `=IFERROR(${wholeIntensityCol}${firstDataRow + index}/${nucleiCol}${firstDataRow + index},"")`,
      ],
    );
    sheet.getRange(`${areaRatioCol}${firstDataRow}:${areaRatioCol}${lastDataRow}`).formulas = rows.map(
      (_row, index) => [`=IFERROR(${areaCol}${firstDataRow + index}/${nucleiCol}${firstDataRow + index},"")`],
    );
    sheet.getRange(`${anchoredAreaRatioCol}${firstDataRow}:${anchoredAreaRatioCol}${lastDataRow}`).formulas =
      rows.map((_row, index) => [
        `=IFERROR(${anchoredAreaCol}${firstDataRow + index}/${nucleiCol}${firstDataRow + index},"")`,
      ]);
  }

  const table = sheet.tables.add(`A1:${lastCol}${lastRow}`, true, tableName);
  table.style = "TableStyleMedium2";
  table.showFilterButton = true;
  sheet.freezePanes.freezeRows(1);
  sheet.getRange(`A1:${lastCol}1`).format = {
    fill: "#17324D",
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
  };
  sheet.getRange("1:1").format.rowHeight = 46;
  sheet.getRange(`A1:${lastCol}${lastRow}`).format.borders = {
    preset: "inside",
    style: "thin",
    color: "#D4DEE8",
  };
  sheet.getRange(`${cellposeIntensityCol}:${wholeRatioCol}`).format.numberFormat = "0.00E+00";
  sheet.getRange(`${areaCol}:${anchoredAreaRatioCol}`).format.numberFormat = "#,##0.0";
  sheet.getRange(`${nucleiCol}:${nucleiCol}`).format.numberFormat = "#,##0";
  sheet.getRange(`A:${lastCol}`).format.autofitColumns();
  sheet.getRange(`A:${lastCol}`).format.wrapText = false;
  sheet.getRange(`A1:${lastCol}1`).format.wrapText = true;
  if (includePlate) {
    sheet.getRange("A:A").format.columnWidth = 14;
    sheet.getRange("B:B").format.columnWidth = 32;
    sheet.getRange("C:C").format.columnWidth = 14;
    sheet.getRange("D:N").format.columnWidth = 18;
  } else {
    sheet.getRange("A:A").format.columnWidth = 32;
    sheet.getRange("B:B").format.columnWidth = 14;
    sheet.getRange("C:M").format.columnWidth = 18;
  }
}

function parseArgs(argv) {
  const parsed = {};
  for (let index = 0; index < argv.length; index += 1) {
    const token = argv[index];
    if (token.startsWith("--")) {
      parsed[token.slice(2)] = argv[index + 1];
      index += 1;
    }
  }
  return parsed;
}

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let index = 0; index < text.length; index += 1) {
    const char = text[index];
    const next = text[index + 1];
    if (quoted) {
      if (char === '"' && next === '"') {
        field += '"';
        index += 1;
      } else if (char === '"') {
        quoted = false;
      } else {
        field += char;
      }
    } else if (char === '"') {
      quoted = true;
    } else if (char === ",") {
      row.push(field);
      field = "";
    } else if (char === "\n") {
      row.push(field);
      rows.push(row);
      row = [];
      field = "";
    } else if (char !== "\r") {
      field += char;
    }
  }
  if (field || row.length) {
    row.push(field);
    rows.push(row);
  }
  const headers = rows.shift() || [];
  return rows
    .filter((values) => values.some((value) => value !== ""))
    .map((values) => Object.fromEntries(headers.map((header, index) => [header, values[index] ?? ""])));
}

function unique(values) {
  return Array.from(new Set(values));
}

function number(value) {
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : null;
}

function fieldSort(left, right) {
  return xyNumber(left.location) - xyNumber(right.location) || left.source_id.localeCompare(right.source_id);
}

function plateSort(left, right) {
  return xyNumber(left) - xyNumber(right) || left.localeCompare(right);
}

function xyNumber(value) {
  const match = String(value).match(/\d+/);
  return match ? Number(match[0]) : 9999;
}

function sanitizeTableName(value) {
  return value.replace(/[^A-Za-z0-9_]/g, "");
}

function colName(index1) {
  let value = index1;
  let name = "";
  while (value > 0) {
    const rem = (value - 1) % 26;
    name = String.fromCharCode(65 + rem) + name;
    value = Math.floor((value - 1) / 26);
  }
  return name;
}
