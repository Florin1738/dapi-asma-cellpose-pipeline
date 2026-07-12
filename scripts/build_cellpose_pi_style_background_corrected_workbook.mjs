#!/usr/bin/env node
import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const args = parseArgs(process.argv.slice(2));
if (!args.summary || !args.output) {
  console.error(
    "usage: build_cellpose_pi_style_background_corrected_workbook.mjs --summary <csv> --output <xlsx>",
  );
  process.exit(2);
}

const summaryPath = path.resolve(args.summary);
const outputPath = path.resolve(args.output);
const records = parseCsv(await fs.readFile(summaryPath, "utf8"));
if (records.length === 0) {
  throw new Error(`No rows found in ${summaryPath}`);
}

const workbook = Workbook.create();
workbook.comments.setSelf({ displayName: "User" });

const readme = workbook.worksheets.add("README");
const plates = unique(records.map((row) => row.plate)).sort(plateSort);
const sheets = [readme];

writeReadme(readme, records, summaryPath);
for (const plate of plates) {
  const sheet = workbook.worksheets.add(plate);
  sheets.push(sheet);
  const rows = records.filter((row) => row.plate === plate).sort(fieldSort);
  writePlateSheet(sheet, rows, plate);
}

for (const sheet of sheets) {
  sheet.showGridLines = false;
}

await fs.mkdir(path.dirname(outputPath), { recursive: true });
const exported = await SpreadsheetFile.exportXlsx(workbook);
await exported.save(outputPath);

const preview = await workbook.render({
  sheetName: plates[0] || "Plate 1",
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
  sheet.getRange("A1:F1").merge();
  sheet.getRange("A1").values = [["Cellpose Retained-Region aSMA/DAPI Summary"]];
  sheet.getRange("A1").format = {
    fill: "#17324D",
    font: { bold: true, color: "#FFFFFF", size: 16 },
  };
  sheet.getRange("A3:B12").values = [
    ["Source CSV", sourceCsv],
    ["Generated", new Date().toISOString()],
    ["Fields", rows.length],
    ["Plates", plateCounts],
    ["CH2", "Alpha smooth muscle actin (aSMA) signal channel"],
    ["CH4", "DAPI nuclei channel"],
    [
      "Cellpose definition",
      "Cellpose means retained Cellpose objects with at least one DAPI-positive nucleus centroid inside the object.",
    ],
    ["Intensity column", "Cellpose retained-region CH2 integrated background-corrected intensity"],
    ["Ratio column", "Cellpose retained-region intensity divided by DAPI-positive nucleus count"],
    [
      "Background correction note",
      "For the current full-plate run, background_value_per_px = 0.0; Cellpose correction means pixels outside the retained-region mask were excluded.",
    ],
  ];
  sheet.getRange("A3:A12").format = {
    fill: "#E8EEF5",
    font: { bold: true },
  };
  sheet.getRange("A3:B12").format.borders = {
    preset: "inside",
    style: "thin",
    color: "#D4DEE8",
  };
  sheet.getRange("A14:F14").merge();
  sheet.getRange("A14").values = [[
    "First four columns match the PI-style workbook structure: LOCATION, aSMA intensity, Nuclei Count, Ratio. The final XY Location column is included only to make well lookup easier.",
  ]];
  sheet.getRange("A14").format = {
    fill: "#FFF7D6",
    font: { bold: true, color: "#5F4100" },
    wrapText: true,
  };
  sheet.getRange("A:A").format.columnWidth = 28;
  sheet.getRange("B:B").format.columnWidth = 110;
}

function writePlateSheet(sheet, rows, plate) {
  const headers = [
    "LOCATION",
    "aSMA intensity",
    "Nuclei Count",
    "Ratio",
    "XY Location",
  ];
  const values = rows.map((row) => [
    row.source_id,
    number(row.dapi_anchored_cellpose_ch2_integrated_background_corrected),
    number(row.dapi_positive_nucleus_count),
    null,
    row.location,
  ]);
  const matrix = [headers, ...values];
  const lastRow = matrix.length;
  sheet.getRange(`A1:E${lastRow}`).values = matrix;
  if (rows.length > 0) {
    sheet.getRange(`D2:D${lastRow}`).formulas = rows.map((_row, index) => [
      `=IFERROR(B${index + 2}/C${index + 2},"")`,
    ]);
  }
  const table = sheet.tables.add(`A1:E${lastRow}`, true, sanitizeTableName(`${plate}CellposeSummary`));
  table.style = "TableStyleMedium2";
  table.showFilterButton = true;
  sheet.freezePanes.freezeRows(1);
  sheet.getRange("A1:E1").format = {
    fill: "#17324D",
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
  };
  sheet.getRange("1:1").format.rowHeight = 42;
  sheet.getRange(`A1:E${lastRow}`).format.borders = {
    preset: "inside",
    style: "thin",
    color: "#D4DEE8",
  };
  sheet.getRange("A:A").format.columnWidth = 34;
  sheet.getRange("B:B").format.columnWidth = 18;
  sheet.getRange("C:D").format.columnWidth = 16;
  sheet.getRange("E:E").format.columnWidth = 14;
  sheet.getRange("B:B").format.numberFormat = "0.00E+00";
  sheet.getRange("C:C").format.numberFormat = "#,##0";
  sheet.getRange("D:D").format.numberFormat = "0.00E+00";
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
