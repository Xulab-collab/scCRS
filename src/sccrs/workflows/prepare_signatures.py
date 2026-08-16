#!/usr/bin/env python3
"""Merge all sheets of Cytokine_Signatures_Human.xlsx into a scorer-ready CSV.

The official Immune Dictionary workbook has one cell type per worksheet.
Reading it with pandas.read_excel(path) selects only the first B_cell sheet.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REQUIRED = {"Celltype_Str", "Cytokine", "Gene_Human", "Avg_log2FC"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    workbook = pd.ExcelFile(args.input)
    tables = []
    for sheet in workbook.sheet_names:
        table = pd.read_excel(workbook, sheet_name=sheet)
        missing = REQUIRED.difference(table.columns)
        if missing:
            raise ValueError(f"Sheet {sheet!r} lacks {sorted(missing)}; columns are {list(table.columns)}")
        table = table.loc[:, ["Celltype_Str", "Cytokine", "Gene_Human", "Avg_log2FC"]].copy()
        table.columns = ["cell_type", "cytokine", "gene", "log2fc"]
        table["gene"] = table["gene"].astype(str).str.strip().str.upper()
        table["log2fc"] = pd.to_numeric(table["log2fc"], errors="coerce")
        table = table[table["gene"].notna() & table["gene"].ne("") & table["gene"].ne("--")].dropna(subset=["log2fc"])
        table["direction"] = table["log2fc"].apply(lambda value: "up" if value >= 0 else "down")
        table["weight"] = table["log2fc"].abs()
        tables.append(table[["cytokine", "cell_type", "gene", "direction", "weight"]])

    result = pd.concat(tables, ignore_index=True).drop_duplicates()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False, encoding="utf-8")
    print(
        f"Wrote {len(result):,} genes across {result.cell_type.nunique()} cell types and "
        f"{result.cytokine.nunique()} cytokines to {args.output}"
    )


if __name__ == "__main__":
    main()
