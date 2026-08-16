#!/usr/bin/env python3
"""Append a patient-level AnnData .obs group column to cytokine-score CSV files.

This does not recalculate cytokine scores or perform control calibration. It
only adds the requested group label to copies of the output tables, making the
results ready for group comparisons and plotting.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd


def get_patient_groups(adata: ad.AnnData, patient_col: str, group_col: str) -> pd.DataFrame:
    if patient_col not in adata.obs or group_col not in adata.obs:
        raise ValueError(
            f"AnnData .obs must contain {patient_col!r} and {group_col!r}. "
            f"Available columns: {list(adata.obs.columns)}"
        )
    metadata = adata.obs[[patient_col, group_col]].copy()
    consistency = metadata.groupby(patient_col, observed=True)[group_col].nunique(dropna=False)
    inconsistent = consistency[consistency != 1]
    if not inconsistent.empty:
        raise ValueError(
            f"Each patient must have exactly one group label. Inconsistent patients: {inconsistent.index.tolist()[:10]}"
        )
    metadata = metadata.drop_duplicates(subset=[patient_col]).astype({patient_col: str})
    metadata.columns = ["patient_id", group_col]
    return metadata


def append_group(path: Path, patient_groups: pd.DataFrame, group_col: str) -> Path:
    table = pd.read_csv(path)
    if "patient_id" not in table.columns:
        raise ValueError(f"{path.name} has no patient_id column.")
    table["patient_id"] = table["patient_id"].astype(str)
    if group_col in table.columns:
        table = table.drop(columns=[group_col])
    result = table.merge(patient_groups, on="patient_id", how="left", validate="many_to_one")
    missing = result[group_col].isna().sum()
    if missing:
        raise ValueError(f"{path.name}: {missing} rows have no matching group label.")
    output = path.with_name(f"{path.stem}_with_{group_col}{path.suffix}")
    result.to_csv(output, index=False)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--patient-col", default="patient_id")
    parser.add_argument("--group-col", default="group")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.h5ad, backed="r")
    patient_groups = get_patient_groups(adata, args.patient_col, args.group_col)
    patient_groups.to_csv(args.outdir / f"patient_groups_from_h5ad_{args.group_col}.csv", index=False)

    tables = [
        "patient_cytokine_scores.csv",
        "cell_type_cytokine_scores.csv",
        "patient_celltype_counts.csv",
    ]
    outputs = []
    for name in tables:
        source = args.outdir / name
        if not source.exists():
            raise FileNotFoundError(f"Expected score table is missing: {source}")
        outputs.append(append_group(source, patient_groups, args.group_col))
    print(f"Wrote patient-group metadata for {len(patient_groups)} patients.")
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
