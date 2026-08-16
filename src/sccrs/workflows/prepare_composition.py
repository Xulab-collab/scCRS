#!/usr/bin/env python3
"""Create patient_celltype_counts.csv directly from an annotated AnnData file."""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, help="Annotated .h5ad file")
    parser.add_argument("--out", required=True, help="Output patient_celltype_counts.csv")
    parser.add_argument("--patient-col", default="patient_id")
    parser.add_argument("--celltype-col", default="cell_type")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.h5ad, backed="r")
    required = [args.patient_col, args.celltype_col]
    missing = [column for column in required if column not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"Missing .obs columns: {missing}. Available columns: {list(adata.obs.columns)}"
        )

    obs = adata.obs[required].copy()
    obs.columns = ["patient_id", "cell_type"]
    obs = obs.dropna().astype(str)
    counts = (
        obs.groupby(["patient_id", "cell_type"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values(["patient_id", "cell_type"])
    )
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output, index=False)
    print(f"Wrote {len(counts)} patient-by-cell-type rows for {counts.patient_id.nunique()} patients to {output}")


if __name__ == "__main__":
    main()
