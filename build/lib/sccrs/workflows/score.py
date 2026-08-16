#!/usr/bin/env python3
"""Patient-level cytokine-response scoring from single-cell RNA-seq.

Research-use prototype based on Cui et al., Nature 2024, Immune Dictionary.
It scores downstream transcriptional-response similarity, not cytokine plasma
concentration or causal cytokine activity.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, rankdata


def clean_label(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def read_table(path: Path) -> pd.DataFrame:
    """Read Excel/CSV, including CSVs saved by Chinese Windows Excel."""
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    decode_error = None
    for encoding in ("utf-8", "utf-8-sig", "gb18030", "gbk"):
        try:
            return pd.read_csv(path, encoding=encoding)
        except UnicodeDecodeError as error:
            decode_error = error
    raise UnicodeDecodeError(
        decode_error.encoding, decode_error.object, decode_error.start, decode_error.end,
        f"Could not decode {path} as UTF-8, UTF-8-BOM, GB18030, or GBK: {decode_error.reason}",
    )


def choose_column(columns: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    names = list(columns)
    normalized = {clean_label(c): c for c in names}
    for candidate in candidates:
        if clean_label(candidate) in normalized:
            return normalized[clean_label(candidate)]
    for candidate in candidates:
        hits = [c for c in names if clean_label(candidate) in clean_label(c)]
        if len(hits) == 1:
            return hits[0]
    return None


def load_signatures(path: Path) -> pd.DataFrame:
    source = read_table(path)
    cytokine = choose_column(source.columns, ["cytokine", "ligand", "stimulation", "condition"])
    cell_type = choose_column(source.columns, ["cell_type", "celltype", "cell type"])
    human_gene = choose_column(source.columns, ["human gene symbol", "human_gene_symbol", "human gene"])
    gene = human_gene or choose_column(source.columns, ["gene_symbol", "gene symbol", "gene", "symbol"])
    direction = choose_column(source.columns, ["direction", "regulation", "sign"])
    logfc = choose_column(source.columns, ["avg_log2fc", "log2fc", "logfc", "log fold change"])
    weight = choose_column(source.columns, ["weight", "abs_log2fc", "effect_size"])
    missing = [name for name, column in {"cytokine": cytokine, "cell_type": cell_type, "gene": gene}.items() if column is None]
    if missing:
        raise ValueError(
            f"Could not identify {missing} in {path.name}. Found columns: {list(source.columns)}. "
            "Use the standardized CSV format documented in example_signature_format.csv."
        )
    if direction is None and logfc is None:
        raise ValueError("Signature table needs a direction (up/down) or log2FC column.")
    result = pd.DataFrame({
        "cytokine": source[cytokine].astype(str).str.strip(),
        "cell_type": source[cell_type].astype(str).str.strip(),
        "gene": source[gene].astype(str).str.strip().str.upper(),
    })
    if direction is not None:
        raw = source[direction].astype(str).str.strip().str.lower()
        result["direction"] = raw.map({
            "up": "up", "upregulated": "up", "positive": "up", "+": "up", "1": "up",
            "down": "down", "downregulated": "down", "negative": "down", "-": "down", "-1": "down",
        })
    else:
        result["direction"] = np.where(pd.to_numeric(source[logfc], errors="coerce") >= 0, "up", "down")
    if weight is not None:
        result["weight"] = pd.to_numeric(source[weight], errors="coerce").abs().fillna(1.0)
    elif logfc is not None:
        result["weight"] = pd.to_numeric(source[logfc], errors="coerce").abs().fillna(1.0)
    else:
        result["weight"] = 1.0
    result = result[result.direction.isin(["up", "down"])].copy()
    result = result.replace({"": np.nan, "nan": np.nan}).dropna(subset=["cytokine", "cell_type", "gene"])
    result = result[result.gene.str.match(r"^[A-Z0-9._-]+$")]
    result["signature_cell_type"] = result.cell_type.map(clean_label)
    return result.drop_duplicates(["cytokine", "signature_cell_type", "gene", "direction"])


def read_cell_type_map(path: Optional[Path]) -> Dict[str, str]:
    if path is None:
        return {}
    table = read_table(path)
    source = choose_column(table.columns, ["patient_cell_type", "source", "from"])
    target = choose_column(table.columns, ["signature_cell_type", "target", "to"])
    if source is None or target is None:
        raise ValueError("Cell-type map must contain patient_cell_type and signature_cell_type columns.")
    return {clean_label(a): clean_label(b) for a, b in zip(table[source], table[target])}


def matrix_and_names(adata: ad.AnnData, layer: Optional[str]):
    if layer is None or layer.lower() == "x":
        return adata.X, adata.var_names
    if layer.lower() == "raw":
        if adata.raw is None:
            raise ValueError("--layer raw was requested but adata.raw is empty.")
        return adata.raw.X, adata.raw.var_names
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' is absent. Available layers: {list(adata.layers.keys())}")
    return adata.layers[layer], adata.var_names


def create_pseudobulk(adata: ad.AnnData, patient_col: str, celltype_col: str, signature_celltype: pd.Series, layer: Optional[str]):
    obs = adata.obs.copy()
    obs["_signature_cell_type"] = signature_celltype.values
    labels = obs[[patient_col, celltype_col, "_signature_cell_type"]].astype(str)
    codes, values = pd.factorize(pd.MultiIndex.from_frame(labels), sort=True)
    membership = sparse.csr_matrix((np.ones(adata.n_obs), (codes, np.arange(adata.n_obs))), shape=(len(values), adata.n_obs))
    counts = np.asarray(membership.sum(axis=1)).ravel().astype(int)
    expression, var_names = matrix_and_names(adata, layer)
    pseudo = membership @ expression
    pseudo = sparse.diags(1.0 / counts) @ pseudo if sparse.issparse(pseudo) else pseudo / counts[:, None]
    groups = values.to_frame(index=False)
    groups.columns = ["patient_id", "cell_type", "signature_cell_type"]
    groups["n_cells"] = counts
    return pseudo, groups, var_names


def make_signature_lookup(signatures: pd.DataFrame, var_names: pd.Index, min_genes: int):
    gene_to_index = {}
    for i, gene in enumerate(pd.Index(var_names).astype(str).str.upper()):
        gene_to_index.setdefault(gene, i)
    lookup, report = {}, []
    for (cell_type, cytokine), table in signatures.groupby(["signature_cell_type", "cytokine"], sort=True):
        entry = {}
        for direction in ("up", "down"):
            subset = table[table.direction == direction]
            pairs = [(gene_to_index[g], float(w)) for g, w in zip(subset.gene, subset.weight) if g in gene_to_index]
            entry[direction] = (
                np.asarray([i for i, _ in pairs], dtype=int),
                np.asarray([w for _, w in pairs], dtype=float),
                len(subset),
            )
        available = len(entry["up"][0]) + len(entry["down"][0])
        report.append({
            "signature_cell_type": cell_type, "cytokine": cytokine,
            "available_genes": available, "signature_genes": entry["up"][2] + entry["down"][2],
            "usable": available >= min_genes and len(entry["up"][0]) > 0,
        })
        lookup[(cell_type, cytokine)] = entry
    return lookup, pd.DataFrame(report)


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.average(values, weights=weights)) if len(values) else np.nan


def score_groups(pseudo, groups: pd.DataFrame, lookup: dict, min_cells: int, min_genes: int) -> pd.DataFrame:
    rows = []
    for index, meta in groups.iterrows():
        if meta.n_cells < min_cells:
            continue
        values = pseudo.getrow(index).toarray().ravel() if sparse.issparse(pseudo) else np.asarray(pseudo[index]).ravel()
        ranks = rankdata(values, method="average") / (len(values) + 1.0)
        for (signature_ct, cytokine), entry in lookup.items():
            if signature_ct != meta.signature_cell_type:
                continue
            up_i, up_w, up_total = entry["up"]
            down_i, down_w, down_total = entry["down"]
            available = len(up_i) + len(down_i)
            if len(up_i) == 0 or available < min_genes:
                continue
            up_score = weighted_mean(ranks[up_i], up_w)
            down_score = weighted_mean(ranks[down_i], down_w) if len(down_i) else 0.5
            coverage = available / max(up_total + down_total, 1)
            rows.append({
                "patient_id": meta.patient_id, "cell_type": meta.cell_type,
                "signature_cell_type": meta.signature_cell_type, "cytokine": cytokine,
                "n_cells": meta.n_cells, "raw_rank_score": 2.0 * (up_score - down_score),
                "up_genes_matched": len(up_i), "down_genes_matched": len(down_i),
                "gene_coverage": coverage, "data_quality": min(1.0, meta.n_cells / 100.0) * coverage,
            })
    return pd.DataFrame(rows)


def patient_metadata(obs: pd.DataFrame, patient_col: str, control_col: Optional[str]) -> pd.DataFrame:
    patients = pd.DataFrame({"patient_id": obs[patient_col].astype(str).unique()})
    if control_col is None:
        return patients
    check = obs.groupby(patient_col, observed=True)[control_col].nunique()
    if (check > 1).any():
        raise ValueError(f"Each patient must have one {control_col} value; inconsistent IDs: {check[check > 1].index.tolist()[:5]}")
    labels = obs.groupby(patient_col, observed=True)[control_col].first().astype(str).reset_index()
    labels.columns = ["patient_id", "cohort"]
    labels.patient_id = labels.patient_id.astype(str)
    return patients.merge(labels, on="patient_id", how="left")


def calibrate_to_controls(scores: pd.DataFrame, metadata: pd.DataFrame, control_label: Optional[str], min_controls: int):
    result = scores.merge(metadata, on="patient_id", how="left")
    result["z_score"] = np.nan
    result["score_0_100"] = 50.0 + 50.0 * result.raw_rank_score.clip(-1, 1)
    result["control_calibrated"] = False
    if control_label is None or "cohort" not in result:
        return result
    is_control = result.cohort.astype(str).eq(str(control_label))
    for _, indices in result.groupby(["cell_type", "cytokine"], sort=False).groups.items():
        indices = list(indices)
        controls = result.loc[indices].loc[is_control.loc[indices], "raw_rank_score"].dropna().to_numpy()
        if len(controls) < min_controls:
            continue
        center = float(np.median(controls))
        scale = float(1.4826 * np.median(np.abs(controls - center)))
        if scale < 1e-6:
            scale = float(np.std(controls, ddof=1))
        if scale < 1e-6:
            continue
        z = (result.loc[indices, "raw_rank_score"] - center) / scale
        result.loc[indices, "z_score"] = z
        result.loc[indices, "score_0_100"] = 100.0 * norm.cdf(z.clip(-4, 4))
        result.loc[indices, "control_calibrated"] = True
    return result


def aggregate_patient_scores(cell_scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (patient, cytokine), table in cell_scores.groupby(["patient_id", "cytokine"], sort=True):
        per_type = table.groupby("cell_type", as_index=False).agg(
            score=("score_0_100", "mean"), n_cells=("n_cells", "first"), quality=("data_quality", "mean")
        )
        rows.append({
            "patient_id": patient, "cytokine": cytokine,
            "overall_score_0_100": float(per_type.score.mean()),
            "composition_weighted_score_0_100": float(np.average(per_type.score, weights=per_type.n_cells)),
            "n_cell_types": len(per_type), "mean_data_quality": float(per_type.quality.mean()),
        })
    result = pd.DataFrame(rows)
    result["patient_rank"] = result.groupby("patient_id").overall_score_0_100.rank(ascending=False, method="min").astype(int)
    return result.sort_values(["patient_id", "patient_rank", "cytokine"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--signatures", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--patient-col", default="patient_id")
    parser.add_argument("--celltype-col", default="cell_type")
    parser.add_argument("--celltype-map", type=Path)
    parser.add_argument("--layer", default="X", help="X (default), raw, or a named AnnData layer.")
    parser.add_argument("--control-col")
    parser.add_argument("--control-label")
    parser.add_argument("--min-cells", type=int, default=30)
    parser.add_argument("--min-genes", type=int, default=5)
    parser.add_argument("--min-controls", type=int, default=5)
    args = parser.parse_args()
    if (args.control_col is None) != (args.control_label is None):
        parser.error("Provide --control-col and --control-label together, or neither.")
    adata = ad.read_h5ad(args.h5ad)
    required = [args.patient_col, args.celltype_col] + ([args.control_col] if args.control_col else [])
    missing = [col for col in required if col not in adata.obs]
    if missing:
        parser.error(f"AnnData .obs is missing {missing}. Available columns: {list(adata.obs.columns)}")
    if adata.n_obs == 0 or adata.n_vars == 0:
        parser.error("AnnData cannot be empty.")
    args.outdir.mkdir(parents=True, exist_ok=True)
    signatures = load_signatures(args.signatures)
    mapping = read_cell_type_map(args.celltype_map)
    # AnnData .obs cell-type columns are often pandas Categoricals.  Convert
    # them to ordinary strings before mapping, because fillna/where cannot
    # insert labels absent from the original categorical levels.
    raw_cell_types = adata.obs[args.celltype_col].astype(str).map(clean_label).astype(object)
    mapped_cell_types = raw_cell_types.map(mapping).astype(object)
    signature_cell_types = mapped_cell_types.where(mapped_cell_types.notna(), raw_cell_types).astype(str)
    pseudo, groups, gene_names = create_pseudobulk(adata, args.patient_col, args.celltype_col, signature_cell_types, args.layer)
    lookup, coverage = make_signature_lookup(signatures, gene_names, args.min_genes)
    coverage.to_csv(args.outdir / "signature_gene_coverage.csv", index=False)
    cell_scores = score_groups(pseudo, groups, lookup, args.min_cells, args.min_genes)
    if cell_scores.empty:
        raise RuntimeError("No score could be calculated. Check cell-type labels, gene symbols, --layer, --min-cells and --min-genes.")
    metadata = patient_metadata(adata.obs, args.patient_col, args.control_col)
    cell_scores = calibrate_to_controls(cell_scores, metadata, args.control_label, args.min_controls)
    patient_scores = aggregate_patient_scores(cell_scores)
    cell_scores.sort_values(["patient_id", "cell_type", "cytokine"]).to_csv(args.outdir / "cell_type_cytokine_scores.csv", index=False)
    patient_scores.to_csv(args.outdir / "patient_cytokine_scores.csv", index=False)
    groups.to_csv(args.outdir / "patient_celltype_counts.csv", index=False)
    print(f"Wrote {len(cell_scores):,} cell-type cytokine scores and {len(patient_scores):,} patient cytokine scores to {args.outdir}")


if __name__ == "__main__":
    main()
