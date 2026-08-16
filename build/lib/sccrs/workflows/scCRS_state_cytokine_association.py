#!/usr/bin/env python3
"""Score T/B-cell states and correlate them with scCRS cytokine-response scores.

This is an exploratory association workflow. It uses patient-by-cell-type
pseudobulk profiles, avoiding cell-level pseudoreplication.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import re
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr

def read_csv(path: str) -> pd.DataFrame:
    err = None
    for enc in ("utf-8", "utf-8-sig", "gb18030", "gbk"):
        try: return pd.read_csv(path, encoding=enc)
        except UnicodeDecodeError as e: err = e
    raise err

def lineage(cell_type: str) -> str | None:
    x = str(cell_type).lower()
    if "b cell" in x or x in {"b", "b cells"}: return "B"
    if "nk" in x or "myeloid" in x: return None
    if "t" in x or "mait" in x: return "T"
    return None

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--cytokine-scores", required=True, help="cell_type_cytokine_scores_with_group.csv")
    p.add_argument("--state-signatures", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--patient-col", default="patient_id")
    p.add_argument("--celltype-col", default="cell_type")
    p.add_argument("--layer", default=None, help="Use raw for adata.raw; otherwise a named layer or omit for X")
    p.add_argument("--min-cells", type=int, default=20)
    p.add_argument("--min-genes", type=int, default=3)
    p.add_argument("--cytokine-score-col", default="raw_rank_score")
    args = p.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    sig = read_csv(args.state_signatures)
    needed = {"lineage", "state", "gene"}
    if not needed.issubset(sig.columns): raise ValueError(f"Signature file needs {needed}")
    sig["gene"] = sig.gene.astype(str).str.upper()
    a = ad.read_h5ad(args.h5ad)
    if args.patient_col not in a.obs or args.celltype_col not in a.obs: raise ValueError("patient/cell-type columns absent from adata.obs")
    if args.layer == "raw":
        if a.raw is None: raise ValueError("adata.raw is absent")
        X, genes = a.raw.X, a.raw.var_names
    elif args.layer:
        X, genes = a.layers[args.layer], a.var_names
    else: X, genes = a.X, a.var_names
    obs = a.obs[[args.patient_col, args.celltype_col]].copy()
    obs.columns = ["patient_id", "cell_type"]
    obs["lineage"] = obs.cell_type.map(lineage)
    obs = obs[obs.lineage.notna()].copy()
    cell_idx = obs.index.to_numpy()
    # AnnData observations are normally a RangeIndex; resolve row positions safely.
    positions = a.obs_names.get_indexer(cell_idx)
    codes, groups = pd.factorize(pd.MultiIndex.from_frame(obs[["patient_id", "cell_type", "lineage"]].astype(str)), sort=True)
    member = sparse.csr_matrix((np.ones(len(codes)), (codes, positions)), shape=(len(groups), a.n_obs))
    n_cells = np.asarray(member.sum(axis=1)).ravel().astype(int)
    pseudo = member @ X
    pseudo = sparse.diags(1 / n_cells) @ pseudo if sparse.issparse(pseudo) else pseudo / n_cells[:, None]
    meta = groups.to_frame(index=False); meta.columns = ["patient_id", "cell_type", "lineage"]; meta["n_cells"] = n_cells
    gene_index = {g:i for i,g in enumerate(pd.Index(genes).astype(str).str.upper())}
    rows = []
    for i, m in meta.iterrows():
        if m.n_cells < args.min_cells: continue
        values = pseudo.getrow(i).toarray().ravel() if sparse.issparse(pseudo) else np.asarray(pseudo[i]).ravel()
        ranks = rankdata(values, method="average") / (len(values) + 1.0)
        for state, tab in sig[sig.lineage.eq(m.lineage)].groupby("state"):
            idx = [gene_index[g] for g in tab.gene if g in gene_index]
            if len(idx) < args.min_genes: continue
            rows.append(dict(patient_id=m.patient_id, cell_type=m.cell_type, lineage=m.lineage,
                state=state, n_cells=m.n_cells, matched_genes=len(idx), signature_genes=len(tab),
                gene_coverage=len(idx)/len(tab), state_rank_score=float(np.mean(ranks[idx]))))
    states = pd.DataFrame(rows)
    states.to_csv(out / "patient_celltype_T_B_state_scores.csv", index=False)
    cyto = read_csv(args.cytokine_scores)
    required = {"patient_id", "cell_type", "cytokine", args.cytokine_score_col}
    if not required.issubset(cyto.columns): raise ValueError(f"Cytokine table needs {required}")
    keep = ["patient_id", "cell_type", "cytokine", args.cytokine_score_col] + [x for x in ["group"] if x in cyto]
    merged = states.merge(cyto[keep], on=["patient_id", "cell_type"], how="inner")
    merged.to_csv(out / "state_cytokine_merged_data.csv", index=False)
    result=[]
    for key, t in merged.groupby(["lineage", "cell_type", "state", "cytokine"], sort=True):
        t=t.dropna(subset=["state_rank_score", args.cytokine_score_col])
        if len(t)<5 or t.state_rank_score.nunique()<2 or t[args.cytokine_score_col].nunique()<2: continue
        rho, pv=spearmanr(t.state_rank_score, t[args.cytokine_score_col])
        result.append(dict(lineage=key[0], cell_type=key[1], state=key[2], cytokine=key[3], n_patients=len(t), spearman_rho=rho, p_value=pv))
    cor=pd.DataFrame(result)
    if len(cor): cor["p_adjust_bh"] = pd.Series(cor.p_value).rank(method="first").to_numpy()/len(cor)*cor.p_value.rank(method="first").to_numpy()*0 + pd.NA
    # Correct BH with pandas-free implementation to preserve monotonicity.
    if len(cor):
        order=np.argsort(cor.p_value.to_numpy()); q=np.empty(len(cor)); q[order]=np.minimum.accumulate((cor.p_value.to_numpy()[order]*len(cor)/np.arange(1,len(cor)+1))[::-1])[::-1]; cor["p_adjust_bh"]=np.minimum(q,1)
    cor.to_csv(out / "T_B_state_cytokine_spearman_correlations.csv", index=False)
    print(f"Wrote {len(states)} state scores and {len(cor)} correlations to {out}")
if __name__ == "__main__": main()
