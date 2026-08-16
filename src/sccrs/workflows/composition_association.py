#!/usr/bin/env python3
"""Associate patient-level scCRS cytokine scores with cell-type proportions.

Use --cytokine to restrict the heatmap to biologically selected cytokines.
The primary endpoint is overall_score_0_100, which is not cell-count weighted.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def read_csv(path: str) -> pd.DataFrame:
    error = None
    for encoding in ("utf-8", "utf-8-sig", "gb18030", "gbk"):
        try:
            return pd.read_csv(path, encoding=encoding)
        except UnicodeDecodeError as exc:
            error = exc
    raise error  # pragma: no cover


def parse_list(items: list[str] | None) -> list[str]:
    return [name.strip() for item in (items or []) for name in item.split(",") if name.strip()]


def bh(p: np.ndarray) -> np.ndarray:
    answer = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    values = p[valid]
    if not len(values):
        return answer
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    restored = np.empty(len(values)); restored[order] = np.minimum(ranked, 1)
    answer[valid] = restored
    return answer


def stars(q: float) -> str:
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""


def draw_heatmap(rho: pd.DataFrame, q: pd.DataFrame, pdf: PdfPages) -> None:
    rows, cols = rho.shape
    # Equal aspect means every tile is a true square; plot stays compact for small selected panels.
    fig, ax = plt.subplots(figsize=(max(6.6, 0.78 * cols + 3.0), max(4.6, 0.60 * rows + 2.2)))
    image = ax.imshow(rho.to_numpy(float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
    ax.set_xticks(np.arange(cols), rho.columns, rotation=45, ha="right", fontsize=12)
    ax.set_yticks(np.arange(rows), rho.index, fontsize=12)
    ax.set_xlabel("Cell-type proportion", fontsize=13, fontweight="bold")
    ax.set_ylabel("Cytokine response", fontsize=13, fontweight="bold")
    ax.set_title("Patient-level cytokine response versus cell-type proportion", fontsize=14, fontweight="bold", pad=12)
    for i in range(rows):
        for j in range(cols):
            label = stars(q.iloc[i, j]) if np.isfinite(q.iloc[i, j]) else ""
            if label:
                ax.text(j, i, label, ha="center", va="center", fontsize=10, color="black")
    ax.set_xticks(np.arange(-.5, cols, 1), minor=True)
    ax.set_yticks(np.arange(-.5, rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=.6)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(image, ax=ax, fraction=.035, pad=.03)
    cbar.set_label("Spearman correlation (rho)", fontsize=12, fontweight="bold")
    cbar.ax.tick_params(labelsize=11)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    fig.savefig(pdf._file.fh.name.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-scores", required=True)
    parser.add_argument("--celltype-counts", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--patient-score-col", default="overall_score_0_100")
    parser.add_argument("--min-patients", type=int, default=5)
    parser.add_argument("--min-cell-proportion", type=float, default=0.0)
    parser.add_argument("--cytokine", action="append", default=None,
                        help="Heatmap cytokine(s); repeat option or use comma-separated names")
    args = parser.parse_args()
    target_cytokines = parse_list(args.cytokine)

    scores, counts = read_csv(args.patient_scores), read_csv(args.celltype_counts)
    required_scores = {"patient_id", "cytokine", args.patient_score_col}
    required_counts = {"patient_id", "cell_type", "n_cells"}
    if missing := required_scores.difference(scores.columns): raise ValueError(f"Missing score columns: {sorted(missing)}")
    if missing := required_counts.difference(counts.columns): raise ValueError(f"Missing count columns: {sorted(missing)}")
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

    counts = counts[["patient_id", "cell_type", "n_cells"]].dropna().copy()
    counts["patient_id"] = counts.patient_id.astype(str); counts["cell_type"] = counts.cell_type.astype(str)
    counts["n_cells"] = pd.to_numeric(counts.n_cells, errors="coerce")
    counts = counts.dropna().groupby(["patient_id", "cell_type"], as_index=False).n_cells.sum()
    counts["total_cells"] = counts.groupby("patient_id").n_cells.transform("sum")
    counts["cell_proportion"] = counts.n_cells / counts.total_cells
    counts = counts.loc[counts.cell_proportion.groupby(counts.cell_type).transform("mean") >= args.min_cell_proportion]
    counts.to_csv(out / "patient_celltype_proportions.csv", index=False)

    keep = ["patient_id", "cytokine", args.patient_score_col] + (["group"] if "group" in scores else [])
    scores = scores[keep].copy(); scores.patient_id = scores.patient_id.astype(str)
    scores = scores.rename(columns={args.patient_score_col: "cytokine_response_score"})
    merged = scores.merge(counts, on="patient_id", how="inner")
    merged.to_csv(out / "patient_cytokine_celltype_proportion_merged_data.csv", index=False)

    result = []
    for (cytokine, cell_type), tab in merged.groupby(["cytokine", "cell_type"], sort=True):
        tab = tab.dropna(subset=["cytokine_response_score", "cell_proportion"])
        if len(tab) < args.min_patients or tab.cytokine_response_score.nunique() < 2 or tab.cell_proportion.nunique() < 2:
            rho, pvalue = np.nan, np.nan
        else:
            rho, pvalue = spearmanr(tab.cytokine_response_score, tab.cell_proportion)
        result.append(dict(cytokine=cytokine, cell_type=cell_type, n_patients=len(tab), spearman_rho=rho,
                           p_value=pvalue, mean_cell_proportion=tab.cell_proportion.mean()))
    result = pd.DataFrame(result)
    result["p_adjust_bh"] = bh(result.p_value.to_numpy(float)); result["significance"] = result.p_adjust_bh.map(stars)
    result.to_csv(out / "cytokine_celltype_proportion_spearman_correlations.csv", index=False)
    rho = result.pivot(index="cytokine", columns="cell_type", values="spearman_rho")
    q = result.pivot(index="cytokine", columns="cell_type", values="p_adjust_bh").reindex_like(rho)
    if target_cytokines:
        absent = [x for x in target_cytokines if x not in rho.index]
        if absent: raise ValueError(f"Cytokines not found: {absent}; available examples: {sorted(rho.index)[:15]}")
        rho = rho.reindex(target_cytokines); q = q.reindex(index=target_cytokines, columns=rho.columns)
    else:
        rho = rho.loc[rho.abs().max(axis=1).sort_values(ascending=False).index]; q = q.reindex_like(rho)
    rho.to_csv(out / "cytokine_celltype_proportion_rho_matrix.csv"); q.to_csv(out / "cytokine_celltype_proportion_BH_qvalue_matrix.csv")
    with PdfPages(out / "cytokine_response_celltype_proportion_correlation_heatmap.pdf") as pdf: draw_heatmap(rho, q, pdf)
    print(f"Wrote {len(result)} tests and a {rho.shape[0]} cytokine x {rho.shape[1]} cell-type heatmap to {out}")


if __name__ == "__main__": main()
