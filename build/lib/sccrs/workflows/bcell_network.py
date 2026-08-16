#!/usr/bin/env python3
"""Build cytokine-response regulatory networks for every B-cell subtype and group.

For each patient 闂?B-cell-subtype pseudobulk, this workflow calculates a
scCRS rank-based cytokine-response score from the supplied cytokine dictionary.
Within each group, it links score-associated target genes to literature-backed
CollecTRI/DoRothEA TF-target edges.  Between groups, it reports the difference
in the target-gene/response Spearman correlation for every prior edge.

Networks are hypothesis-generating: a high score denotes transcriptional
similarity to an in-vitro cytokine-response signature, not measured cytokine
abundance or proof of a direct in-vivo regulatory interaction.
"""
from __future__ import annotations

import argparse
import itertools
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages
import networkx as nx
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, rankdata, spearmanr


def bh(values: np.ndarray) -> np.ndarray:
    result = np.full(len(values), np.nan, dtype=float)
    keep = np.isfinite(values)
    if not keep.any():
        return result
    x = values[keep]
    order = np.argsort(x)
    adjusted = x[order] * len(x) / np.arange(1, len(x) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty(len(x)); restored[order] = np.minimum(adjusted, 1.0)
    result[np.where(keep)[0][order]] = restored
    return result


def safe_name(x: object) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(x)).strip("_") or "unnamed"


def read_table(path: str) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    last_error = None
    for encoding in ("utf-8", "utf-8-sig", "gb18030", "gbk"):
        try:
            return pd.read_csv(path, sep=None, engine="python", encoding=encoding)
        except (UnicodeDecodeError, pd.errors.ParserError) as exc:
            last_error = exc
    raise last_error


def select_matrix(adata: ad.AnnData, layer: str | None):
    if layer == "raw":
        if adata.raw is None:
            raise ValueError("--layer raw requested but .raw is absent.")
        return adata.raw.X, pd.Index(adata.raw.var_names).astype(str)
    if layer:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' is not present in .layers.")
        return adata.layers[layer], pd.Index(adata.var_names).astype(str)
    return adata.X, pd.Index(adata.var_names).astype(str)


def is_b_cell(label: object) -> bool:
    text = str(label).lower()
    return bool(re.search(r"(^|[ _-])b([ _-]|$)|b[ _-]?cell", text)) or any(token in text for token in ("plasma", "plasmablast", "antibody-secreting", "asc"))


def load_signatures(path: str, signature_celltype: str, cytokines: str | None) -> pd.DataFrame:
    sig = read_table(path)
    required = {"cytokine", "cell_type", "gene", "direction"}
    missing = required.difference(sig.columns)
    if missing:
        raise ValueError(f"Signature table lacks columns {sorted(missing)}; found {list(sig.columns)}")
    sig = sig.copy()
    sig["cytokine"] = sig["cytokine"].astype(str)
    sig["cell_type"] = sig["cell_type"].astype(str)
    sig["gene"] = sig["gene"].astype(str).str.upper().str.strip()
    sig["direction"] = sig["direction"].astype(str).str.lower().str.strip()
    sig["weight"] = pd.to_numeric(sig.get("weight", 1.0), errors="coerce").fillna(1.0).abs()
    sig = sig.loc[(sig["cell_type"].str.lower() == signature_celltype.lower()) & sig["direction"].isin(["up", "down"]) & sig["gene"].ne("")]
    if cytokines:
        wanted = {x.strip() for x in cytokines.split(",") if x.strip()}
        sig = sig.loc[sig["cytokine"].isin(wanted)]
        absent = wanted.difference(sig["cytokine"].unique())
        if absent:
            print(f"Warning: cytokines absent from the selected signature cell type: {sorted(absent)}")
    if sig.empty:
        raise ValueError(f"No usable signatures for cell_type='{signature_celltype}'.")
    return sig.drop_duplicates(["cytokine", "gene", "direction"])


def load_prior(path: str) -> pd.DataFrame:
    prior = read_table(path)
    source = next((c for c in ["source_genesymbol", "source", "tf", "TF"] if c in prior.columns), None)
    target = next((c for c in ["target_genesymbol", "target", "gene", "Target"] if c in prior.columns), None)
    if source is None or target is None:
        raise ValueError(f"TF-target prior needs source and target columns; found {list(prior.columns)}")
    source_values = prior.loc[:, source]
    target_values = prior.loc[:, target]
    if isinstance(source_values, pd.DataFrame): source_values = source_values.iloc[:, 0]
    if isinstance(target_values, pd.DataFrame): target_values = target_values.iloc[:, 0]
    out = prior.loc[:, ~prior.columns.duplicated()].copy().drop(columns=["tf", "target"], errors="ignore")
    out["tf"] = source_values.astype(str).str.upper().to_numpy()
    out["target"] = target_values.astype(str).str.upper().to_numpy()
    out = out.loc[~out["tf"].str.contains("COMPLEX|_", regex=True, na=False)]
    for col in ("sources", "references", "n_references"):
        if col not in out: out[col] = np.nan
    return out.drop_duplicates(["tf", "target"])


def patient_subtype_pseudobulk(adata, patient_col, celltype_col, group_col, layer, min_cells):
    obs = adata.obs.copy()
    selected = obs[celltype_col].map(is_b_cell).to_numpy()
    if not selected.any():
        raise ValueError("No B-cell labels detected. Use B-cell labels such as 'naive B cell' or 'B_Transitional'.")
    obs = obs.loc[selected, [patient_col, celltype_col, group_col]].copy()
    obs[patient_col] = obs[patient_col].astype(str)
    obs[celltype_col] = obs[celltype_col].astype(str)
    obs[group_col] = obs[group_col].astype(str)
    key = obs[patient_col] + "\x1f" + obs[celltype_col]
    codes, names = pd.factorize(key, sort=True)
    positions = np.flatnonzero(selected)
    member = sparse.csr_matrix((np.ones(len(codes)), (codes, positions)), shape=(len(names), adata.n_obs))
    n_cells = np.asarray(member.sum(axis=1)).ravel().astype(int)
    X, genes = select_matrix(adata, layer)
    pseudo = member @ X
    pseudo = sparse.diags(1 / n_cells) @ pseudo if sparse.issparse(pseudo) else pseudo / n_cells[:, None]
    meta = obs.groupby(key, sort=True).agg(**{
        "patient_id": (patient_col, "first"), "cell_type": (celltype_col, "first"), "group": (group_col, lambda x: x.mode().iat[0]),
    }).reindex(names).reset_index(drop=True)
    meta["n_cells"] = n_cells
    keep = n_cells >= min_cells
    return pseudo[keep], meta.loc[keep].reset_index(drop=True), genes


def score_signature(values: np.ndarray, ranks: np.ndarray, gene_index: dict[str, int], signature: pd.DataFrame):
    up = signature.loc[signature.direction.eq("up")]
    down = signature.loc[signature.direction.eq("down")]
    def weighted(frame):
        idx = [gene_index[g] for g in frame.gene if g in gene_index]
        weights = frame.loc[frame.gene.isin(gene_index), "weight"].to_numpy(float)
        if not idx: return np.full(values.shape[0], np.nan), 0
        return np.average(ranks[:, idx], axis=1, weights=weights), len(idx)
    up_score, n_up = weighted(up); down_score, n_down = weighted(down)
    if n_up == 0: return np.full(values.shape[0], np.nan), 0, 0
    if n_down == 0: down_score = np.full(values.shape[0], .5)
    return 2 * (up_score - down_score), n_up, n_down


def gene_associations(values: np.ndarray, genes: np.ndarray, score: np.ndarray) -> pd.DataFrame:
    rows = []
    for j, gene in enumerate(genes):
        x = values[:, j]
        if np.unique(x).size < 2 or np.unique(score).size < 2:
            rho, p = np.nan, np.nan
        else:
            rho, p = spearmanr(x, score)
        rows.append((gene, rho, p, float(np.mean(x))))
    result = pd.DataFrame(rows, columns=["gene", "rho", "p_value", "mean_expression"])
    result["q_value"] = bh(result["p_value"].to_numpy(float))
    return result


def make_edges(prior: pd.DataFrame, stats: pd.DataFrame, cytokine: str, celltype: str, group: str, min_abs_rho: float, max_q: float, max_edges: int):
    targets = stats.rename(columns={"gene": "target", "rho": "target_rho", "p_value": "target_p_value", "q_value": "target_q_value", "mean_expression": "target_mean_expression"})
    tfs = stats.rename(columns={"gene": "tf", "rho": "tf_rho", "p_value": "tf_p_value", "q_value": "tf_q_value", "mean_expression": "tf_mean_expression"})
    edges = prior.merge(targets, on="target", how="inner").merge(tfs[["tf", "tf_rho", "tf_q_value", "tf_mean_expression"]], on="tf", how="inner")
    edges = edges.loc[edges["target_rho"].abs().ge(min_abs_rho) & edges["target_q_value"].le(max_q)].copy()
    edges["cytokine"] = cytokine; edges["cell_type"] = celltype; edges["group"] = group
    references = pd.to_numeric(edges.get("n_references", 0), errors="coerce").fillna(0).clip(upper=5)
    edges["edge_priority"] = edges["target_rho"].abs() * (1 + references / 5)
    return edges.sort_values("edge_priority", ascending=False).head(max_edges).reset_index(drop=True)


def plot_network(edges: pd.DataFrame, title: str, pdf: PdfPages, difference: bool = False, max_plot_edges: int = 35):
    """Draw a deliberately sparse, layered, Illustrator-editable network page.

    Full edge tables are retained in CSV/GraphML.  Restricting the displayed
    edges avoids the unreadable label and edge overlap that occurs in dense
    regulatory networks.
    """
    if edges.empty:
        return
    rank_col = "delta_rho" if difference else "edge_priority"
    if rank_col not in edges:
        rank_col = "target_rho"
    shown = edges.loc[edges[rank_col].abs().nlargest(max_plot_edges).index].copy()
    graph = nx.DiGraph()
    for row in shown.itertuples(index=False):
        value = float(getattr(row, "target_rho", getattr(row, "delta_rho", 0.0)))
        graph.add_edge(row.tf, row.target, weight=max(abs(value), 0.05))
    if not graph.nodes:
        return
    tfs = sorted({r.tf for r in shown.itertuples(index=False)})
    targets = sorted(set(graph.nodes).difference(tfs))
    # Graphviz dot provides non-overlapping, left-to-right TF -> target layers.
    # A deterministic bipartite fallback keeps the script dependency-free.
    graph.graph.update(rankdir="LR", nodesep="0.55", ranksep="2.0")
    try:
        pos = nx.nx_pydot.graphviz_layout(graph, prog="dot")
    except Exception:
        pos = {}
        for i, node in enumerate(tfs): pos[node] = (0, -i)
        for i, node in enumerate(targets): pos[node] = (3, -i)
    height = max(7.5, 0.30 * max(len(tfs), len(targets), 1))
    fig, ax = plt.subplots(figsize=(13.5, height)); ax.axis("off")
    if difference:
        values = shown.groupby("target")["delta_rho"].median().to_dict()
        colors = ["#c53b53" if values.get(n, 0) >= 0 else "#3c8dbc" for n in targets]
        target_label = "Target (red: stronger in comparison; blue: stronger in reference)"
    else:
        values = shown.groupby("target")["target_rho"].median().to_dict()
        colors = ["#c53b53" if values.get(n, 0) >= 0 else "#3c8dbc" for n in targets]
        target_label = "Target (red: positive; blue: negative response association)"
    nx.draw_networkx_edges(
        graph, pos, ax=ax, alpha=.38, arrows=True, arrowsize=9,
        width=[0.7 + 2.4 * graph[u][v]["weight"] for u, v in graph.edges],
        edge_color="#555555", connectionstyle="arc3,rad=0.04",
    )
    nx.draw_networkx_nodes(graph, pos, nodelist=tfs, node_shape="s", node_size=540,
                           node_color="#374f8b", edgecolors="white", linewidths=.8, ax=ax, label="TF")
    nx.draw_networkx_nodes(graph, pos, nodelist=targets, node_size=260,
                           node_color=colors, edgecolors="white", linewidths=.7, ax=ax, label=target_label)
    nx.draw_networkx_labels(graph, pos, font_size=7.2, font_family="DejaVu Sans", ax=ax,
                            bbox={"facecolor": "white", "alpha": .72, "edgecolor": "none", "pad": .15})
    ax.set_title(f"{title}\nDisplayed top {len(shown)} edges; complete network is provided in CSV/GraphML", fontsize=12, fontweight="bold", pad=14)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout(pad=1.0); pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
def differential_edges(reference: pd.DataFrame, comparison: pd.DataFrame, group_a: str, group_b: str, n_a: int, n_b: int):
    a = reference[["tf", "target", "target_rho", "target_q_value"]].rename(columns={"target_rho": "rho_reference", "target_q_value": "q_reference"})
    b = comparison[["tf", "target", "target_rho", "target_q_value"]].rename(columns={"target_rho": "rho_comparison", "target_q_value": "q_comparison"})
    merged = a.merge(b, on=["tf", "target"], how="outer")
    merged[["rho_reference", "rho_comparison"]] = merged[["rho_reference", "rho_comparison"]].fillna(0.0)
    merged["delta_rho"] = merged["rho_comparison"] - merged["rho_reference"]
    # Fisher approximation is used only as an edge-ranking/comparison statistic.
    za = np.arctanh(merged["rho_reference"].clip(-.999, .999)); zb = np.arctanh(merged["rho_comparison"].clip(-.999, .999))
    se = np.sqrt(1 / (n_a - 3) + 1 / (n_b - 3))
    merged["delta_z"] = (zb - za) / se
    merged["delta_p_value"] = 2 * norm.sf(np.abs(merged["delta_z"]))
    merged["delta_q_value"] = bh(merged["delta_p_value"].to_numpy(float))
    merged["reference_group"] = group_a; merged["comparison_group"] = group_b
    return merged.sort_values("delta_rho", key=lambda x: x.abs(), ascending=False)


def main():
    # Keep all plot elements vector-based and retain editable TrueType text in Adobe Illustrator.
    rcParams["pdf.fonttype"] = 42
    rcParams["ps.fonttype"] = 42
    rcParams["pdf.use14corefonts"] = False
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--h5ad", required=True); parser.add_argument("--signatures", required=True)
    parser.add_argument("--tf-target-prior", required=True); parser.add_argument("--outdir", required=True)
    parser.add_argument("--patient-col", default="patient_id"); parser.add_argument("--celltype-col", default="cell_type"); parser.add_argument("--group-col", default="group")
    parser.add_argument("--layer", default="raw"); parser.add_argument("--signature-celltype", default="B_cell", help="Dictionary cell_type to apply to B-cell subtypes.")
    parser.add_argument("--cytokines", default=None, help="Comma-separated subset; default analyses every cytokine in the B_cell dictionary.")
    parser.add_argument("--groups", default=None, help="Comma-separated groups to retain; default retains all.")
    parser.add_argument("--pooled", action="store_true", help="Pool all retained patients for one network per cytokine and B-cell subtype; disables group comparisons.")
    parser.add_argument("--compare", default=None, help="Optional pairs such as 'Normal:HGG,Non-HGG:HGG'; default compares all possible pairs.")
    parser.add_argument("--min-cells", type=int, default=30); parser.add_argument("--min-patients", type=int, default=5)
    parser.add_argument("--min-abs-rho", type=float, default=.50); parser.add_argument("--max-q", type=float, default=.10); parser.add_argument("--max-edges", type=int, default=80)
    parser.add_argument("--make-pdf", action="store_true", help="Write multi-page PDF network atlases; use --cytokines first for a focused inspection.")
    parser.add_argument("--plot-max-edges", type=int, default=35, help="Maximum edges shown on each PDF page; all retained edges remain in CSV/GraphML.")
    args = parser.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad)
    for col in (args.patient_col, args.celltype_col, args.group_col):
        if col not in adata.obs: raise ValueError(f"Missing .obs column '{col}'. Available: {list(adata.obs.columns)}")
    sig = load_signatures(args.signatures, args.signature_celltype, args.cytokines)
    prior = load_prior(args.tf_target_prior)
    pseudo, meta, gene_names = patient_subtype_pseudobulk(adata, args.patient_col, args.celltype_col, args.group_col, args.layer, args.min_cells)
    values = pseudo.toarray() if sparse.issparse(pseudo) else np.asarray(pseudo)
    genes = gene_names.astype(str).str.upper().to_numpy(); gene_index = {g: i for i, g in enumerate(genes)}
    ranks = np.apply_along_axis(rankdata, 1, values, method="average") / (values.shape[1] + 1.0)
    if args.groups:
        wanted = {x.strip() for x in args.groups.split(",") if x.strip()}; keep = meta["group"].isin(wanted).to_numpy()
        values, ranks, meta = values[keep], ranks[keep], meta.loc[keep].reset_index(drop=True)
    if args.pooled:
        meta = meta.copy()
        meta["original_group"] = meta["group"]
        meta["group"] = "All_patients"
        args.compare = None
    meta.to_csv(out / "patient_Bcell_subtype_pseudobulk_metadata.csv", index=False)
    score_rows, coverage_rows, all_edges, all_stats, insufficient = [], [], [], [], []
    score_matrix = {}
    for cytokine, signature in sig.groupby("cytokine", sort=True):
        score, n_up, n_down = score_signature(values, ranks, gene_index, signature)
        matched = n_up + n_down
        coverage_rows.append({"cytokine": cytokine, "signature_cell_type": args.signature_celltype, "matched_up_genes": n_up, "matched_down_genes": n_down, "total_signature_genes": len(signature), "matched_fraction": matched / len(signature)})
        score_matrix[cytokine] = score
        for i, row in meta.iterrows(): score_rows.append({"pseudobulk_index": i, "patient_id": row.patient_id, "cell_type": row.cell_type, "group": row.group, "n_cells": row.n_cells, "cytokine": cytokine, "response_rank_score": score[i]})
    scores = pd.DataFrame(score_rows); scores.to_csv(out / "patient_Bcell_subtype_cytokine_response_scores.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(out / "cytokine_signature_gene_coverage.csv", index=False)
    # Derive independent group-specific networks.
    for (cytokine, celltype, group), subset in scores.groupby(["cytokine", "cell_type", "group"], sort=True):
        idx = subset["pseudobulk_index"].to_numpy(int)
        if len(idx) < args.min_patients:
            insufficient.append({"cytokine": cytokine, "cell_type": celltype, "group": group, "n_patients": len(idx), "reason": "below_min_patients"}); continue
        stats = gene_associations(values[idx], genes, subset["response_rank_score"].to_numpy(float))
        stats["cytokine"] = cytokine; stats["cell_type"] = celltype; stats["group"] = group; stats["n_patients"] = len(idx)
        all_stats.append(stats)
        edges = make_edges(prior, stats, cytokine, celltype, group, args.min_abs_rho, args.max_q, args.max_edges)
        if not edges.empty: all_edges.append(edges)
    gene_stats = pd.concat(all_stats, ignore_index=True) if all_stats else pd.DataFrame()
    edges = pd.concat(all_edges, ignore_index=True) if all_edges else pd.DataFrame()
    gene_stats.to_csv(out / "all_group_subtype_cytokine_gene_response_associations.csv", index=False)
    edges.to_csv(out / "all_group_subtype_cytokine_TF_target_edges.csv", index=False)
    pd.DataFrame(insufficient).to_csv(out / "networks_not_tested_insufficient_patients.csv", index=False)
    diff_tables = []
    if not edges.empty:
        available = scores.groupby(["cytokine", "cell_type", "group"])["patient_id"].nunique().to_dict()
        for (cytokine, celltype), e in edges.groupby(["cytokine", "cell_type"], sort=True):
            gs = sorted(e["group"].unique())
            requested = [tuple(x.split(":", 1)) for x in args.compare.split(",")] if args.compare else list(itertools.combinations(gs, 2))
            for g1, g2 in requested:
                a, b = e.loc[e.group.eq(g1)], e.loc[e.group.eq(g2)]
                n1, n2 = available.get((cytokine, celltype, g1), 0), available.get((cytokine, celltype, g2), 0)
                if a.empty or b.empty or n1 < args.min_patients or n2 < args.min_patients or n1 <= 3 or n2 <= 3: continue
                d = differential_edges(a, b, g1, g2, n1, n2); d["cytokine"] = cytokine; d["cell_type"] = celltype; d["n_reference"] = n1; d["n_comparison"] = n2
                diff_tables.append(d)
    diffs = pd.concat(diff_tables, ignore_index=True) if diff_tables else pd.DataFrame()
    diffs.to_csv(out / "Bcell_subtype_cytokine_network_group_differences.csv", index=False)
    if args.make_pdf and not edges.empty:
        with PdfPages(out / "Bcell_subtype_cytokine_group_network_atlas.pdf") as pdf:
            for (cytokine, celltype, group), e in edges.groupby(["cytokine", "cell_type", "group"], sort=True):
                plot_network(e, f"{cytokine} response | {celltype} | {group}", pdf, max_plot_edges=args.plot_max_edges)
        if not diffs.empty:
            with PdfPages(out / "Bcell_subtype_cytokine_network_group_difference_atlas.pdf") as pdf:
                for (cytokine, celltype, ref, comp), d in diffs.groupby(["cytokine", "cell_type", "reference_group", "comparison_group"], sort=True):
                    selected = d.loc[d["delta_q_value"].le(args.max_q)].head(args.max_edges)
                    if not selected.empty: plot_network(selected, f"{cytokine} | {celltype} | {comp} vs {ref}", pdf, difference=True, max_plot_edges=args.plot_max_edges)
    print(f"Completed {len(scores)} response scores, {len(edges)} group-specific regulatory edges, and {len(diffs)} differential edges. Results: {out}")

if __name__ == "__main__":
    main()
