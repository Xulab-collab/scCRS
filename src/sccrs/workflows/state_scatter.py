#!/usr/bin/env python3
"""Plot selected scCRS cytokine-response versus T/B-cell state associations.

The unit of analysis is one patient within one cell type.  Use the merged
table produced by scCRS_state_cytokine_association_v2.py, rather than
cell-level values, to avoid pseudoreplication.
"""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


GROUP_COLORS = {
    "Normal": "#74ADD1",
    "Non-HGG": "#7FBF7B",
    "Non_High": "#7FBF7B",
    "Low globulin": "#7FBF7B",
    "HGG": "#F46D43",
    "High": "#F46D43",
    "High globulin": "#F46D43",
}


def read_csv(path: str) -> pd.DataFrame:
    last_error = None
    for encoding in ("utf-8", "utf-8-sig", "gb18030", "gbk"):
        try:
            return pd.read_csv(path, encoding=encoding)
        except UnicodeDecodeError as error:
            last_error = error
    raise last_error  # pragma: no cover


def values(argument: list[str]) -> list[str]:
    """Permit either repeated options or comma-separated option values."""
    return [item.strip() for part in argument for item in part.split(",") if item.strip()]


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")


def configure_editable_font(font_family: str, font_file: str | None) -> str:
    """Register a system font with Matplotlib, including fontconfig-only fonts."""
    candidates = [Path(font_file)] if font_file else []
    if not candidates:
        try:
            query = subprocess.run(
                ["fc-match", "-f", "%{file}\n", font_family],
                check=True, capture_output=True, text=True,
            )
            candidates = [Path(line.strip()) for line in query.stdout.splitlines() if line.strip()]
        except (OSError, subprocess.SubprocessError):
            candidates = []
    for candidate in candidates:
        if candidate.is_file():
            font_manager.fontManager.addfont(str(candidate))
            resolved = font_manager.FontProperties(fname=str(candidate)).get_name()
            plt.rcParams.update({
                "font.family": resolved,
                "font.sans-serif": [resolved, font_family, "Arial", "DejaVu Sans"],
                "pdf.fonttype": 42,
                "ps.fonttype": 42,
                "pdf.use14corefonts": False,
                "svg.fonttype": "none",
            })
            print(f"Using editable font: {resolved} ({candidate})")
            return resolved
    plt.rcParams.update({
        "font.family": font_family,
        "font.sans-serif": [font_family, "Arial", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "pdf.use14corefonts": False,
        "svg.fonttype": "none",
    })
    print(f"Warning: could not register '{font_family}'. Use --font-file /absolute/path/to/arial.ttf for an exact font.")
    return font_family


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Draw selected patient-level cytokine-response/state scatter plots."
    )
    parser.add_argument("--data", required=True,
                        help="state_cytokine_merged_data.csv from the association analysis")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--cytokine", action="append", required=True,
                        help="Cytokine name; repeat option or give comma-separated names")
    parser.add_argument("--cell-type", action="append", required=True,
                        help="Cell-type name; repeat option or give comma-separated names")
    parser.add_argument("--state", action="append", required=True,
                        help="State name; repeat option or give comma-separated names")
    parser.add_argument("--cytokine-score-col", default="raw_rank_score")
    parser.add_argument("--color-by", default="group",
                        help="Column used for point colour; set to empty string to disable")
    parser.add_argument("--point-size", type=float, default=72)
    parser.add_argument("--font-family", default="Arial", help="Preferred PDF/figure font.")
    parser.add_argument("--font-file", default=None, help="Optional absolute .ttf/.otf path; overrides fontconfig lookup.")
    args = parser.parse_args()

    # Register the font in Matplotlib itself, then embed TrueType in PDF.
    configure_editable_font(args.font_family, args.font_file)

    cytokines, cell_types, states = values(args.cytokine), values(args.cell_type), values(args.state)
    data = read_csv(args.data)
    required = {"patient_id", "cell_type", "cytokine", "state", "state_rank_score", args.cytokine_score_col}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Input table is missing columns: {sorted(missing)}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    selected = data.loc[
        data["cytokine"].astype(str).isin(cytokines)
        & data["cell_type"].astype(str).isin(cell_types)
        & data["state"].astype(str).isin(states)
    ].copy()
    if selected.empty:
        raise ValueError(
            "No rows matched. Check exact names with: "
            "cut -d, -f... or open state_cytokine_merged_data.csv."
        )

    selected.to_csv(outdir / "selected_state_cytokine_scatter_data.csv", index=False)
    records: list[dict] = []
    keys = ["cell_type", "state", "cytokine"]
    for key, table in selected.groupby(keys, sort=True):
        table = table.dropna(subset=["state_rank_score", args.cytokine_score_col]).copy()
        if len(table) < 3 or table["state_rank_score"].nunique() < 2 or table[args.cytokine_score_col].nunique() < 2:
            rho, p_value = np.nan, np.nan
        else:
            rho, p_value = spearmanr(table["state_rank_score"], table[args.cytokine_score_col])
        records.append(dict(cell_type=key[0], state=key[1], cytokine=key[2],
                            n_patients=len(table), spearman_rho=rho, p_value=p_value))
    statistics = pd.DataFrame(records)
    statistics.to_csv(outdir / "selected_state_cytokine_scatter_statistics.csv", index=False)

    combinations = list(selected.groupby(keys, sort=True))
    n_panels = len(combinations)
    n_cols = min(3, n_panels)
    n_rows = int(np.ceil(n_panels / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.3 * n_cols, 4.8 * n_rows), squeeze=False)
    for ax, (key, table) in zip(axes.ravel(), combinations):
        table = table.dropna(subset=["state_rank_score", args.cytokine_score_col]).copy()
        colour_column = args.color_by if args.color_by and args.color_by in table.columns else None
        if colour_column:
            levels = list(pd.unique(table[colour_column].astype(str)))
            for level in levels:
                subset = table.loc[table[colour_column].astype(str).eq(level)]
                ax.scatter(subset["state_rank_score"], subset[args.cytokine_score_col],
                           s=args.point_size, color=GROUP_COLORS.get(level, "#666666"),
                           edgecolor="white", linewidth=0.7, alpha=0.9, label=level)
            ax.legend(title=colour_column, frameon=False, fontsize=11, title_fontsize=11, loc="best")
        else:
            ax.scatter(table["state_rank_score"], table[args.cytokine_score_col],
                       s=args.point_size, color="#3B75AF", edgecolor="white", linewidth=0.7, alpha=0.9)
        x = table["state_rank_score"].to_numpy(dtype=float)
        y = table[args.cytokine_score_col].to_numpy(dtype=float)
        if len(table) >= 3 and np.unique(x).size > 1:
            slope, intercept = np.polyfit(x, y, 1)
            x_line = np.linspace(x.min(), x.max(), 100)
            ax.plot(x_line, slope * x_line + intercept, color="#222222", linewidth=1.4, zorder=0)
        row = statistics.loc[(statistics.cell_type == key[0]) & (statistics.state == key[1]) & (statistics.cytokine == key[2])].iloc[0]
        rho_text = "NA" if pd.isna(row.spearman_rho) else f"{row.spearman_rho:.2f}"
        ax.set_title(f"{key[0]} | {key[1]} | {key[2]}", fontsize=14, fontweight="bold", pad=10)
        # Retain only the effect-size annotation; P values and sample counts are in the CSV output.
        ax.text(0.03, 0.97, f"Spearman r = {rho_text}",
                transform=ax.transAxes, va="top", ha="left", fontsize=13, fontweight="bold")
        ax.set_xlabel("Cell-state rank score", fontsize=13, fontweight="bold")
        ax.set_ylabel(f"scCRS {key[2]} response score", fontsize=13, fontweight="bold")
        ax.tick_params(axis="both", labelsize=11, width=1.1)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(False)
    for ax in axes.ravel()[n_panels:]:
        ax.set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / "selected_state_cytokine_scatter_plots.pdf", bbox_inches="tight")
    fig.savefig(outdir / "selected_state_cytokine_scatter_plots.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {n_panels} panel(s) to {outdir}")


if __name__ == "__main__":
    main()
