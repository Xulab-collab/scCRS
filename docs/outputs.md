# Outputs

## Core scoring

- `cell_type_cytokine_scores.csv`: patient × annotated cell type × cytokine signed rank score, coverage and optional control calibration.
- `patient_cytokine_scores.csv`: unweighted and composition-weighted patient aggregate response scores.
- `patient_celltype_counts.csv`: patient × cell-type counts.
- `signature_gene_coverage.csv`: matched/available signature genes.

## Publication-ready visualisations

- Patient heatmap: row-z-score matrix, tidy plotting table, patient-group annotation, and vector PDF.
- Group-by-cell-type response comparison: multi-page vector PDF, exact input plotting data and BH-adjusted pairwise Wilcoxon table.

## Associations

- Composition: long merged data, all Spearman tests, BH q values, rho/q matrices and a PDF heatmap.
- State: per-cell state scores, merged cytokine/state data, correlations and selected scatterplots.

## Network atlas

- Patient × B-subtype response score table; complete gene–response associations; retained TF→target edges; group-difference edges; one editable network PDF per B-cell subtype (all cytokines and groups) when requested.

## Trajectory

- Per-cell pseudotime, branch and UCell score table; patient × bin table; GAM curves and statistics; PCRS branch summary; selected response-enhancing branch table; principal-graph curve coordinates; terminal diagnostics; UMAP PDF/SVG files.

All CSV tables should be retained with the figures to make filtering and thresholds auditable.
