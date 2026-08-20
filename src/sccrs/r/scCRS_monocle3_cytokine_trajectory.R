#!/usr/bin/env Rscript
# Post-Monocle3 scCRS cytokine-response trajectory analysis.
# Uses an existing, already ordered monocle3 cell_data_set. It does NOT rerun
# preprocess_cds(), learn_graph(), or reduce_dimension(). A Seurat object is
# optional and is used only as the source of the original RNA counts/UMAP.

suppressPackageStartupMessages({
  library(optparse); library(monocle3); library(Seurat); library(SingleCellExperiment)
  library(UCell); library(igraph); library(mgcv); library(dplyr); library(tidyr); library(ggplot2)
})

opts <- list(
  make_option("--cds", type="character", help="Existing monocle3 cell_data_set .rds or .RData"),
  make_option("--seurat", type="character", default=NULL, help="Optional original Seurat .rds/.RData; supplies RNA counts and original UMAP"),
  make_option("--assay", type="character", default="RNA"),
  make_option("--signatures", type="character", help="scCRS cytokine dictionary CSV"),
  make_option("--outdir", type="character", default="scCRS_monocle3_cytokine_trajectory"),
  make_option("--celltype-col", dest="celltype_col", type="character", default="cell_type"),
  make_option("--patient-col", dest="patient_col", type="character", default="patient_id"),
  make_option("--root-label", dest="root_label", type="character", help="Exact Naive-B root label; used to anchor both paths"),
  make_option("--branch1-label", dest="branch1_label", type="character", help="Exact terminal label of branch 1"),
  make_option("--branch2-label", dest="branch2_label", type="character", help="Exact terminal label of branch 2"),
  make_option("--branch1-name", dest="branch1_name", type="character", default="Branch 1", help="Display name for branch 1 in output tables and figures"),
  make_option("--branch2-name", dest="branch2_name", type="character", default="Branch 2", help="Display name for branch 2 in output tables and figures"),
  make_option("--branch1-path-labels", dest="branch1_path_labels", type="character", default=NULL, help="Optional comma-separated cell types retained for branch 1; root and terminal are always retained"),
  make_option("--branch2-path-labels", dest="branch2_path_labels", type="character", default=NULL, help="Optional comma-separated cell types retained for branch 2; root and terminal are always retained"),  make_option("--signature-celltype", dest="signature_celltype", type="character", default="B_cell"),
  make_option("--cytokines", type="character", default=NULL),
  make_option("--vector-cytokines", dest="vector_cytokines", type="character", default=NULL),
  make_option("--umap-cytokines", dest="umap_cytokines", type="character", default="vector", help="Comma-separated cytokines to render on UMAP, 'vector' (default; use --vector-cytokines), or 'all'."),
  make_option("--show-local-response-vectors", dest="show_local_response_vectors", action="store_true", default=FALSE, help="Also show short local response-slope arrows. Default: FALSE; use long lineage curves only."),
  make_option("--trajectory-line-points", dest="trajectory_line_points", type="integer", default=60, help="Number of points used to smooth each long UMAP lineage curve"),
  make_option("--min-cells-per-patient-bin", dest="min_cells_per_patient_bin", type="integer", default=10),
  make_option("--n-pseudotime-bins", dest="n_pseudotime_bins", type="integer", default=10),
  make_option("--min-patients", dest="min_patients", type="integer", default=5),
  make_option("--max-rank", dest="max_rank", type="integer", default=1500),
  make_option("--umap-point-size", dest="umap_point_size", type="double", default=0.50, help="Point size for UMAP response maps"),
  make_option("--umap-alpha", dest="umap_alpha", type="double", default=0.80, help="Point opacity for UMAP response maps"),
  make_option("--color-q-low", dest="color_q_low", type="double", default=0.00, help="Lower response-score quantile for UMAP color limits"),
  make_option("--color-q-high", dest="color_q_high", type="double", default=0.995, help="Upper response-score quantile for UMAP color limits; values above are saturated"),
  make_option("--umap-color-transform", dest="umap_color_transform", type="character", default="sqrt", help="Color transform: sqrt, identity, or log10"),
  make_option("--make-response-field", dest="make_response_field", action="store_true", default=FALSE, help="Export UMAP kernel-smoothed response-potential fields. Arrows indicate local score-gradient direction, not a causal perturbation force."),
  make_option("--field-grid", dest="field_grid", type="integer", default=50, help="Grid resolution per UMAP axis for response-potential fields"),
  make_option("--field-bandwidth", dest="field_bandwidth", type="double", default=NULL, help="Gaussian smoothing bandwidth in UMAP units; default is 8% of the mean UMAP span"),
  make_option("--field-max-cells", dest="field_max_cells", type="integer", default=6000, help="Maximum cells used to estimate each response field (random subsample for speed)"),
  make_option("--field-min-density", dest="field_min_density", type="double", default=0.08, help="Keep arrows only where relative kernel density is at least this value"),
  make_option("--field-arrow-stride", dest="field_arrow_stride", type="integer", default=4, help="Plot every Nth eligible grid location as a response-field arrow"),
  make_option("--make-density-vector-map", dest="make_density_vector_map", action="store_true", default=FALSE, help="Export UMAP kernel-smoothed mean-response heatmaps with local response-gradient arrows. The historical option name is retained for compatibility; this remains descriptive, not a causal perturbation prediction."),
  make_option("--heatmap-min-density", dest="heatmap_min_density", type="double", default=0.05, help="Minimum relative cell density required to display a mean-response heatmap pixel; low-support UMAP regions are masked to prevent boundary extrapolation artifacts."),
  make_option("--heatmap-bandwidth-factor", dest="heatmap_bandwidth_factor", type="double", default=0.05, help="When --field-bandwidth is not supplied, heatmap Gaussian bandwidth as a fraction of mean UMAP span. Smaller values yield tighter contours."),
  make_option("--heatmap-highlight-quantile", dest="heatmap_highlight_quantile", type="double", default=0.95, help="Overlay cells at or above this within-map raw UCell-score quantile so small high-response populations remain visible.")
)
opt <- parse_args(OptionParser(option_list=opts))
required <- c("cds", "signatures", "root_label", "branch1_label", "branch2_label")
if (any(vapply(required, function(x) is.null(opt[[x]]) || !nzchar(opt[[x]]), logical(1)))) stop("Required: --cds --signatures --root-label --branch1-label --branch2-label")
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

read_r_object <- function(path, wanted) {
  if (grepl("\\.RData$|\\.rda$", path, ignore.case=TRUE)) {
    e <- new.env(parent=emptyenv()); load(path, envir=e); xs <- mget(ls(e), envir=e)
    hits <- Filter(function(x) inherits(x, wanted), xs)
    if (!length(hits)) stop("No ", wanted, " object found in ", path)
    return(hits[[1]])
  }
  x <- readRDS(path)
  if (!inherits(x, wanted)) stop(path, " is not a ", wanted, " object")
  x
}
read_csv_encoding <- function(path) {
  for (enc in c("UTF-8", "UTF-8-BOM", "GB18030", "GBK")) {
    x <- try(read.csv(path, check.names=FALSE, fileEncoding=enc, stringsAsFactors=FALSE), silent=TRUE)
    if (!inherits(x, "try-error")) return(x)
  }
  stop("Cannot read signature table")
}
mode_value <- function(x) names(sort(table(x), decreasing=TRUE))[1]
closest_vertex <- function(cds, cell_ids) {
  graph <- principal_graph(cds)[["UMAP"]]
  vertex_names <- igraph::V(graph)$name
  raw <- as.character(principal_graph_aux(cds)[["UMAP"]]$pr_graph_cell_proj_closest_vertex[cell_ids, 1])
  # Different Monocle3 versions store closest vertices as Y_12, 12, or a
  # vertex index. Match each value independently; never add Y_ twice.
  out <- raw
  unresolved <- !out %in% vertex_names
  out[unresolved] <- paste0("Y_", raw[unresolved])
  unresolved <- !out %in% vertex_names
  numeric_index <- suppressWarnings(as.integer(raw[unresolved]))
  valid_index <- !is.na(numeric_index) & numeric_index >= 1 & numeric_index <= length(vertex_names)
  if (any(valid_index)) out[which(unresolved)[valid_index]] <- vertex_names[numeric_index[valid_index]]
  unresolved <- !out %in% vertex_names
  if (any(unresolved)) {
    stop("Could not map closest principal-graph vertices. Examples from cds: ",
         paste(unique(raw[unresolved])[seq_len(min(5, sum(unresolved)))], collapse=", "),
         "; examples in graph: ", paste(head(vertex_names, 5), collapse=", "))
  }
  out
}
safe_name <- function(x) {
  # Preserve common cytokine Greek letters in filenames. The former ASCII-only
  # sanitizer mapped beta, gamma, epsilon and kappa all to "_", overwriting
  # their individual PDF/SVG exports (for example, IFN-beta and IFN-gamma).
  from <- c("\u03b1", "\u03b2", "\u03b3", "\u03b4", "\u03b5", "\u03ba", "\u03bb")
  to <- c("alpha", "beta", "gamma", "delta", "epsilon", "kappa", "lambda")
  for (i in seq_along(from)) x <- gsub(from[i], to[i], x, fixed=TRUE)
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

cds <- read_r_object(opt$cds, "cell_data_set")
if (is.null(principal_graph(cds)[["UMAP"]]) || is.null(principal_graph_aux(cds)[["UMAP"]]$pr_graph_cell_proj_closest_vertex)) {
  stop("The supplied cds has no learned UMAP principal graph. Run learn_graph() first.")
}
ptime <- pseudotime(cds)
if (!any(is.finite(ptime))) stop("The supplied cds has no finite pseudotime. Run order_cells() with a Naive-B root first.")
for (col in c(opt$celltype_col, opt$patient_col)) if (!col %in% colnames(colData(cds))) stop("cds colData lacks: ", col)
cell_ids <- colnames(cds)
celltype <- as.character(colData(cds)[[opt$celltype_col]])
patient <- as.character(colData(cds)[[opt$patient_col]])
for (label in c(opt$root_label, opt$branch1_label, opt$branch2_label)) if (!any(celltype == label)) stop("Exact cell-type label absent from cds: ", label)

# Optional Seurat object provides original counts and the original UMAP layout.
if (!is.null(opt$seurat)) {
  seu <- read_r_object(opt$seurat, "Seurat")
  if (!opt$assay %in% Assays(seu)) stop("Seurat assay absent: ", opt$assay)
  missing_cells <- setdiff(cell_ids, colnames(seu))
  if (length(missing_cells)) stop(length(missing_cells), " cells in cds are absent from Seurat; ensure cell barcodes match.")
  X <- GetAssayData(seu, assay=opt$assay, slot="counts")[, cell_ids, drop=FALSE]
  if ("umap" %in% Reductions(seu)) {
    U <- Embeddings(seu, reduction="umap")[cell_ids, 1:2, drop=FALSE]
  } else U <- reducedDims(cds)[["UMAP"]][cell_ids, 1:2, drop=FALSE]
} else {
  X <- SingleCellExperiment::counts(cds)[, cell_ids, drop=FALSE]
  U <- reducedDims(cds)[["UMAP"]][cell_ids, 1:2, drop=FALSE]
}
if (is.null(X) || !nrow(X)) stop("No count matrix found. Supply --seurat with the raw RNA assay.")

# Derive two root-to-terminal paths on the EXISTING Monocle3 graph.
g <- principal_graph(cds)[["UMAP"]]
root_cells <- cell_ids[celltype == opt$root_label]
end1_cells <- cell_ids[celltype == opt$branch1_label]; end2_cells <- cell_ids[celltype == opt$branch2_label]
root_node <- mode_value(closest_vertex(cds, root_cells)); end1_node <- mode_value(closest_vertex(cds, end1_cells)); end2_node <- mode_value(closest_vertex(cds, end2_cells))
path_nodes <- function(a, b) {
  vertex_path <- shortest_paths(g, from=a, to=b, output="vpath")$vpath[[1]]
  # as_ids() retrieves the vertex 'name' attribute, unlike as.character()
  # which may return internal numerical IDs in some igraph releases.
  as.character(igraph::as_ids(vertex_path))
}
path1 <- path_nodes(root_node, end1_node); path2 <- path_nodes(root_node, end2_node)
if (!length(path1) || !length(path2)) stop("A specified terminal label is not connected to the root in the existing principal graph.")
closest <- closest_vertex(cds, cell_ids)
valid_vertices <- igraph::V(g)$name
for (item in list(closest=closest, path1=path1, path2=path2)) {
  bad <- unique(item[!item %in% valid_vertices])
  if (length(bad)) stop("Invalid principal-graph vertex mapping: ", paste(head(bad, 5), collapse=", "),
                        ". Available graph vertices include: ", paste(head(valid_vertices, 5), collapse=", "))
}
d1 <- apply(distances(g, v=closest, to=path1), 1, min); d2 <- apply(distances(g, v=closest, to=path2), 1, min)
branch <- ifelse(d1 < d2, "branch1", "branch2")
branch[closest %in% intersect(path1, path2)] <- "shared_trunk"
branch[!is.finite(ptime)] <- NA_character_
split_labels <- function(x) trimws(unlist(strsplit(x, ",", fixed=TRUE)))
allowed_labels <- list(
  branch1=unique(c(opt$root_label, opt$branch1_label, if (is.null(opt$branch1_path_labels)) character() else split_labels(opt$branch1_path_labels))),
  branch2=unique(c(opt$root_label, opt$branch2_label, if (is.null(opt$branch2_path_labels)) character() else split_labels(opt$branch2_path_labels)))
)
branch_display <- c(branch1=opt$branch1_name, branch2=opt$branch2_name, shared_trunk="Shared root/trunk")

# Coordinates for the existing Monocle3 principal-graph nodes. The preferred
# source is dp_mst; if unavailable, node positions are estimated from cells
# projected to each node. This keeps the UMAP curve anchored to the graph path.
graph_node_coordinates <- function(cds, graph, U, closest_nodes) {
  nodes <- igraph::V(graph)$name
  # First use the actual displayed UMAP: cells projected to a graph node provide
  # its coordinate in the same space used by the plot (including a Seurat UMAP).
  projected <- lapply(nodes, function(v) {
    ii <- which(closest_nodes==v)
    if (!length(ii)) return(NULL)
    data.frame(node=v, x=median(U[ii,1],na.rm=TRUE), y=median(U[ii,2],na.rm=TRUE), source="cell_projection", stringsAsFactors=FALSE)
  })
  out <- bind_rows(projected)
  # If a graph node has no projected cells, fall back to Monocle3's dp_mst layout.
  missing <- setdiff(nodes, out$node)
  aux <- principal_graph_aux(cds)[["UMAP"]]; mst <- aux$dp_mst
  if (length(missing) && !is.null(mst) && is.matrix(mst) && nrow(mst) >= 2) {
    if (!is.null(colnames(mst)) && any(colnames(mst) %in% missing)) {
      keep <- colnames(mst) %in% missing
      out <- bind_rows(out, data.frame(node=colnames(mst)[keep], x=as.numeric(mst[1,keep]), y=as.numeric(mst[2,keep]), source="dp_mst_fallback", stringsAsFactors=FALSE))
    } else if (!is.null(rownames(mst)) && any(rownames(mst) %in% missing) && ncol(mst) >= 2) {
      keep <- rownames(mst) %in% missing
      out <- bind_rows(out, data.frame(node=rownames(mst)[keep], x=as.numeric(mst[keep,1]), y=as.numeric(mst[keep,2]), source="dp_mst_fallback", stringsAsFactors=FALSE))
    }
  }
  out[is.finite(out$x) & is.finite(out$y),,drop=FALSE]
}
node_xy <- graph_node_coordinates(cds, g, U, closest)

# Signed UCell score: UCell supports GENE+ and GENE- in one response signature.
sig <- read_csv_encoding(opt$signatures)
need <- c("cytokine", "cell_type", "gene", "direction")
if (!all(need %in% colnames(sig))) stop("Signature file must contain: ", paste(need, collapse=", "))
sig <- sig[tolower(sig$cell_type) == tolower(opt$signature_celltype), , drop=FALSE]
sig$gene <- toupper(trimws(sig$gene)); sig$direction <- tolower(trimws(sig$direction))
if (!is.null(opt$cytokines)) sig <- sig[sig$cytokine %in% trimws(unlist(strsplit(opt$cytokines, ","))), , drop=FALSE]
all_genes <- toupper(rownames(X))
coverage <- sig %>% group_by(cytokine, direction) %>% summarise(total_genes=n(), matched_genes=sum(gene %in% all_genes), .groups="drop")
write.csv(coverage, file.path(opt$outdir, "cytokine_signature_gene_coverage.csv"), row.names=FALSE)
features <- lapply(split(sig, sig$cytokine), function(z) {
  z <- z[z$gene %in% all_genes & z$direction %in% c("up", "down"), , drop=FALSE]
  c(paste0(z$gene[z$direction == "up"], "+"), paste0(z$gene[z$direction == "down"], "-"))
})
features <- features[lengths(features) > 0]; cytokines <- names(features)
if (!length(features)) stop("No signature genes matched the expression matrix.")
u_scores <- as.data.frame(UCell::ScoreSignatures_UCell(X, features=features, maxRank=opt$max_rank))
rownames(u_scores) <- cell_ids
score_col <- vapply(cytokines, function(x) {
  hit <- c(paste0(x, "_UCell"), x); hit <- hit[hit %in% colnames(u_scores)]
  if (!length(hit)) stop("UCell output does not contain: ", x); hit[1]
}, character(1))

meta <- data.frame(cell_id=cell_ids, patient_id=patient, cell_type=celltype, pseudotime=as.numeric(ptime), branch=branch, branch_name=unname(branch_display[branch]), UMAP_1=U[,1], UMAP_2=U[,2], stringsAsFactors=FALSE)
for (cy in cytokines) meta[[cy]] <- u_scores[cell_ids, score_col[[cy]]]
write.csv(meta, file.path(opt$outdir, "Bcell_monocle3_pseudotime_and_cytokine_scores.csv"), row.names=FALSE)
saveRDS(cds, file.path(opt$outdir, "Bcell_monocle3_cds_with_scCRS.rds"))

# Patient 闂?pseudotime-bin GAM analysis, separately for each branch.
fit_one <- function(cytokine, branch_name) {
  z <- meta[((meta$branch == branch_name & meta$cell_type %in% allowed_labels[[branch_name]]) | meta$cell_type == opt$root_label) & is.finite(meta$pseudotime), , drop=FALSE]
  if (length(unique(z$patient_id)) < opt$min_patients) return(NULL)
  br <- unique(quantile(z$pseudotime, probs=seq(0,1,length.out=opt$n_pseudotime_bins+1), na.rm=TRUE)); if (length(br) < 4) return(NULL)
  z$bin <- cut(z$pseudotime, breaks=br, include.lowest=TRUE)
  b <- z %>% group_by(patient_id, bin) %>% summarise(pseudotime=median(pseudotime), score=mean(.data[[cytokine]]), n_cells=n(), .groups="drop") %>% filter(n_cells >= opt$min_cells_per_patient_bin)
  # Numerical QC before GAM: finite, nonconstant response and at least two
  # pseudotime bins per patient are required for a patient random effect.
  b <- b %>% filter(is.finite(pseudotime), is.finite(score)) %>% group_by(patient_id) %>% filter(n() >= 2) %>% ungroup()
  if (nrow(b) < 12 || length(unique(b$pseudotime)) < 4 || length(unique(b$patient_id)) < opt$min_patients || sd(b$score) < 1e-8) return(NULL)
  b$patient_id <- droplevels(factor(as.character(b$patient_id)))
  k <- min(5, length(unique(b$pseudotime))-1, floor(nrow(b)/3))
  if (k < 3) return(NULL)
  # Primary model: smooth pseudotime plus a patient random intercept.
  fit <- tryCatch(gam(score ~ s(pseudotime, k=k) + s(patient_id, bs="re"), data=b, method="REML"), error=function(e) NULL)
  model_type <- "patient_random_effect"
  # Some older mgcv builds fail to optimize an RE penalty for sparse bins.
  # Preserve patient adjustment with a fixed-effect fallback and label it.
  if (is.null(fit)) {
    fit <- tryCatch(gam(score ~ s(pseudotime, k=k) + patient_id, data=b, method="REML"), error=function(e) NULL)
    model_type <- "patient_fixed_effect_fallback"
  }
  if (is.null(fit)) return(NULL)
  gr <- data.frame(pseudotime=seq(quantile(z$pseudotime,.05), quantile(z$pseudotime,.95), length.out=80), patient_id=factor(levels(b$patient_id)[1], levels=levels(b$patient_id)))
  pr <- if (model_type == "patient_random_effect") predict(fit, newdata=gr, se.fit=TRUE, exclude="s(patient_id)") else predict(fit, newdata=gr, se.fit=TRUE)
  gr$fit <- as.numeric(pr$fit); gr$se <- as.numeric(pr$se.fit); gr$lower <- gr$fit-1.96*gr$se; gr$upper <- gr$fit+1.96*gr$se
  gr$slope <- c(diff(gr$fit)/diff(gr$pseudotime), NA_real_)
  st <- data.frame(cytokine=cytokine, branch=branch_name, model_type=model_type, n_patients=length(unique(b$patient_id)), n_patient_bins=nrow(b), spearman_rho=cor(b$pseudotime,b$score,method="spearman"), GAM_p_value=summary(fit)$s.table[1,"p-value"], PCRS=(tail(gr$fit,1)-gr$fit[1])/sd(b$score))
  list(stats=st, bins=transform(b, cytokine=cytokine, branch=branch_name), curve=transform(gr, cytokine=cytokine, branch=branch_name))
}
res <- list(); k <- 1
for (cy in cytokines) for (b in c("branch1","branch2")) { x <- fit_one(cy,b); if (!is.null(x)) {res[[k]] <- x; k <- k+1} }
if (!length(res)) stop("No branch model could be fitted; inspect patient counts and --min-cells-per-patient-bin.")
stats <- bind_rows(lapply(res, `[[`, "stats")); stats$branch_name <- unname(branch_display[stats$branch]); stats$BH_FDR <- p.adjust(stats$GAM_p_value, method="BH")
bins <- bind_rows(lapply(res, `[[`, "bins")); curves <- bind_rows(lapply(res, `[[`, "curve")); bins$branch_name <- unname(branch_display[bins$branch]); curves$branch_name <- unname(branch_display[curves$branch])
write.csv(stats, file.path(opt$outdir,"cytokine_branch_pseudotime_association_statistics.csv"),row.names=FALSE)
write.csv(bins,file.path(opt$outdir,"patient_pseudotime_bin_cytokine_scores.csv"),row.names=FALSE)
write.csv(curves,file.path(opt$outdir,"cytokine_branch_GAM_curves.csv"),row.names=FALSE)
branch_summary <- stats %>% select(cytokine,branch,PCRS,BH_FDR) %>% pivot_wider(names_from=branch,values_from=c(PCRS,BH_FDR)) %>% mutate(PCRS_branch_difference=PCRS_branch1-PCRS_branch2)
write.csv(branch_summary,file.path(opt$outdir,"cytokine_branch_specificity_summary.csv"),row.names=FALSE)

# UMAP overlays. Long curves show the user-defined, existing Monocle3 lineage
# direction. Optional short arrows show the local GAM response slope.
make_arrows <- function(cytokine, b) {
  cv <- curves[curves$cytokine==cytokine & curves$branch==b,]; z <- meta[((meta$branch==b & meta$cell_type %in% allowed_labels[[b]]) | meta$cell_type==opt$root_label) & is.finite(meta$pseudotime),]
  if (!nrow(cv) || nrow(z)<30) return(data.frame())
  a <- cv[round(seq(4,nrow(cv)-4,length.out=10)),]; radius <- diff(range(z$pseudotime))/16
  a$x <- vapply(a$pseudotime,function(t) median(z$UMAP_1[abs(z$pseudotime-t)<=radius],na.rm=TRUE),numeric(1)); a$y <- vapply(a$pseudotime,function(t) median(z$UMAP_2[abs(z$pseudotime-t)<=radius],na.rm=TRUE),numeric(1)); a <- a[complete.cases(a),]
  if(nrow(a)<3) return(data.frame()); dx <- c(diff(a$x),NA); dy <- c(diff(a$y),NA); nm <- sqrt(dx^2+dy^2); mag <- abs(a$slope)/max(abs(a$slope),na.rm=TRUE); step <- .09*max(diff(range(meta$UMAP_1)),diff(range(meta$UMAP_2)))*mag; sg <- sign(a$slope); sg[!is.finite(sg)] <- 0
  data.frame(cytokine=cytokine,branch=b,x=a$x,y=a$y,xend=a$x+sg*step*dx/nm,yend=a$y+sg*step*dy/nm,slope=a$slope)
}

make_trajectory_curve <- function(b, n_points=60) {
  nodes <- if (b=="branch1") path1 else path2
  z <- node_xy[match(nodes,node_xy$node), c("node","x","y"), drop=FALSE]
  z <- z[complete.cases(z),,drop=FALSE]
  if (nrow(z) < 2) return(data.frame())
  # Parameterise by arc length instead of cell-density quantiles. Thus a highly
  # abundant Naive-B population cannot pull a lineage curve away from its terminal.
  d <- c(0, cumsum(sqrt(diff(z$x)^2 + diff(z$y)^2)))
  if (max(d) <= 0) return(data.frame())
  tt <- d/max(d)
  keep <- c(TRUE, diff(tt) > 1e-8)
  z <- z[keep,,drop=FALSE]; tt <- tt[keep]
  pred_t <- seq(0,1,length.out=max(20,n_points))
  if (nrow(z) >= 4) {
    sx <- tryCatch(smooth.spline(tt,z$x,spar=.55), error=function(e) NULL)
    sy <- tryCatch(smooth.spline(tt,z$y,spar=.55), error=function(e) NULL)
    if (!is.null(sx) && !is.null(sy)) {
      xx <- predict(sx,pred_t)$y; yy <- predict(sy,pred_t)$y
    } else { xx <- approx(tt,z$x,xout=pred_t,rule=2)$y; yy <- approx(tt,z$y,xout=pred_t,rule=2)$y }
  } else { xx <- approx(tt,z$x,xout=pred_t,rule=2)$y; yy <- approx(tt,z$y,xout=pred_t,rule=2)$y }
  # Preserve exact root and terminal graph-node coordinates despite smoothing.
  xx[1] <- z$x[1]; yy[1] <- z$y[1]; xx[length(xx)] <- z$x[nrow(z)]; yy[length(yy)] <- z$y[nrow(z)]
  data.frame(branch=b, order=seq_along(pred_t), node_start=z$node[1], node_end=z$node[nrow(z)], x=xx, y=yy)
}

# Keep the user-specified vector list intact and audit it explicitly. A cytokine
# is only unavailable when its requested B-cell signature has no matched gene;
# a missing GAM fit must never suppress its UMAP / response-field export.
if (is.null(opt$vector_cytokines)) {
  vector_cytokines <- stats %>% group_by(cytokine) %>% summarise(x=max(abs(PCRS)), .groups="drop") %>% arrange(desc(x)) %>% slice_head(n=6) %>% pull(cytokine)
} else {
  vector_cytokines <- trimws(unlist(strsplit(opt$vector_cytokines, ",", fixed=TRUE)))
}
vector_cytokines <- unique(vector_cytokines[nzchar(vector_cytokines)])
if (!length(vector_cytokines)) stop("No cytokines were supplied through --vector-cytokines.")
missing_vector_signature <- setdiff(vector_cytokines, cytokines)
if (length(missing_vector_signature)) warning("These --vector-cytokines have no usable matched signature and cannot be plotted: ", paste(missing_vector_signature, collapse=", "))
available_vector_cytokines <- intersect(vector_cytokines, cytokines)
# By default every available --vector-cytokines member is exported. Do not use
# the GAM/PCRS table as a filter: maps are still valid descriptive score maps
# when a branch model cannot be fitted.
if (tolower(trimws(opt$umap_cytokines)) %in% c("vector", "selected")) {
  umap_cytokines <- available_vector_cytokines
} else if (tolower(trimws(opt$umap_cytokines)) %in% c("all", "*")) {
  umap_cytokines <- cytokines
} else {
  umap_cytokines <- unique(trimws(unlist(strsplit(opt$umap_cytokines, ",", fixed=TRUE))))
  missing_umap <- setdiff(umap_cytokines, cytokines)
  if (length(missing_umap)) warning("Requested UMAP cytokines not available after signature filtering: ", paste(missing_umap, collapse=", "))
  umap_cytokines <- intersect(umap_cytokines, cytokines)
}
if (!length(umap_cytokines)) stop("No requested cytokines have usable matched signatures for UMAP export.")
export_manifest <- data.frame(
  cytokine=vector_cytokines,
  requested_in_vector_cytokines=TRUE,
  signature_score_available=vector_cytokines %in% cytokines,
  included_in_default_umap_export=vector_cytokines %in% umap_cytokines,
  GAM_path_stat_available=vector_cytokines %in% stats$cytokine,
  stringsAsFactors=FALSE
)
write.csv(export_manifest, file.path(opt$outdir, "UMAP_cytokine_export_manifest.csv"), row.names=FALSE)
# One comparison figure for the selected cytokines: both branch GAM curves in each panel.
curve_plot <- curves %>% filter(cytokine %in% vector_cytokines) %>% group_by(cytokine, branch) %>% mutate(normalized_pseudotime=(pseudotime-min(pseudotime))/(max(pseudotime)-min(pseudotime))) %>% ungroup()
if (nrow(curve_plot)) {
  p_curve <- ggplot(curve_plot, aes(normalized_pseudotime, fit, color=branch, fill=branch)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=.16, linewidth=0, color=NA) + geom_line(linewidth=.85) +
    facet_wrap(~cytokine, scales="free_y", ncol=3) +
    scale_color_manual(values=c(branch1="#D55E00", branch2="#0072B2"), labels=c(branch1=opt$branch1_name, branch2=opt$branch2_name), name="Terminal path") +
    scale_fill_manual(values=c(branch1="#D55E00", branch2="#0072B2"), labels=c(branch1=opt$branch1_name, branch2=opt$branch2_name), name="Terminal path") +
    labs(x="Normalized Monocle3 pseudotime", y="UCell cytokine-response score", title="Branch-specific cytokine-response dynamics") +
    theme_classic(base_size=12) + theme(legend.position="top")
  ggsave(file.path(opt$outdir,"branch_cytokine_response_comparison.pdf"), p_curve, width=10, height=max(5.5,3*ceiling(length(unique(curve_plot$cytokine))/3)), device=cairo_pdf)
  ggsave(file.path(opt$outdir,"branch_cytokine_response_comparison.svg"), p_curve, width=10, height=max(5.5,3*ceiling(length(unique(curve_plot$cytokine))/3)))
}

vdir <- file.path(opt$outdir,"UMAP_response_aligned_vectors"); dir.create(vdir,showWarnings=FALSE); all_arrows <- list()
long_curves <- bind_rows(make_trajectory_curve("branch1",opt$trajectory_line_points), make_trajectory_curve("branch2",opt$trajectory_line_points))
write.csv(long_curves,file.path(opt$outdir,"UMAP_principal_graph_lineage_curve_coordinates.csv"),row.names=FALSE)
# Diagnose whether the selected terminal graph node lies near the intended cell-type island.
terminal_diag <- bind_rows(
  data.frame(branch="branch1", terminal_label=opt$branch1_label, terminal_node=end1_node, n_terminal_cells=sum(celltype==opt$branch1_label), n_finite_pseudotime=sum(celltype==opt$branch1_label & is.finite(ptime))),
  data.frame(branch="branch2", terminal_label=opt$branch2_label, terminal_node=end2_node, n_terminal_cells=sum(celltype==opt$branch2_label), n_finite_pseudotime=sum(celltype==opt$branch2_label & is.finite(ptime)))) %>%
  rowwise() %>% mutate(terminal_cell_x=median(U[celltype==terminal_label,1],na.rm=TRUE), terminal_cell_y=median(U[celltype==terminal_label,2],na.rm=TRUE), graph_node_x=node_xy$x[match(terminal_node,node_xy$node)], graph_node_y=node_xy$y[match(terminal_node,node_xy$node)], terminal_node_distance=sqrt((terminal_cell_x-graph_node_x)^2+(terminal_cell_y-graph_node_y)^2)) %>% ungroup()
write.csv(terminal_diag,file.path(opt$outdir,"UMAP_terminal_path_diagnostics.csv"),row.names=FALSE)
line_colors <- c(branch1="#D55E00",branch2="#0072B2")
# Select exactly one branch per cytokine: the path with the largest absolute PCRS.
# PCRS sign controls the arrow direction, so it points toward increasing response.
best_path <- stats %>% filter(is.finite(PCRS)) %>% group_by(cytokine) %>%
  slice_max(order_by=abs(PCRS), n=1, with_ties=FALSE) %>%
  transmute(cytokine, selected_branch=branch, PCRS, response_direction=ifelse(PCRS >= 0, "increasing_pseudotime", "decreasing_pseudotime")) %>% ungroup()
write.csv(best_path,file.path(opt$outdir,"UMAP_selected_response_enhancing_branch.csv"),row.names=FALSE)
all_pdf <- file.path(vdir,"selected_cytokines_UMAP_response_enhancing_trajectory.pdf")
grDevices::cairo_pdf(all_pdf,width=8.7,height=7.0)
for (cy in umap_cytokines) {
  ar <- bind_rows(make_arrows(cy,"branch1"),make_arrows(cy,"branch2")); all_arrows[[cy]] <- ar
  plot_cells <- meta[is.finite(meta$pseudotime),]
  color_limits <- as.numeric(quantile(plot_cells[[cy]], probs=c(opt$color_q_low, opt$color_q_high), na.rm=TRUE, names=FALSE))
  if (!all(is.finite(color_limits)) || color_limits[2] <= color_limits[1]) color_limits <- range(plot_cells[[cy]], na.rm=TRUE)
  p <- ggplot(plot_cells,aes(UMAP_1,UMAP_2,color=.data[[cy]]))+
    geom_point(size=opt$umap_point_size,alpha=opt$umap_alpha)+
    scale_color_viridis_c(name="UCell cytokine\nresponse",limits=color_limits,oob=scales::squish,trans=opt$umap_color_transform)
  selected <- best_path[best_path$cytokine==cy,,drop=FALSE]
  has_response_trajectory <- FALSE
  trajectory_legend_label <- "Monocle3 root-to-terminal\npseudotime branch (inferred)"
  if (nrow(selected)) {
    b <- selected$selected_branch[1]
    lc <- long_curves[long_curves$branch==b,,drop=FALSE]
    if (nrow(lc) >= 3) {
      p <- p + geom_path(data=lc,aes(x=x,y=y,linetype=trajectory_legend_label),inherit.aes=FALSE,color=line_colors[[b]],linewidth=1.35,alpha=.98,show.legend=TRUE)
      has_response_trajectory <- TRUE
      tip <- lc[nrow(lc),,drop=FALSE]; prev <- lc[nrow(lc)-1,,drop=FALSE]
      p <- p + geom_segment(data=data.frame(x=prev$x,y=prev$y,xend=tip$x,yend=tip$y),aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,color=line_colors[[b]],linewidth=1.35,arrow=arrow(length=grid::unit(3.4,"mm")))
    }
    display_branch <- if (b=="branch1") opt$branch1_name else opt$branch2_name
    direction_text <- if (selected$PCRS[1] >= 0) "increasing response with pseudotime" else "increasing response toward the root"
    subtitle_text <- paste0("Observed Monocle3 branch: ",display_branch," (largest |PCRS| across the two paths; ",direction_text,"). Arrow follows root-to-terminal pseudotime; response increase or decrease is reported separately by PCRS.")
  } else {
    subtitle_text <- "No eligible GAM/PCRS result for either path; no response-enhancing curve was drawn."
  }
  if (isTRUE(opt$show_local_response_vectors) && nrow(ar)) p <- p + geom_segment(data=ar,aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,color="black",linewidth=.45,arrow=arrow(length=grid::unit(1.8,"mm")))
  p <- p + coord_equal()+theme_classic(base_size=12)+labs(
    title=paste0(cy,": B-cell cytokine-response map"),
    subtitle=paste0(subtitle_text," Color: ",opt$umap_color_transform," transform; upper values saturated at q=",opt$color_q_high,"."))
  if (has_response_trajectory) {
    p <- p + scale_linetype_manual(name="Trajectory overlay", values=setNames("solid",trajectory_legend_label)) +
      guides(linetype=guide_legend(override.aes=list(color=line_colors[[b]],linewidth=1.35),order=1), color=guide_colorbar(order=2)) +
      theme(legend.position="bottom", legend.box="vertical")
  }
  ggsave(file.path(vdir,paste0(safe_name(cy),"_UMAP_response_aligned_vector.pdf")),p,width=8.5,height=6.8,device=cairo_pdf)
  ggsave(file.path(vdir,paste0(safe_name(cy),"_UMAP_response_aligned_vector.svg")),p,width=8.5,height=6.8)
  print(p)
}
grDevices::dev.off()
write.csv(bind_rows(all_arrows),file.path(opt$outdir,"UMAP_response_aligned_vector_coordinates.csv"),row.names=FALSE)

# Optional descriptive response-potential field. This is deliberately not a
# perturbation-prediction model: it estimates the local gradient of a Gaussian
# kernel-smoothed single-cell UCell score in the displayed UMAP coordinates.
# Therefore arrows identify where response-associated transcription increases
# in the embedding, not cellular motion, a physical force, or causal cytokine action.
make_response_field <- function(z, cytokine, n_grid=50, bandwidth=NULL, max_cells=6000,
                                min_density=.08, arrow_stride=4) {
  z <- z[is.finite(z$UMAP_1) & is.finite(z$UMAP_2) & is.finite(z[[cytokine]]),
         c("UMAP_1", "UMAP_2", cytokine), drop=FALSE]
  if (nrow(z) < 50) return(NULL)
  # Define grid limits from all displayed cells before sampling. Sampling only
  # speeds kernel estimation and must never contract the plotted UMAP extent.
  xr <- range(z$UMAP_1); yr <- range(z$UMAP_2)
  if (nrow(z) > max_cells) {
    set.seed(20260818)
    z <- z[sample.int(nrow(z), max_cells), , drop=FALSE]
  }
  span_x <- diff(xr); span_y <- diff(yr)
  if (!is.finite(span_x) || !is.finite(span_y) || span_x <= 0 || span_y <= 0) return(NULL)
  n_grid <- max(20, as.integer(n_grid))
  xg <- seq(xr[1], xr[2], length.out=n_grid)
  yg <- seq(yr[1], yr[2], length.out=n_grid)
  bw <- if (is.null(bandwidth) || !is.finite(bandwidth) || bandwidth <= 0) {
    .08 * mean(c(span_x, span_y))
  } else bandwidth
  grid_df <- expand.grid(UMAP_1=xg, UMAP_2=yg, KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  grid_df$ix <- match(grid_df$UMAP_1, xg); grid_df$iy <- match(grid_df$UMAP_2, yg)
  score <- z[[cytokine]]
  # Rank weighting makes the density panel show where high-response cells are
  # concentrated, while avoiding domination by the absolute UCell score scale.
  score_rank <- rank(score, ties.method="average") / (length(score) + 1)
  estimate_one <- function(x, y) {
    w <- exp(-((z$UMAP_1-x)^2 + (z$UMAP_2-y)^2) / (2*bw^2))
    sw <- sum(w)
    c(
      potential=if (sw > 0) sum(w*score)/sw else NA_real_,
      density=mean(w),
      response_density=mean(w*score_rank)
    )
  }
  est <- t(vapply(seq_len(nrow(grid_df)), function(i) estimate_one(grid_df$UMAP_1[i], grid_df$UMAP_2[i]), numeric(3)))
  grid_df$potential <- est[, "potential"]; grid_df$density <- est[, "density"]; grid_df$response_density <- est[, "response_density"]
  grid_df$relative_density <- grid_df$density / max(grid_df$density, na.rm=TRUE)
  grid_df$relative_response_density <- grid_df$response_density / max(grid_df$response_density, na.rm=TRUE)
  pmat <- matrix(grid_df$potential, nrow=length(xg), ncol=length(yg))
  gx <- matrix(NA_real_, nrow=length(xg), ncol=length(yg)); gy <- gx
  if (length(xg) > 2) {
    gx[2:(length(xg)-1), ] <- (pmat[3:length(xg), ] - pmat[1:(length(xg)-2), ]) / (xg[3] - xg[1])
    gx[1, ] <- (pmat[2, ]-pmat[1, ]) / (xg[2]-xg[1]); gx[length(xg), ] <- (pmat[length(xg), ]-pmat[length(xg)-1, ]) / (xg[length(xg)]-xg[length(xg)-1])
  }
  if (length(yg) > 2) {
    gy[, 2:(length(yg)-1)] <- (pmat[, 3:length(yg)] - pmat[, 1:(length(yg)-2)]) / (yg[3] - yg[1])
    gy[, 1] <- (pmat[, 2]-pmat[, 1]) / (yg[2]-yg[1]); gy[, length(yg)] <- (pmat[, length(yg)]-pmat[, length(yg)-1]) / (yg[length(yg)]-yg[length(yg)-1])
  }
  grid_df$grad_x <- as.vector(gx); grid_df$grad_y <- as.vector(gy)
  grid_df$magnitude <- sqrt(grid_df$grad_x^2 + grid_df$grad_y^2)
  keep <- is.finite(grid_df$magnitude) & is.finite(grid_df$relative_density) & grid_df$relative_density >= min_density
  qmag <- if (any(keep)) as.numeric(quantile(grid_df$magnitude[keep], .90, na.rm=TRUE)) else NA_real_
  if (!is.finite(qmag) || qmag <= 0) qmag <- 1
  arrow_length <- .05 * max(span_x, span_y)
  grid_df$field_dx <- ifelse(keep, arrow_length * pmin(grid_df$magnitude/qmag, 1) * grid_df$grad_x/pmax(grid_df$magnitude, 1e-12), NA_real_)
  grid_df$field_dy <- ifelse(keep, arrow_length * pmin(grid_df$magnitude/qmag, 1) * grid_df$grad_y/pmax(grid_df$magnitude, 1e-12), NA_real_)
  grid_df$xend <- grid_df$UMAP_1 + grid_df$field_dx; grid_df$yend <- grid_df$UMAP_2 + grid_df$field_dy
  arrows <- grid_df[keep & (grid_df$ix %% max(1, arrow_stride) == 1) & (grid_df$iy %% max(1, arrow_stride) == 1), , drop=FALSE]
  list(grid=grid_df, arrows=arrows, bandwidth=bw, n_cells=nrow(z))
}

if (isTRUE(opt$make_response_field)) {
  fdir <- file.path(opt$outdir, "UMAP_response_potential_fields"); dir.create(fdir, showWarnings=FALSE)
  field_grids <- list()
  all_field_pdf <- file.path(fdir, "selected_cytokines_UMAP_response_potential_fields.pdf")
  grDevices::cairo_pdf(all_field_pdf, width=8.7, height=7.0)
  for (cy in umap_cytokines) {
    plot_cells <- meta[is.finite(meta$pseudotime), , drop=FALSE]
    fld <- make_response_field(plot_cells, cy, n_grid=opt$field_grid, bandwidth=opt$field_bandwidth,
                               max_cells=opt$field_max_cells, min_density=opt$field_min_density,
                               arrow_stride=opt$field_arrow_stride)
    if (is.null(fld)) { warning("Response field skipped for ", cy, ": too few valid cells."); next }
    fld$grid$cytokine <- cy; fld$grid$bandwidth <- fld$bandwidth; fld$grid$n_cells_used <- fld$n_cells
    field_grids[[cy]] <- fld$grid
    color_limits <- as.numeric(quantile(plot_cells[[cy]], probs=c(opt$color_q_low, opt$color_q_high), na.rm=TRUE, names=FALSE))
    if (!all(is.finite(color_limits)) || color_limits[2] <= color_limits[1]) color_limits <- range(plot_cells[[cy]], na.rm=TRUE)
    p <- ggplot(plot_cells, aes(UMAP_1, UMAP_2, color=.data[[cy]])) +
      geom_point(size=opt$umap_point_size, alpha=opt$umap_alpha) +
      scale_color_viridis_c(name="UCell cytokine\nresponse", limits=color_limits, oob=scales::squish, trans=opt$umap_color_transform) +
      geom_segment(data=fld$arrows, aes(x=UMAP_1, y=UMAP_2, xend=xend, yend=yend, linetype="Local response-potential gradient"),
                   inherit.aes=FALSE, color="#1A1A1A", linewidth=.42, alpha=.82, arrow=arrow(length=grid::unit(1.5, "mm")), show.legend=TRUE)
    selected <- best_path[best_path$cytokine == cy, , drop=FALSE]
    if (nrow(selected)) {
      b <- selected$selected_branch[1]; lc <- long_curves[long_curves$branch == b, , drop=FALSE]
      if (nrow(lc) >= 3) {
        p <- p + geom_path(data=lc, aes(x=x, y=y, linetype="Pseudotime branch with strongest cytokine-response association"),
                           inherit.aes=FALSE, color=line_colors[[b]], linewidth=1.35, show.legend=TRUE)
        tip <- lc[nrow(lc), , drop=FALSE]; prev <- lc[nrow(lc)-1, , drop=FALSE]
        p <- p + geom_segment(data=data.frame(x=prev$x, y=prev$y, xend=tip$x, yend=tip$y), aes(x=x, y=y, xend=xend, yend=yend),
                              inherit.aes=FALSE, color=line_colors[[b]], linewidth=1.35, arrow=arrow(length=grid::unit(3.4, "mm")))
      }
    }
    p <- p + coord_equal() + theme_classic(base_size=12) +
      scale_linetype_manual(name="UMAP overlay", values=c("Local response-potential gradient"="solid", "Pseudotime branch with strongest cytokine-response association"="solid")) +
      guides(linetype=guide_legend(override.aes=list(color=c("#1A1A1A", "#D55E00"), linewidth=c(.42, 1.35)), order=1), color=guide_colorbar(order=2)) +
      theme(legend.position="bottom", legend.box="vertical") +
      labs(title=paste0(cy, ": UMAP cytokine-response potential field"),
           subtitle=paste0("Arrows follow local increases in a Gaussian kernel-smoothed UCell response score (bandwidth=", signif(fld$bandwidth, 3), "). Descriptive only: not a predicted perturbation, physical force, or causal effect."))
    ggsave(file.path(fdir, paste0(safe_name(cy), "_UMAP_response_potential_field.pdf")), p, width=8.5, height=6.8, device=cairo_pdf)
    ggsave(file.path(fdir, paste0(safe_name(cy), "_UMAP_response_potential_field.svg")), p, width=8.5, height=6.8)
    print(p)
  }
  grDevices::dev.off()
  if (length(field_grids)) write.csv(bind_rows(field_grids), file.path(opt$outdir, "UMAP_response_potential_field_grid.csv"), row.names=FALSE)
}

# Squidiff Fig. 4c-inspired display: a continuous UMAP density surface for
# high-response cells, plus arrows from the local gradient of the smoothed
# response score. The density is observational and does not show generated or
# counterfactual cells; the label intentionally avoids causal language.
if (isTRUE(opt$make_density_vector_map)) {
  ddir <- file.path(opt$outdir, "UMAP_mean_response_heatmap_vector_maps"); dir.create(ddir, showWarnings=FALSE)
  mean_response_grids <- list(); mean_response_score_diagnostics <- list()
  all_density_pdf <- file.path(ddir, "selected_cytokines_UMAP_response_mean_heatmap_vector_maps.pdf")
  grDevices::cairo_pdf(all_density_pdf, width=8.7, height=7.0)
  for (cy in umap_cytokines) {
    plot_cells <- meta[is.finite(meta$pseudotime), , drop=FALSE]
    heatmap_bandwidth <- if (!is.null(opt$field_bandwidth) && is.finite(opt$field_bandwidth) && opt$field_bandwidth > 0) {
      opt$field_bandwidth
    } else {
      opt$heatmap_bandwidth_factor * mean(c(diff(range(plot_cells$UMAP_1, na.rm=TRUE)), diff(range(plot_cells$UMAP_2, na.rm=TRUE))))
    }
    fld <- make_response_field(plot_cells, cy, n_grid=opt$field_grid, bandwidth=heatmap_bandwidth,
                               max_cells=opt$field_max_cells, min_density=opt$field_min_density,
                               arrow_stride=opt$field_arrow_stride)
    if (is.null(fld)) { warning("Mean-response heatmap skipped for ", cy, ": too few valid cells."); next }
    fld$grid$cytokine <- cy; fld$grid$bandwidth <- fld$bandwidth; fld$grid$n_cells_used <- fld$n_cells
    mean_response_grids[[cy]] <- fld$grid
    # A kernel-smoothed mean is unstable where its denominator is almost zero.
    # Mask unsupported grid pixels instead of extrapolating into empty UMAP space.
    support_cutoff <- max(opt$field_min_density, opt$heatmap_min_density)
    fld$grid$potential_in_support <- ifelse(fld$grid$relative_density >= support_cutoff, fld$grid$potential, NA_real_)
    # Plot the kernel-normalized local mean response, not a response-weighted
    # cell count. Unsupported UMAP space is explicitly set to the lowest
    # supported score so it remains purple rather than being extrapolated.
    potential_floor <- as.numeric(quantile(fld$grid$potential_in_support, probs=opt$color_q_low, na.rm=TRUE, names=FALSE))
    if (!is.finite(potential_floor)) potential_floor <- min(fld$grid$potential, na.rm=TRUE)
    fld$grid$potential_for_plot <- ifelse(is.finite(fld$grid$potential_in_support), fld$grid$potential_in_support, potential_floor)
    mean_response_grids[[cy]] <- fld$grid
    potential_limits <- as.numeric(quantile(fld$grid$potential_in_support, probs=c(opt$color_q_low, opt$color_q_high), na.rm=TRUE, names=FALSE))
    if (!all(is.finite(potential_limits)) || potential_limits[2] <= potential_limits[1]) potential_limits <- range(fld$grid$potential_in_support, na.rm=TRUE)
    potential_transform <- if (potential_limits[1] < 0) "identity" else opt$umap_color_transform
    individual_limits <- as.numeric(quantile(plot_cells[[cy]], probs=c(opt$color_q_low, opt$color_q_high), na.rm=TRUE, names=FALSE))
    if (!all(is.finite(individual_limits)) || individual_limits[2] <= individual_limits[1]) individual_limits <- range(plot_cells[[cy]], na.rm=TRUE)
    individual_transform <- if (individual_limits[1] < 0) "identity" else opt$umap_color_transform
    # Use an exact top fraction rather than score >= quantile: many cytokine
    # signatures have zero-inflated UCell values, for which a quantile cutoff
    # of zero would otherwise recolour nearly every cell dark purple.
    finite_order <- order(-plot_cells[[cy]], plot_cells$cell_id, na.last=NA)
    n_highlight <- max(1L, min(length(finite_order), ceiling(length(finite_order) * (1 - opt$heatmap_highlight_quantile))))
    highlight_cells <- plot_cells[finite_order[seq_len(n_highlight)], , drop=FALSE]
    mean_response_score_diagnostics[[cy]] <- data.frame(
      cytokine=cy, n_cells=nrow(plot_cells), n_highlight=nrow(highlight_cells),
      zero_fraction=mean(abs(plot_cells[[cy]]) < 1e-12, na.rm=TRUE),
      score_q50=as.numeric(quantile(plot_cells[[cy]], .50, na.rm=TRUE)),
      score_q95=as.numeric(quantile(plot_cells[[cy]], .95, na.rm=TRUE)),
      score_max=max(plot_cells[[cy]], na.rm=TRUE), stringsAsFactors=FALSE)
    # In-panel callouts make explicit which visual layer each color bar describes.
    bg_values <- ifelse(is.finite(fld$grid$potential_in_support), fld$grid$potential_in_support, -Inf)
    bg_target <- fld$grid[which.max(bg_values), , drop=FALSE]
    point_target <- highlight_cells[which.max(highlight_cells[[cy]]), , drop=FALSE]
    x_span <- diff(range(fld$grid$UMAP_1, na.rm=TRUE)); y_span <- diff(range(fld$grid$UMAP_2, na.rm=TRUE))
    x_min <- min(fld$grid$UMAP_1, na.rm=TRUE); x_max <- max(fld$grid$UMAP_1, na.rm=TRUE)
    y_min <- min(fld$grid$UMAP_2, na.rm=TRUE); y_max <- max(fld$grid$UMAP_2, na.rm=TRUE)
    layer_callouts <- bind_rows(
      data.frame(layer="Background heatmap", label="Background heatmap\nLocal mean UCell response", x=x_min + .20*x_span, y=y_max - .08*y_span, xend=bg_target$UMAP_1, yend=bg_target$UMAP_2, color="#FFFFFF", stringsAsFactors=FALSE),
      data.frame(layer="Individual cells", label=paste0("Colored dots\nTop ", round((1-opt$heatmap_highlight_quantile)*100), "% individual-response cells"), x=x_max - .21*x_span, y=y_min + .09*y_span, xend=point_target$UMAP_1, yend=point_target$UMAP_2, color="#FDE725", stringsAsFactors=FALSE)
    )
    p <- ggplot() +
      geom_raster(data=fld$grid, aes(UMAP_1, UMAP_2, fill=potential_for_plot), interpolate=TRUE, alpha=.98) +
      geom_contour(data=fld$grid, aes(UMAP_1, UMAP_2, z=potential_for_plot), inherit.aes=FALSE,
                   color="white", linewidth=.28, bins=7, alpha=.72) +
      geom_point(data=plot_cells, aes(UMAP_1, UMAP_2), inherit.aes=FALSE, color="grey82", size=max(.28, opt$umap_point_size * .45), alpha=.58) +
      geom_point(data=highlight_cells, aes(UMAP_1, UMAP_2, color=.data[[cy]]), inherit.aes=FALSE, size=max(.45, opt$umap_point_size * 1.35), alpha=.95) +
      geom_segment(data=fld$arrows, aes(x=UMAP_1, y=UMAP_2, xend=xend, yend=yend, linetype="Cytokine-response field direction"),
                   inherit.aes=FALSE, color="#141414", linewidth=.46, alpha=.88,
                   arrow=arrow(length=grid::unit(1.6, "mm")), show.legend=TRUE)
    selected <- best_path[best_path$cytokine == cy, , drop=FALSE]
    if (nrow(selected)) {
      b <- selected$selected_branch[1]; lc <- long_curves[long_curves$branch == b, , drop=FALSE]
      if (nrow(lc) >= 3) {
        p <- p + geom_path(data=lc, aes(x=x, y=y, linetype="Pseudotime branch with strongest cytokine-response association"),
                           inherit.aes=FALSE, color=line_colors[[b]], linewidth=1.45, show.legend=TRUE)
        tip <- lc[nrow(lc), , drop=FALSE]; prev <- lc[nrow(lc)-1, , drop=FALSE]
        p <- p + geom_segment(data=data.frame(x=prev$x, y=prev$y, xend=tip$x, yend=tip$y), aes(x=x, y=y, xend=xend, yend=yend),
                              inherit.aes=FALSE, color=line_colors[[b]], linewidth=1.45, arrow=arrow(length=grid::unit(3.5, "mm")))
      }
    }
    p <- p +
      geom_segment(data=layer_callouts, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE,
                   linewidth=.48, linetype="dashed", color=layer_callouts$color, arrow=arrow(length=grid::unit(2.0, "mm")), show.legend=FALSE) +
      geom_label(data=layer_callouts, aes(x=x, y=y, label=label), inherit.aes=FALSE,
                 fill="white", color=layer_callouts$color, alpha=.92, fontface="bold", size=3.05, label.size=.32, show.legend=FALSE) +
      coord_equal(expand=FALSE) + theme_classic(base_size=12) +
      scale_fill_viridis_c(name="Kernel-smoothed\nmean UCell response", limits=potential_limits, breaks=potential_limits, oob=scales::squish, trans=potential_transform, na.value="#440154") +
      scale_color_viridis_c(name=paste0("Individual UCell response\n(top ", round((1-opt$heatmap_highlight_quantile)*100), "% cells)"), limits=individual_limits, breaks=individual_limits, oob=scales::squish, trans=individual_transform) +
      scale_linetype_manual(name="UMAP overlay", values=c("Cytokine-response field direction"="solid", "Pseudotime branch with strongest cytokine-response association"="solid")) +
      guides(linetype=guide_legend(override.aes=list(color=c("#141414", "#D55E00"), linewidth=c(.46, 1.45)), order=1), fill=guide_colorbar(order=2), color=guide_colorbar(order=3)) +
      theme(legend.position="bottom", legend.box="vertical") +
      labs(title=paste0(cy, ": UMAP mean-response heatmap and vector map"),
           subtitle=paste0("Background is the kernel-normalized local mean UCell response (sum K*score / sum K; not multiplied by cell abundance; bandwidth=", signif(fld$bandwidth, 3), "). Unsupported UMAP space is shown at the lowest response color. Black arrows show the same score gradient; colored points mark exactly the top ", round((1-opt$heatmap_highlight_quantile)*100), "% individual response scores. This is descriptive, not a predicted perturbation or causal force."),
           x="UMAP_1", y="UMAP_2")
    ggsave(file.path(ddir, paste0(safe_name(cy), "_UMAP_mean_response_heatmap_vector_map.pdf")), p, width=8.5, height=6.8, device=cairo_pdf)
    ggsave(file.path(ddir, paste0(safe_name(cy), "_UMAP_mean_response_heatmap_vector_map.svg")), p, width=8.5, height=6.8)
    print(p)
  }
  grDevices::dev.off()
  if (length(mean_response_grids)) write.csv(bind_rows(mean_response_grids), file.path(opt$outdir, "UMAP_mean_response_heatmap_vector_grid.csv"), row.names=FALSE)
  if (length(mean_response_score_diagnostics)) write.csv(bind_rows(mean_response_score_diagnostics), file.path(opt$outdir, "UMAP_mean_response_heatmap_score_diagnostics.csv"), row.names=FALSE)
}
cat("Completed post-Monocle3 scCRS trajectory analysis: ",opt$outdir,"\n",sep="")

