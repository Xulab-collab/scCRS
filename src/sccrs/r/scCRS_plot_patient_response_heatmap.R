#!/usr/bin/env Rscript
# Patient-level scCRS cytokine-response heatmap.
suppressPackageStartupMessages({ library(optparse); library(dplyr); library(pheatmap); library(Cairo) })
opts <- list(
  make_option('--input', type='character', help='patient_cytokine_scores_with_group.csv'),
  make_option('--outdir', type='character', default='sccrs_patient_heatmap'),
  make_option('--score-col', dest='score_col', type='character', default='overall_score_0_100'),
  make_option('--patient-col', dest='patient_col', type='character', default='patient_id'),
  make_option('--group-col', dest='group_col', type='character', default='group'),
  make_option('--cytokines', type='character', default=NULL, help='Optional comma-separated cytokines'),
  make_option('--group-order', dest='group_order', type='character', default=NULL, help='Optional comma-separated group order'),
  make_option('--width', type='double', default=8.0), make_option('--height', type='double', default=13.0)
)
a <- parse_args(OptionParser(option_list=opts)); if (is.null(a$input)) stop('Required: --input')
dir.create(a$outdir, recursive=TRUE, showWarnings=FALSE)
x <- read.csv(a$input, check.names=FALSE, stringsAsFactors=FALSE)
need <- c(a$patient_col,a$group_col,'cytokine',a$score_col); miss <- setdiff(need,names(x)); if(length(miss)) stop('Missing columns: ',paste(miss,collapse=', '))
x <- x %>% transmute(patient_id=as.character(.data[[a$patient_col]]), group=as.character(.data[[a$group_col]]), cytokine=as.character(cytokine), score=as.numeric(.data[[a$score_col]])) %>% filter(complete.cases(.))
if (!is.null(a$cytokines)) { keep <- trimws(strsplit(a$cytokines,',',fixed=TRUE)[[1]]); x <- filter(x, cytokine %in% keep) }
ann <- distinct(x,patient_id,group); if(anyDuplicated(ann$patient_id)) stop('Each patient must have exactly one group.')
if (!is.null(a$group_order)) { lev <- trimws(strsplit(a$group_order,',',fixed=TRUE)[[1]]); ann$group <- factor(ann$group,levels=lev) } else ann$group <- factor(ann$group)
ann <- arrange(ann,group,patient_id); mat <- with(x,tapply(score,list(cytokine,patient_id),median,na.rm=TRUE)); mat <- mat[sort(rownames(mat)),ann$patient_id,drop=FALSE]
if(anyNA(mat)) stop('Missing patient-by-cytokine combinations; filter inputs or complete the score table before plotting.')
z <- t(scale(t(mat))); z[!is.finite(z)] <- 0; z <- pmax(pmin(z,2.5),-2.5)
write.csv(z,file.path(a$outdir,'row_zscore_cytokine_by_patient_matrix.csv'))
tidy <- as.data.frame(as.table(z)); names(tidy) <- c('cytokine','patient_id','row_zscore'); tidy$group <- as.character(ann$group[match(tidy$patient_id,ann$patient_id)]); write.csv(tidy,file.path(a$outdir,'row_zscore_cytokine_by_patient_plotting_data.csv'),row.names=FALSE); write.csv(ann,file.path(a$outdir,'row_zscore_heatmap_patient_group_annotation.csv'),row.names=FALSE)
pal <- setNames(c('#6BAED6','#74C476','#FB6A4A','#9E9AC8','#FDAE6B')[seq_along(levels(ann$group))],levels(ann$group)); annotation_col <- data.frame(Group=ann$group,row.names=ann$patient_id)
CairoPDF(file.path(a$outdir,'scCRS_patient_cytokine_response_row_zscore_heatmap.pdf'),width=a$width,height=a$height)
pheatmap(z,color=colorRampPalette(c('#2166AC','#F7F7F7','#B2182B'))(120),breaks=seq(-2.5,2.5,length.out=121),cluster_rows=TRUE,cluster_cols=FALSE,annotation_col=annotation_col,annotation_colors=list(Group=pal),cellwidth=12,cellheight=12,border_color='grey85',fontsize_row=10.5,fontsize_col=9.5,treeheight_row=40,treeheight_col=0,angle_col=45,main='Patient cytokine-response heatmap: row-standardized scores')
dev.off(); cat('Wrote patient heatmap to ',a$outdir,'\n',sep='')
