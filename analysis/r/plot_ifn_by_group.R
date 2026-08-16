suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(Cairo)
})

input_path <- "F:/1000G-20130502-hg19/cell_type_cytokine_scores_with_group.csv"
output_dir <- "C:/Users/PC_XH/Documents/Codex/2026-08-11/dictionary-of-immune-responses-to-cytokines/outputs/cytokine_response_figures"
output_pdf <- file.path(output_dir, "pSS_IFN_response_by_group_violin_plots.pdf")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scores <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
group_levels <- c("Normal", "Non_High", "High")
group_labels <- c("Normal", "Low globulin", "High globulin")
cell_levels <- c("B cells", "CD4+ T", "CD8+ T", "NK", "\u03b3\u03b4T", "ILCs", "myeloid cell", "MAIT cells")
ifn_factors <- c("IFN-\u03b11", "IFN-\u03b2", "IFN-\u03b3", "IFN-\u03b5", "IFN-\u03ba", "IFN-\u03bb2")
group_palette <- c("Normal" = "#6BAED6", "Low globulin" = "#74C476", "High globulin" = "#FB6A4A")

scores <- scores %>%
  filter(cytokine %in% ifn_factors, cell_type %in% cell_levels) %>%
  mutate(
    group = factor(group, levels = group_levels, labels = group_labels),
    cell_type = factor(cell_type, levels = cell_levels),
    cytokine = factor(cytokine, levels = ifn_factors)
  )

if (nrow(scores) == 0L) stop("No IFN score rows were found. Check the input column labels and encodings.")

coverage <- scores %>%
  group_by(cytokine, cell_type) %>%
  summarise(n_patients = n_distinct(patient_id), .groups = "drop")
write.csv(coverage, file.path(output_dir, "IFN_response_coverage_by_cell_type.csv"), row.names = FALSE)

make_ifn_plot <- function(factor_name) {
  data <- scores %>% filter(cytokine == factor_name)
  present_cells <- unique(as.character(data$cell_type))
  missing_cells <- setdiff(cell_levels, present_cells)
  y_limits <- range(data$raw_rank_score, na.rm = TRUE)
  y_pad <- max(0.08, diff(y_limits) * 0.08)
  y_annotate <- mean(y_limits)
  note_data <- data.frame(
    cell_type = factor(missing_cells, levels = cell_levels),
    group = factor(rep("Low globulin", length(missing_cells)), levels = group_labels),
    raw_rank_score = rep(y_annotate, length(missing_cells)),
    note = rep("No dictionary\nsignature", length(missing_cells))
  )

  ggplot(data, aes(x = group, y = raw_rank_score, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.62, colour = "grey35", linewidth = 0.3, na.rm = TRUE) +
    geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", colour = "black", linewidth = 0.35, na.rm = TRUE) +
    geom_jitter(aes(colour = group), width = 0.11, height = 0, size = 1.45, alpha = 0.85, show.legend = FALSE, na.rm = TRUE) +
    stat_summary(fun = median, geom = "point", shape = 95, size = 5, colour = "black", na.rm = TRUE) +
    geom_text(data = note_data, aes(x = group, y = raw_rank_score, label = note), inherit.aes = FALSE, colour = "grey45", size = 3) +
    facet_wrap(~cell_type, ncol = 4, drop = FALSE) +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    scale_colour_manual(values = group_palette, drop = FALSE) +
    coord_cartesian(ylim = c(y_limits[1] - y_pad, y_limits[2] + y_pad)) +
    labs(
      title = paste0(as.character(factor_name), " response across groups and immune cell types"),
      subtitle = "Each point is one patient; white boxes show interquartile range and black dash marks the median",
      x = NULL, y = "Raw rank score (uncalibrated)",
      caption = "Groups: Normal n=3, Low globulin n=7, High globulin n=11. No significance stars are shown because the Normal group is small."
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 8),
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
}

CairoPDF(output_pdf, width = 14, height = 8.8, onefile = TRUE)
for (factor_name in ifn_factors) print(make_ifn_plot(factor_name))
dev.off()
cat("Wrote", output_pdf, "with", length(ifn_factors), "pages.\n")
