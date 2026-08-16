# Use the available overall-score column and a Cairo PDF device so Greek cytokine labels render correctly.
source_path <- "C:/Users/PC_XH/Documents/Codex/2026-08-11/dictionary-of-immune-responses-to-cytokines/work/create_patient_group_violin_and_heatmap.R"
workflow <- readLines(source_path, warn = FALSE)
workflow <- gsub("raw_rank_score", "overall_score_0_100", workflow, fixed = TRUE)
workflow <- gsub("Mean raw rank score", "Mean overall response score (0-100)", workflow, fixed = TRUE)
workflow <- gsub("Raw rank score (uncalibrated)", "Overall response score (0-100; uncalibrated)", workflow, fixed = TRUE)
workflow <- gsub("raw rank scores", "overall response scores (0-100)", workflow, fixed = TRUE)
workflow <- gsub("grDevices::pdf(heatmap_pdf", "CairoPDF(heatmap_pdf", workflow, fixed = TRUE)
eval(parse(text = workflow), envir = globalenv())
