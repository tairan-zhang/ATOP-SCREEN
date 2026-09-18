suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    library(shiny)
})
for (path in list.files("R", pattern = "[.]R$", recursive = TRUE, full.names = TRUE)) source(path)
for (path in list.files(".", pattern = "[.]R$", recursive = TRUE, full.names = TRUE)) parse(path)

data <- data.frame(
    Gene = c("A", "B", "C"), GenePositiveScore = c(2, 0.5, 1),
    GeneNegativeScore = c(-1, -2, -0.5), P_positive = c(0.01, NA, 0),
    P_negative = c(0.3, 0.001, NA), GenePositiveRank = c(1, 3, 2),
    GeneNegativeRank = c(2, 1, 3)
)
volcano <- generate_volcano_plot(data)
stopifnot(inherits(volcano, "ggplot"), all(is.finite(volcano$data$P)), all(volcano$data$P > 0))
invisible(ggplot_build(volcano))
ranking <- generate_ranking_plot(data)
invisible(ggplot_build(ranking))
data$P_positive <- data$P_negative <- 0
stopifnot(all(is.finite(generate_volcano_plot(data)$data$P)))
data$P_positive <- data$P_negative <- NA_real_
stopifnot(inherits(tryCatch(generate_volcano_plot(data), error = identity), "error"))

counts <- data.frame(Gene = rep("A", 3), gRNA = letters[1:3], treatment = 1:3, control = 3:1)
paired <- generate_sgrna_paired_plot(counts, "A", "gRNA", "Gene", "treatment", "control")
invisible(ggplot_build(paired))
ranked <- prepare_ranked_gene_list_gsea(
    data.frame(Gene = c("A", "A", "B", "C"), score = c(1, 2, Inf, -1)),
    "score", "Gene", NULL
)
stopifnot(identical(ranked, c(A = 2, C = -1)))
stopifnot(is_valid_hex_or_name("#12ABCD"), !is_valid_hex_or_name("not-a-color"))
stopifnot(validate_hex_inputs(character(), character(), NULL))
stopifnot(select_col(character(), "Gene") %>% is.null())
cat("PASS: module loading, R syntax, volcano/ranking/paired plots, GSEA ranking, and UI helpers.\n")
