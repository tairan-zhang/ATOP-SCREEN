#' Build a ranked positive or negative ATOP gene-score chart.
generate_ranking_plot <- function(data,
                                  type = "positive",
                                  top_n = 10,
                                  gradient_low = if (type == "negative") "#BFD7E6" else "#EFC8CC",
                                  gradient_high = if (type == "negative") "#00609D" else "#C02635",
                                  base_font_size = 6.8,
                                  plot_title = NULL,
                                  known_genes = NULL,
                                  font_family = "Helvetica") {
    type <- match.arg(type, c("positive", "negative"))
    top_n <- max(1L, as.integer(top_n))
    columns <- if (type == "positive") {
        c(rank = "GenePositiveRank", score = "GenePositiveScore")
    } else {
        c(rank = "GeneNegativeRank", score = "GeneNegativeScore")
    }
    atop_require_columns(data, c("Gene", unname(columns)), paste(type, "gene ranking data"))
    plot_data <- dplyr::filter(data, !grepl("^hsa-mir-", .data$Gene, ignore.case = TRUE))
    plot_data <- dplyr::transmute(
        dplyr::slice_head(dplyr::arrange(plot_data, .data[[columns[["rank"]]]]), n = top_n),
        gene = .data$Gene,
        rank = .data[[columns[["rank"]]]],
        score = .data[[columns[["score"]]]],
        fill_score = abs(.data[[columns[["score"]]]]),
        is_known = .data$Gene %in% (known_genes %||% character())
    )
    if (nrow(plot_data) == 0L) stop("No genes are available for the ranking plot.", call. = FALSE)
    plot_data$gene <- factor(plot_data$gene, levels = rev(plot_data$gene))
    title <- plot_title %||% paste("Top", type, "gene scores")
    span <- max(abs(plot_data$score), .Machine$double.eps)
    reverse_axis <- type == "negative" && all(plot_data$score <= 0)
    plot_data$position <- if (reverse_axis) -plot_data$score else plot_data$score
    plot_data$rank_label <- paste0("#", plot_data$rank, ifelse(plot_data$is_known, "*", ""))
    plot_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$position, y = .data$gene, fill = .data$fill_score)) +
        ggplot2::geom_col(width = 0.68) +
        ggplot2::geom_text(ggplot2::aes(label = .data$rank_label), x = 0, hjust = 0,
            size = atop_font_size_mm(base_font_size), family = font_family, fontface = "plain") +
        ggplot2::geom_text(ggplot2::aes(x = .data$position + ifelse(.data$position < 0, -1, 1) * span * 0.025,
                hjust = ifelse(.data$position < 0, 1, 0), label = sprintf("%.3f", .data$score)), size = atop_font_size_mm(base_font_size), family = font_family) +
        ggplot2::scale_fill_gradient(low = gradient_low, high = gradient_high, guide = "none") +
        ggplot2::scale_x_continuous(limits = c(min(0, plot_data$position) * 1.27, max(0, plot_data$position) * 1.27 + span * 0.02), expand = c(0, 0),
            labels = if (reverse_axis) function(x) format(-x, trim = TRUE) else ggplot2::waiver()) +
        ggplot2::labs(title = title, x = paste("Gene", type, "score"), y = NULL) +
        ggplot2::theme_classic(base_size = base_font_size, base_family = font_family) +
        atop_plot_text_theme(base_font_size, font_family) +
        ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = base_font_size, face = "plain"),
            panel.grid.major.x = ggplot2::element_line(colour = "grey88", linewidth = 0.35))
    attr(plot_obj, "source_data") <- plot_data
    plot_obj
}
