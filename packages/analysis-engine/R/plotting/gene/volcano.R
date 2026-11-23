#' Build the ATOP gene-score volcano plot.
#'
#' @return A ggplot object. The public signature is retained for the desktop
#' adapter and saved-project compatibility.
generate_volcano_plot <- function(data,
                                  p_value_cutoff = 0.05,
                                  score_cutoff = 0.2,
                                  label_genes = NULL,
                                  col_up = "#C02635",
                                  col_down = "#00609D",
                                  col_ns = "#CECECE",
                                  base_font_size = 7,
                                  plot_title = "Global gene-level evidence",
                                  show_count_labels = FALSE,
                                  font_family = "Helvetica") {
    atop_require_columns(
        data,
        c("Gene", "GeneNegativeScore", "P_negative", "GenePositiveScore", "P_positive"),
        "Gene summary"
    )
    if (!is.finite(p_value_cutoff) || p_value_cutoff <= 0 || p_value_cutoff > 1) {
        stop("p_value_cutoff must be between 0 and 1.", call. = FALSE)
    }

    negative <- dplyr::transmute(
        dplyr::filter(data, .data$GeneNegativeScore < 0),
        Gene = .data$Gene, Score = .data$GeneNegativeScore, P = .data$P_negative
    )
    positive <- dplyr::transmute(
        dplyr::filter(data, .data$GenePositiveScore > 0),
        Gene = .data$Gene, Score = .data$GenePositiveScore, P = .data$P_positive
    )
    plot_data <- dplyr::filter(dplyr::bind_rows(negative, positive), is.finite(.data$Score), is.finite(.data$P))
    if (nrow(plot_data) == 0L) stop("No finite gene scores and P values are available for the volcano plot.", call. = FALSE)

    positive_p <- plot_data$P[plot_data$P > 0]
    p_floor <- if (length(positive_p) == 0L) .Machine$double.xmin else min(positive_p) / 10
    plot_data$P[plot_data$P <= 0] <- p_floor
    plot_data$group <- "NS"
    plot_data$group[plot_data$Score > score_cutoff & plot_data$P < p_value_cutoff] <- "Up"
    plot_data$group[plot_data$Score < -score_cutoff & plot_data$P < p_value_cutoff] <- "Down"
    plot_data$group <- factor(plot_data$group, levels = c("Down", "NS", "Up"))
    y_limit <- max(-log10(plot_data$P), na.rm = TRUE) * 1.2
    label_data <- if (length(label_genes %||% character()) > 0L) {
        dplyr::filter(plot_data, .data$Gene %in% label_genes)
    } else plot_data[FALSE, , drop = FALSE]

    palette <- c(Down = col_down, NS = col_ns, Up = col_up)
    background <- vapply(palette, function(x) {
        rgb <- grDevices::col2rgb(x) / 255
        grDevices::rgb(t(rgb + 0.6 * (1 - rgb)))
    }, character(1))
    plot_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Score, y = -log10(.data$P), colour = .data$group)) +
        ggplot2::geom_point(size = 0.8, stroke = 0) +
        ggplot2::geom_hline(yintercept = -log10(p_value_cutoff), linetype = 2, linewidth = 0.2, colour = "#999999") +
        ggplot2::scale_colour_manual(values = background, name = "Group") +
        ggplot2::scale_y_continuous(limits = c(0, max(y_limit, -log10(p_value_cutoff) * 1.2, 1)), expand = c(0, 0)) +
        ggplot2::labs(x = "Gene score", y = "−log₁₀(P)", title = plot_title) +
        ggplot2::theme_classic(base_size = base_font_size, base_family = font_family) +
        atop_plot_text_theme(base_font_size, font_family) +
        ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal",
            legend.key.width = grid::unit(3, "mm"))
    if (score_cutoff > 0) {
        plot_obj <- plot_obj + ggplot2::geom_vline(xintercept = c(-score_cutoff, score_cutoff), linetype = 2, linewidth = 0.2, colour = "#999999")
    }
    if (nrow(label_data) > 0L) {
        plot_obj <- plot_obj + ggplot2::geom_point(data = label_data,
            colour = unname(palette[as.character(label_data$group)]), size = 1.5, stroke = 0, show.legend = FALSE) +
            ggrepel::geom_text_repel(
            data = label_data, ggplot2::aes(label = .data$Gene), inherit.aes = TRUE,
            seed = 123, size = atop_font_size_mm(base_font_size), family = font_family, colour = "#000000",
            min.segment.length = 0, max.overlaps = Inf, segment.colour = "#4D4D4D"
        )
    }
    if (isTRUE(show_count_labels)) {
        plot_obj <- plot_obj + ggplot2::annotate(
            "text", x = max(plot_data$Score) * 0.8, y = y_limit * 0.9,
            label = paste0("Up = ", sum(plot_data$group == "Up")), colour = col_up,
            family = font_family, size = atop_font_size_mm(base_font_size)
        ) + ggplot2::annotate(
            "text", x = min(plot_data$Score) * 0.8, y = y_limit * 0.9,
            label = paste0("Down = ", sum(plot_data$group == "Down")), colour = col_down,
            family = font_family, size = atop_font_size_mm(base_font_size)
        )
    }
    plot_obj
}
