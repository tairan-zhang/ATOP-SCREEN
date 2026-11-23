#' Build the GSEA lollipop summary used by the pathway results view.
atop_gsea_legend_breaks <- function(values) {
    values <- values[is.finite(values)]
    if (length(values) == 0L) return(NULL)
    domain <- range(values)
    candidates <- unique(signif(as.numeric(stats::quantile(values, c(0, 0.5, 1), names = FALSE)), 2))
    candidates[candidates >= domain[[1]] & candidates <= domain[[2]]]
}

generate_gsea_lollipop_plot <- function(gsea_data_input_df,
                                        analysis_type_label,
                                        num_pathways_to_plot,
                                        plot_title_main,
                                        bar_fill_low_hex = "#EFC8CC",
                                        bar_fill_high_hex = "#C02635",
                                        circle_fill_low_hex = "#FFFFFF",
                                        circle_fill_high_hex = "#C02635",
                                        base_font_size = 6.5,
                                        bar_line_width = 3.4,
                                        point_size_min = 2.4,
                                        point_size_max = 5.4,
                                        point_stroke = 0.32,
                                        font_family = "Helvetica",
                                        pathway_wrap_width = 36,
                                        positive_color = NULL,
                                        negative_color = NULL) {
    if (is.null(gsea_data_input_df) || !nrow(gsea_data_input_df)) stop("No pathway results to plot.")
    table <- gsea_data_input_df
    q_column <- intersect(c("qvalue", "qvalues", "p.adjust", "padj", "FDR"), names(table))
    id_column <- intersect(c("ID", "Description"), names(table))
    if (!length(q_column) || !length(id_column) || !"NES" %in% names(table)) {
        stop("Pathway results require an identifier, NES and a q-value or adjusted P column.")
    }
    table$SourcePathwayID <- as.character(table[[id_column[1L]]])
    table$QValue <- as.numeric(table[[q_column[1L]]])
    table <- table[is.finite(table$QValue) & table$QValue >= 0 & table$QValue <= 1 & is.finite(table$NES), , drop = FALSE]
    table$Direction <- ifelse(table$NES >= 0, "Positive", "Negative")
    if (tolower(analysis_type_label) %in% c("positive", "negative")) {
        table <- table[tolower(table$Direction) == tolower(analysis_type_label), , drop = FALSE]
    }
    table <- table[order(table$QValue, -abs(table$NES), table$SourcePathwayID), , drop = FALSE]
    table <- head(table[!duplicated(table$SourcePathwayID), , drop = FALSE], max(1L, as.integer(num_pathways_to_plot)))
    if (!nrow(table)) stop("No pathways remain in the selected direction.")
    table$effect <- -log10(pmax(table$QValue, .Machine$double.xmin)) * ifelse(table$Direction == "Negative", -1, 1)
    display <- if ("Description" %in% names(table)) as.character(table$Description) else table$SourcePathwayID
    missing_label <- is.na(display) | !nzchar(display)
    display[missing_label] <- table$SourcePathwayID[missing_label]
    table$label <- stringr::str_wrap(gsub("_", " ", display, fixed = TRUE), width = pathway_wrap_width)
    lines <- lengths(strsplit(table$label, "\n", fixed = TRUE))
    centres <- cumsum(lines) - lines / 2
    table$row <- max(centres) - centres + 1
    palette <- c(Positive = positive_color %||% bar_fill_low_hex,
                 Negative = negative_color %||% "#BFD7E6")
    tick_label <- if (q_column[1L] %in% c("qvalue", "qvalues")) "q" else "adjusted P"
    ggplot2::ggplot(table, ggplot2::aes(y = .data$row, x = .data$effect, fill = .data$Direction)) +
        ggplot2::geom_col(orientation = "y", width = 0.72) +
        ggplot2::geom_vline(xintercept = 0, linewidth = atop_axis_linewidth) +
        ggplot2::geom_point(ggplot2::aes(size = abs(.data$NES)), shape = 21, stroke = 0) +
        ggplot2::geom_text(ggplot2::aes(label = .data$label, hjust = ifelse(.data$effect < 0, 1, 0)),
            x = 0, size = atop_font_size_mm(base_font_size), family = font_family, lineheight = 0.9) +
        ggplot2::scale_fill_manual(values = palette, guide = "none") +
        ggplot2::scale_size_continuous(name = "|NES|", range = c(point_size_min, point_size_max),
            breaks = atop_gsea_legend_breaks(abs(table$NES)),
            guide = ggplot2::guide_legend(override.aes = list(fill = "grey65"))) +
        ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(add = 0.7)) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.05)) +
        ggplot2::labs(x = paste0("Signed −log10(", tick_label, ")"), y = NULL, title = plot_title_main) +
        ggplot2::theme_classic(base_size = base_font_size, base_family = font_family) +
        atop_plot_text_theme(base_font_size, font_family) +
        ggplot2::theme(plot.title.position = "plot", legend.position = "bottom", legend.key.height = grid::unit(3, "mm"),
            axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
}
