#' Compare normalized read counts for a gene's sgRNAs across two conditions.
generate_sgrna_paired_plot <- function(normalized_counts_df,
                                      target_gene_id,
                                      sgrna_id_col,
                                      gene_col,
                                      condition_1_col,
                                      condition_2_col,
                                      condition_1_label = NULL,
                                      condition_2_label = NULL,
                                      color_condition_1 = "#00609D",
                                      color_condition_2 = "#C02635",
                                      connecting_line_color = "#7F7F7F",
                                      connecting_line_size = 0.5,
                                      connecting_line_alpha = 0.7,
                                      point_size = 2.5,
                                      point_alpha = 0.8,
                                      plot_title = NULL,
                                      y_axis_label = "Normalized read count",
                                      base_font_size = 8,
                                      font_family = "Helvetica") {
    atop_require_columns(normalized_counts_df, c(sgrna_id_col, gene_col, condition_1_col, condition_2_col), "Normalized count data")
    rows <- dplyr::filter(normalized_counts_df, .data[[gene_col]] == target_gene_id)
    if (nrow(rows) == 0L) stop("Target gene was not found in normalized counts: ", target_gene_id, call. = FALSE)
    rows <- dplyr::transmute(
        rows, sgrna = .data[[sgrna_id_col]],
        condition_1 = as.numeric(.data[[condition_1_col]]), condition_2 = as.numeric(.data[[condition_2_col]])
    )
    rows <- stats::na.omit(rows)
    if (nrow(rows) == 0L) stop("The selected gene has no complete paired sgRNA measurements.", call. = FALSE)
    labels <- c(condition_1_label %||% condition_1_col, condition_2_label %||% condition_2_col)
    long_data <- tidyr::pivot_longer(rows, c("condition_1", "condition_2"), names_to = "group", values_to = "value")
    long_data$group <- factor(long_data$group, levels = c("condition_1", "condition_2"), labels = labels)
    ggplot2::ggplot(long_data, ggplot2::aes(x = .data$group, y = .data$value, colour = .data$group)) +
        ggplot2::geom_line(ggplot2::aes(group = .data$sgrna), colour = connecting_line_color, linewidth = connecting_line_size, alpha = connecting_line_alpha) +
        ggplot2::geom_point(size = point_size, alpha = point_alpha) +
        ggplot2::scale_colour_manual(name = NULL, values = stats::setNames(c(color_condition_1, color_condition_2), labels)) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.10))) +
        ggplot2::labs(y = y_axis_label, x = NULL, title = plot_title %||% paste("sgRNA paired plot for", target_gene_id)) +
        ggpubr::theme_pubr(base_size = base_font_size, base_family = font_family) +
        atop_plot_text_theme(base_font_size, font_family) +
        ggplot2::theme(legend.position = "top", legend.key.width = grid::unit(3, "mm"), plot.margin = ggplot2::margin(3, 4, 3, 4))
}
