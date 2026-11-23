# ATOP-SCREEN running-plot typography and panel composition.

simplify_gsea_running_curve_data_atop <- function(data) {
    required_columns <- c("x", "runningScore", "position")
    if (!is.data.frame(data) || !all(required_columns %in% names(data)) || nrow(data) < 3L) return(data)
    groups <- if ("id" %in% names(data)) split(seq_len(nrow(data)), data$id, drop = TRUE) else list(seq_len(nrow(data)))
    retained <- unlist(lapply(groups, function(indices) {
        positions <- suppressWarnings(as.numeric(data$position[indices]))
        hits <- which(is.finite(positions) & positions != 0)
        local_indices <- unique(c(1L, hits - 1L, hits, hits + 1L, length(indices)))
        indices[local_indices[local_indices >= 1L & local_indices <= length(indices)]]
    }), use.names = FALSE)
    data[sort(unique(retained)), , drop = FALSE]
}

style_gsea_running_plot_atop <- function(plot_object,
                                         direction = c("positive", "negative"),
                                         pathway_rows = NULL,
                                         subplot_type = 3,
                                         pathway_count = 1,
                                         legend_position = c(0.85, 0.8),
                                         base_font_size = 11,
                                         show_nes = TRUE,
                                         show_pvalue = TRUE,
                                         pvalue_type = c("p.adjust", "pvalue"),
                                         ranked_gene_list = NULL,
                                         show_rank_explanation = TRUE,
                                         rank_positive_label = "Top-ranked genes (positively correlated)",
                                         rank_negative_label = "Bottom-ranked genes (negatively correlated)",
                                         font_family = "Helvetica", plot_title = NULL) {
    direction <- match.arg(direction)
    pvalue_type <- match.arg(pvalue_type)
    subplot_type <- max(1L, min(3L, as.integer(subplot_type)))
    pathway_count <- max(1L, as.integer(pathway_count))
    annotation_size <- atop_font_size_mm(base_font_size)

    curve_colour <- if (direction == "positive") "#C02635" else "#00609D"
    panel_title <- if (!is.null(pathway_rows) && nrow(pathway_rows) > 0L) {
        pathway_values <- if ("OriginalDescription" %in% names(pathway_rows)) {
            as.character(pathway_rows$OriginalDescription)
        } else if ("Description" %in% names(pathway_rows)) {
            as.character(pathway_rows$Description)
        } else {
            as.character(pathway_rows$ID)
        }
        missing_values <- is.na(pathway_values) | !nzchar(pathway_values)
        pathway_values[missing_values] <- as.character(pathway_rows$ID)[missing_values]
        stringr::str_wrap(
            paste(format_gsea_pathway_name_atop(pathway_values, width = NULL), collapse = "  /  "),
            width = 72
        )
    } else {
        "Selected pathway enrichment"
    }

    if (!is.null(plot_title)) panel_title <- plot_title
    top <- if (inherits(plot_object, "patchwork")) plot_object[[1]] else plot_object
    for (layer_index in seq_along(top$layers)) {
        if (inherits(top$layers[[layer_index]]$geom, "GeomLine") && inherits(top$layers[[layer_index]]$data, "waiver")) {
            top$layers[[layer_index]]$data <- simplify_gsea_running_curve_data_atop(top$data)
        }
    }
    top <- top +
        ggplot2::labs(title = NULL) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.10))) +
        ggplot2::theme(
            text = ggplot2::element_text(family = font_family, colour = "#000000"),
            plot.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_text(size = base_font_size, colour = "#000000"),
            axis.title = ggplot2::element_text(size = base_font_size, colour = "#000000"),
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            axis.line.x = atop_axis_line(),
            axis.line.y.left = atop_axis_line(),
            axis.ticks = atop_axis_line(),
            legend.position = if (pathway_count > 1L) legend_position else "none",
            legend.title = ggplot2::element_text(size = base_font_size, face = "plain", colour = "#000000"),
            legend.text = ggplot2::element_text(size = base_font_size, colour = "#000000"),
            legend.key.height = grid::unit(2.4, "mm"),
            legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
            plot.margin = ggplot2::margin(2, 3, 0, 3)
        )

    for (layer_index in seq_along(top$layers)) {
        if (inherits(top$layers[[layer_index]]$geom, "GeomHline")) {
            top$layers[[layer_index]]$aes_params$linewidth <- atop_axis_linewidth
        }
        if (inherits(top$layers[[layer_index]]$geom, "GeomLine")) {
            if (pathway_count == 1L) {
                top$layers[[layer_index]]$mapping$colour <- NULL
                top$layers[[layer_index]]$aes_params$colour <- curve_colour
            }
            top$layers[[layer_index]]$aes_params$linewidth <- 0.62
        }
        if (inherits(top$layers[[layer_index]]$geom, "GeomTextRepel") || inherits(top$layers[[layer_index]]$geom, "GeomText")) {
            top$layers[[layer_index]]$aes_params$size <- annotation_size
            top$layers[[layer_index]]$aes_params$family <- font_family
        }
    }

    hits <- NULL
    ranked <- NULL
    if (subplot_type >= 2L && inherits(plot_object, "patchwork")) {
        hits <- plot_object[[2]] + ggplot2::theme(
            text = ggplot2::element_text(family = font_family, colour = "#000000"),
            axis.text = ggplot2::element_text(size = base_font_size, colour = "#000000"),
            axis.title = ggplot2::element_text(size = base_font_size, colour = "#000000"),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            axis.line.x = atop_axis_line(),
            axis.line.y.left = atop_axis_line(),
            axis.ticks = atop_axis_line(),
            plot.margin = if (subplot_type == 2L) ggplot2::margin(0, 3, 2, 3) else ggplot2::margin(0, 3, 0, 3)
        )
        for (i in seq_along(hits$layers)) {
            if (inherits(hits$layers[[i]]$geom, "GeomSegment")) {
                hits$layers[[i]]$aes_params$linewidth <- atop_axis_linewidth
                hits$layers[[i]]$aes_params$lineend <- "butt"
            }
        }
    }
    if (subplot_type == 3L && inherits(plot_object, "patchwork")) {
        ranked <- plot_object[[3]] + ggplot2::theme(
            text = ggplot2::element_text(family = font_family, colour = "#000000"),
            axis.text = ggplot2::element_text(size = base_font_size, colour = "#000000"),
            axis.title = ggplot2::element_text(size = base_font_size, colour = "#000000"),
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            axis.line.x = atop_axis_line(),
            axis.line.y.left = atop_axis_line(),
            axis.ticks = atop_axis_line(),
            plot.margin = ggplot2::margin(0, 3, 2, 3)
        )

        for (i in seq_along(ranked$layers)) {
            if (inherits(ranked$layers[[i]]$geom, "GeomHline")) {
                ranked$layers[[i]]$aes_params$linewidth <- atop_axis_linewidth
            }
        }
        rank_values <- suppressWarnings(as.numeric(ranked_gene_list))
        rank_values <- rank_values[is.finite(rank_values)]
        if (isTRUE(show_rank_explanation) && length(rank_values) > 1L) {
            rank_range <- range(rank_values, na.rm = TRUE)
            rank_span <- diff(rank_range)
            if (!is.finite(rank_span) || rank_span <= 0) rank_span <- max(abs(rank_range), 1)
            rank_count <- length(rank_values)
            zero_candidates <- which(rank_values <= 0)
            zero_rank <- if (length(zero_candidates) > 0L) zero_candidates[[1]] else which.min(abs(rank_values))
            zero_in_first_half <- zero_rank <= rank_count / 2
            zero_offset <- rank_span * 0.08
            zero_y <- if (zero_in_first_half) {
                max(rank_range[[1]] + rank_span * 0.12, -zero_offset)
            } else {
                min(rank_range[[2]] - rank_span * 0.12, zero_offset)
            }

            ranked <- ranked +
                ggplot2::annotate(
                    "text", x = 1 + rank_count * 0.018,
                    y = rank_range[[2]] - rank_span * 0.08,
                    label = rank_positive_label, hjust = 0, vjust = 1,
                    family = font_family, fontface = "plain",
                    size = annotation_size, colour = "#C02635"
                ) +
                ggplot2::annotate(
                    "text", x = rank_count - rank_count * 0.018,
                    y = rank_range[[1]] + rank_span * 0.08,
                    label = rank_negative_label, hjust = 1, vjust = 0,
                    family = font_family, fontface = "plain",
                    size = annotation_size, colour = "#00609D"
                ) +
                ggplot2::annotate(
                    "label", x = zero_rank, y = zero_y,
                    label = paste0("Zero cross at ", format(zero_rank, big.mark = ",")),
                    family = font_family, size = annotation_size,
                    linewidth = 0, fill = "white", colour = "#555555", vjust = 0.5
                )
        }
    }

    subtitle <- NULL
    if (!is.null(pathway_rows) && nrow(pathway_rows) > 0L && (isTRUE(show_nes) || isTRUE(show_pvalue))) {
        p_label <- if (pvalue_type == "pvalue") "P" else "adjusted P"
        stat_lines <- vapply(seq_len(nrow(pathway_rows)), function(row_index) {
            stat_parts <- character()
            if (isTRUE(show_nes)) {
                stat_parts <- c(stat_parts, sprintf("NES %.2f", as.numeric(pathway_rows$NES[[row_index]])))
            }
            if (isTRUE(show_pvalue)) {
                stat_parts <- c(
                    stat_parts,
                    paste(p_label, format_gsea_p_atop(as.numeric(pathway_rows[[pvalue_type]][[row_index]])))
                )
            }
            paste(stat_parts, collapse = "  |  ")
        }, character(1))
        subtitle <- paste(stat_lines, collapse = "\n")
    }

    combined <- if (subplot_type == 1L) {
        patchwork::wrap_plots(top)
    } else if (subplot_type == 2L) {
        (top / hits) + patchwork::plot_layout(heights = c(0.78, 0.22))
    } else {
        (top / hits / ranked) + patchwork::plot_layout(heights = c(0.58, 0.16, 0.26))
    }

    has_subtitle <- !is.null(subtitle) && nzchar(subtitle)
    title_grobs <- list(grid::textGrob(
        panel_title,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(if (has_subtitle) 0.62 else 0.5, "npc"),
        gp = grid::gpar(
            fontfamily = font_family,
            fontsize = base_font_size,
            fontface = "plain",
            col = "#000000"
        )
    ))
    if (has_subtitle) {
        title_grobs[[2]] <- grid::textGrob(
            subtitle,
            x = grid::unit(0.5, "npc"),
            y = grid::unit(0.25, "npc"),
            gp = grid::gpar(fontfamily = font_family, fontsize = base_font_size, col = "#000000")
        )
    }
    title_panel <- patchwork::wrap_elements(full = do.call(grid::grobTree, title_grobs), clip = FALSE)
    (title_panel / combined) + patchwork::plot_layout(
        heights = c(if (has_subtitle) 0.105 else 0.08, if (has_subtitle) 0.895 else 0.92)
    )
}
