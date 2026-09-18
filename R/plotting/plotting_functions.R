# SPDX-License-Identifier: GPL-3.0-or-later
# plotting_functions.R
# Part of ATOP-SCREEN
#
# Copyright (C) 2026 Tairan Zhang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

atop_palette <- c("#104e8b", "#376b9e", "#5f89b1", "#afc3d8", "#d7e1eb",
                  "#f2dada", "#e5b5b5", "#d89090", "#b22222")

atop_pathway_colors <- function(n) {
    rep(atop_palette[c(1, 9, 3, 8, 2, 7, 4, 6, 5)], length.out = n)
}

publication_device <- function(format) {
    if (!identical(format, "pdf")) return(format)
    if (identical(Sys.info()[["sysname"]], "Darwin")) {
        function(filename, width, height, ...) {
            grDevices::quartz(file = filename, type = "pdf", width = width, height = height, ...)
        }
    } else grDevices::cairo_pdf
}

style_publication_plot <- function(plot, base_size = 8) {
    if (inherits(plot, "patchwork")) {
        for (i in seq_len(length(plot))) {
            plot[[i]] <- style_publication_plot(plot[[i]], base_size)
        }
        return(plot)
    }
    for (i in seq_along(plot$layers)) {
        layer <- plot$layers[[i]]
        if (any(grepl("Text|Label", class(layer$geom)))) {
            layer$mapping$colour <- NULL
            layer$mapping$family <- NULL
            layer$aes_params$colour <- "black"
            layer$aes_params$family <- "Arial"
            layer$aes_params$fontface <- "plain"
            plot$layers[[i]] <- layer
        }
    }
    style <- ggplot2::theme(
        text = ggplot2::element_text(family = "Arial", size = base_size, colour = "black"),
        axis.text = ggplot2::element_text(size = base_size, colour = "black"),
        axis.title = ggplot2::element_text(size = base_size, colour = "black"),
        plot.title = ggplot2::element_text(size = base_size, colour = "black", face = "bold"),
        plot.subtitle = ggplot2::element_text(size = base_size, colour = "black"),
        plot.caption = ggplot2::element_text(size = base_size, colour = "black"),
        strip.text = ggplot2::element_text(size = base_size, colour = "black"),
        legend.text = ggplot2::element_text(size = base_size, colour = "black"),
        legend.title = ggplot2::element_text(size = base_size, colour = "black"),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA))
    for (name in names(plot$theme)) {
        if (inherits(plot$theme[[name]], "element_text")) {
            plot$theme[[name]]$family <- "Arial"
            plot$theme[[name]]$colour <- "black"
        }
    }
    plot + style
}

add_gsea_statistics <- function(plot, object, ids, x, y, base_size) {
    results <- object@result[match(ids, object@result$ID), , drop = FALSE]
    format_p <- function(p) ifelse(p < 0.001, "< 0.001", formatC(p, format = "g", digits = 3))
    labels <- paste0("NES: ", formatC(results$NES, format = "f", digits = 2),
                     "\nP value: ", format_p(results$pvalue),
                     "\nAdjusted P value: ", format_p(results$p.adjust))
    if (length(ids) > 1) labels <- paste(results$Description, labels, sep = "\n")
    panel <- if (inherits(plot, "patchwork")) plot[[1]] else plot
    ranges <- ggplot2::ggplot_build(panel)$layout$panel_params[[1]]
    panel <- panel + ggplot2::annotate("label",
        x = ranges$x.range[1] + x * diff(ranges$x.range),
        y = ranges$y.range[1] + y * diff(ranges$y.range),
        label = paste(labels, collapse = "\n\n"),
        hjust = if (x > 0.5) 1 else 0, vjust = if (y > 0.5) 1 else 0,
        family = "Arial", colour = "black", fill = "white", linewidth = 0,
        size = base_size / ggplot2::.pt)
    if (inherits(plot, "patchwork")) {
        plot[[1]] <- panel
        plot
    } else panel
}

generate_gsea_lollipop_plot_shiny <- function(
    gsea_data_input_df,
    analysis_type_label,
    num_pathways_to_plot,
    plot_title_main,
    bar_fill_low_hex = "#d7e1eb",
    bar_fill_high_hex = "#104e8b",
    circle_fill_low_hex = "#f2dada",
    circle_fill_high_hex = "#b22222",
    base_font_size = 8,
    bar_width_ratio = 0.3) {
    if (is.null(gsea_data_input_df) || nrow(gsea_data_input_df) == 0) {
        stop(paste("GSEA result data for", analysis_type_label, "is empty or NULL. Cannot generate plot."))
    }
    if (!"Description" %in% names(gsea_data_input_df) && "ID" %in% names(gsea_data_input_df)) {
        gsea_data_input_df <- gsea_data_input_df %>% rename(Description = ID)
    } else if (!"Description" %in% names(gsea_data_input_df)) {
        stop("GSEA data must have 'Description' or 'ID' column.")
    }
    required_gsea_cols <- c("NES", "pvalue", "p.adjust", "setSize")
    missing_gsea_cols <- setdiff(required_gsea_cols, names(gsea_data_input_df))
    if (length(missing_gsea_cols) > 0) {
        stop(paste("GSEA data missing required columns:", paste(missing_gsea_cols, collapse = ", ")))
    }

    gsea_data_filtered <- NULL
    if (grepl("Positive", analysis_type_label, ignore.case = TRUE)) {
        gsea_data_filtered <- gsea_data_input_df %>%
            filter(NES > 0) %>%
            arrange(desc(abs(NES)), p.adjust) %>%
            head(num_pathways_to_plot)
    } else if (grepl("Negative", analysis_type_label, ignore.case = TRUE)) {
        gsea_data_filtered <- gsea_data_input_df %>%
            filter(NES < 0) %>%
            arrange(desc(abs(NES)), p.adjust) %>%
            head(num_pathways_to_plot)
    } else {
        gsea_data_filtered <- gsea_data_input_df %>%
            arrange(desc(abs(NES)), p.adjust) %>%
            head(num_pathways_to_plot)
    }

    if (is.null(gsea_data_filtered) || nrow(gsea_data_filtered) == 0) {
        stop(paste("No pathways for", analysis_type_label, "after filtering (NES sign & top N). No plot."))
    }

    gsea_data_filtered <- gsea_data_filtered %>%
        mutate(
            Description = factor(Description, levels = rev(unique(.$Description))),
            negLog10PAdjust = -log10(p.adjust + .Machine$double.xmin),
            negLog10PValue = -log10(pvalue + .Machine$double.xmin),
            absNES = abs(NES),
            DescriptionWrapped = stringr::str_wrap(as.character(Description), width = 50)
        ) %>%
        mutate(DescriptionWrapped = factor(DescriptionWrapped, levels = stringr::str_wrap(levels(Description), width = 50)))

    plot_subtitle <- paste("Top", nrow(gsea_data_filtered), analysis_type_label, "Enriched Pathways (by |NES|)")

    plot_obj <- ggplot(gsea_data_filtered, aes(y = DescriptionWrapped)) +
        geom_col(aes(x = NES, fill = absNES), width = bar_width_ratio, alpha = 0.7) +
        scale_fill_gradient(name = "|NES|\n(Bar Color)", low = bar_fill_low_hex, high = bar_fill_high_hex) +
        ggnewscale::new_scale_fill() +
        geom_point(aes(x = NES, size = setSize, fill = negLog10PValue), shape = 21, color = "black", stroke = 0.6) +
        scale_size_continuous(name = "Set Size\n(Circle Size)", range = c(3, 10)) +
        scale_fill_gradient(name = "-log10(p-value)\n(Circle Fill)", low = circle_fill_low_hex, high = circle_fill_high_hex) +
        labs(
            x = "Normalized Enrichment Score (NES)", y = NULL,
            title = plot_title_main,
            subtitle = plot_subtitle
        ) +
        theme_minimal(base_size = base_font_size) +
        theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_line(linetype = "dashed", color = "gray80"),
            axis.text.y = element_text(size = rel(0.95), hjust = 1),
            axis.text.x = element_text(size = rel(0.95)),
            axis.title.x = element_text(size = rel(1.05), face = "bold", margin = margin(t = 10)),
            plot.title = element_text(hjust = 0.5, size = rel(1.3), face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = rel(1.1), margin = margin(b = 15)),
            legend.position = "right",
            legend.title = element_text(size = rel(0.85), face = "bold"),
            legend.text = element_text(size = rel(0.8)),
            legend.key.size = unit(0.5, "lines"),
            legend.spacing.y = unit(0.1, "cm")
        ) +
        geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5)

    return(style_publication_plot(plot_obj, base_font_size))
}

generate_sgrna_paired_plot <- function(
    normalized_counts_df,
    target_gene_id,
    sgrna_id_col,
    gene_col,
    condition_1_col,
    condition_2_col,
    condition_1_label = NULL,
    condition_2_label = NULL,
    color_condition_1 = "#104e8b",
    color_condition_2 = "#b22222",
    connecting_line_color = "#afc3d8",
    connecting_line_size = 0.5,
    connecting_line_alpha = 0.7,
    point_size = 2.5,
    point_alpha = 0.8,
    plot_title = NULL,
    y_axis_label = "sgRNA Normalized Read Count",
    base_font_size = 8) {
    if (!target_gene_id %in% normalized_counts_df[[gene_col]]) {
        stop(paste("Target gene '", target_gene_id, "' not found in the provided data under column '", gene_col, "'.", sep = ""))
    }
    if (!sgrna_id_col %in% names(normalized_counts_df)) stop(paste("sgRNA ID column '", sgrna_id_col, "' not found."))
    if (!condition_1_col %in% names(normalized_counts_df)) stop(paste("Condition 1 column '", condition_1_col, "' not found."))
    if (!condition_2_col %in% names(normalized_counts_df)) stop(paste("Condition 2 column '", condition_2_col, "' not found."))

    gene_sgrna_data_subset <- normalized_counts_df %>%
        filter(!!sym(gene_col) == target_gene_id) %>%
        select(all_of(sgrna_id_col), all_of(condition_1_col), all_of(condition_2_col))

    if (nrow(gene_sgrna_data_subset) == 0) {
        stop(paste("No sgRNA data found for gene:", target_gene_id))
    }

    gene_sgrna_data_subset[[condition_1_col]] <- as.numeric(gene_sgrna_data_subset[[condition_1_col]])
    gene_sgrna_data_subset[[condition_2_col]] <- as.numeric(gene_sgrna_data_subset[[condition_2_col]])

    gene_sgrna_data_subset <- gene_sgrna_data_subset[
        complete.cases(gene_sgrna_data_subset[[condition_1_col]]) &
            complete.cases(gene_sgrna_data_subset[[condition_2_col]]),
    ]

    if (nrow(gene_sgrna_data_subset) == 0) {
        stop(paste("No valid (non-NA in both conditions) sgRNA data for plotting gene:", target_gene_id, "after NA removal."))
    }

    long_df_for_plot <- gene_sgrna_data_subset %>%
        rename(unique_sgrna_identifier_for_plot = !!sym(sgrna_id_col)) %>%
        pivot_longer(
            cols = c(all_of(condition_1_col), all_of(condition_2_col)),
            names_to = "Group",
            values_to = "Value"
        )

    long_df_for_plot$Group <- factor(long_df_for_plot$Group, levels = c(condition_1_col, condition_2_col))

    actual_label_1 <- condition_1_label %||% condition_1_col
    actual_label_2 <- condition_2_label %||% condition_2_col

    group_colors_map <- setNames(c(color_condition_1, color_condition_2), c(condition_1_col, condition_2_col))
    group_labels_map <- setNames(c(actual_label_1, actual_label_2), c(condition_1_col, condition_2_col))

    plot_final_title <- plot_title %||% target_gene_id

    p <- ggplot(long_df_for_plot, aes(x = Group, y = Value, color = Group)) +
        geom_line(aes(group = unique_sgrna_identifier_for_plot),
            color = connecting_line_color,
            linewidth = connecting_line_size,
            alpha = connecting_line_alpha
        ) +
        geom_point(size = point_size, alpha = point_alpha) +
        scale_color_manual(values = group_colors_map, labels = group_labels_map, name = "Condition") +
        scale_x_discrete(labels = group_labels_map) +
        labs(
            y = y_axis_label,
            x = NULL,
            title = plot_final_title
        ) +
        ggpubr::theme_pubr(base_size = base_font_size) +
        theme(
            plot.title = element_text(hjust = 0.5, size = rel(1.1)),
            legend.position = "none"

        )

    return(style_publication_plot(p, base_font_size))
}

.generate_gsea_single_from_script_logic <- function(
    gsea_s4_object,
    pathway_id,
    highlighted_genes,
    highlight_colors,
    base_font_size,
    subplot_type,
    add_pval,
    pval_x,
    pval_y) {
    tryCatch(
        {
            p <- GseaVis::gseaNb(
                object = gsea_s4_object,
                geneSetID = pathway_id,
                curveCol = atop_palette[1],
                subPlot = subplot_type,
                addPval = FALSE,
                pvalX = pval_x,
                pvalY = pval_y,
                addGene = highlighted_genes,
                htCol = highlight_colors,
                rankCol = atop_palette[c(1, 5, 9)], base_size = base_font_size,
                segCol = atop_palette[9], geneCol = "black",
                geneSize = base_font_size / ggplot2::.pt, pvalSize = base_font_size / ggplot2::.pt
            )

            if (add_pval) p <- add_gsea_statistics(p, gsea_s4_object, pathway_id, pval_x, pval_y, base_font_size)
            p <- p + ggplot2::theme(
                text = ggplot2::element_text(size = base_font_size),
                axis.text = ggplot2::element_text(size = base_font_size),
                axis.title = ggplot2::element_text(size = base_font_size),
                plot.title = ggplot2::element_text(size = base_font_size),
                legend.text = ggplot2::element_text(size = base_font_size),
                legend.title = ggplot2::element_text(size = base_font_size)
            )
            return(style_publication_plot(p, base_font_size))
        },
        error = function(e) {
            stop(paste("Failed to generate GSEA plot for pathway '", pathway_id, "': ", e$message))
        }
    )
}

generate_gsea_multi_pathway_plot <- function(
    gsea_s4_object,
    selected_pathway_ids,
    pathway_colors = NULL,
    highlighted_genes = NULL,
    highlight_colors = c("#104e8b", "#b22222"),
    term_width = 20,
    legend_position = c(0.85, 0.8),
    subplot_type = 2,
    add_pval = TRUE,
    pval_x = 0.02,
    pval_y = 0.04,
    base_font_size = 8
    ) {
    if (is.null(gsea_s4_object) || !inherits(gsea_s4_object, "gseaResult")) {
        stop("Input is not a valid gseaResult S4 object.")
    }
    if (is.null(selected_pathway_ids) || length(selected_pathway_ids) == 0) {
        stop("Please select at least one pathway.")
    }
    if (!requireNamespace("GseaVis", quietly = TRUE)) {
        stop("Package 'GseaVis' is required. Please install it.")
    }

    if (length(selected_pathway_ids) == 1) {
        return(.generate_gsea_single_from_script_logic(
            gsea_s4_object = gsea_s4_object,
            pathway_id = selected_pathway_ids[1],
            highlighted_genes = highlighted_genes,
            highlight_colors = highlight_colors,
            base_font_size = base_font_size,
            subplot_type = subplot_type,
            add_pval = add_pval,
            pval_x = pval_x,
            pval_y = pval_y
        ))
    } else {
        final_colors <- pathway_colors
        if (is.null(final_colors)) {
            final_colors <- atop_pathway_colors(length(selected_pathway_ids))
        } else if (length(final_colors) < length(selected_pathway_ids)) {
            final_colors <- rep(final_colors, length.out = length(selected_pathway_ids))
        }

        tryCatch(
            {
                p <- GseaVis::gseaNb(
                    object = gsea_s4_object,
                    termWidth = term_width,
                    legend.position = legend_position,
                    geneSetID = selected_pathway_ids,
                    curveCol = final_colors,
                    subPlot = subplot_type,
                    addPval = FALSE,
                    addGene = highlighted_genes,
                    htCol = highlight_colors,
                rankCol = atop_palette[c(1, 5, 9)], base_size = base_font_size,
                segCol = atop_palette[9], geneCol = "black",
                geneSize = base_font_size / ggplot2::.pt, pvalSize = base_font_size / ggplot2::.pt,
                    pvalX = pval_x,
                    pvalY = pval_y
                )
                if (add_pval) p <- add_gsea_statistics(p, gsea_s4_object, selected_pathway_ids, pval_x, pval_y, base_font_size)
                p <- p + ggplot2::theme(
                    text = ggplot2::element_text(size = base_font_size),
                    axis.text = ggplot2::element_text(size = base_font_size),
                    axis.title = ggplot2::element_text(size = base_font_size),
                    plot.title = ggplot2::element_text(size = base_font_size),
                    legend.text = ggplot2::element_text(size = base_font_size),
                    legend.title = ggplot2::element_text(size = base_font_size)
                )
                return(style_publication_plot(p, base_font_size))
            },
            error = function(e) {
                stop(paste("Failed to generate GSEA plot using GseaVis:", e$message))
            }
        )
    }
}

generate_volcano_plot <- function(data,
                                  p_value_cutoff = 0.05,
                                  score_cutoff = 0.2,
                                  label_genes = NULL,
                                  col_up = "#b22222",
                                  col_down = "#104e8b",
                                  col_ns = "#d7e1eb",
                                  base_font_size = 8,
                                  plot_title = "Volcano Plot") {
    required_cols <- c("Gene", "GeneNegativeScore", "P_negative", "GenePositiveScore", "P_positive")
    if (!all(required_cols %in% colnames(data))) {
        stop(paste("Data missing required columns:", paste(setdiff(required_cols, colnames(data)), collapse = ", ")))
    }

    data_neg <- data %>%
        filter(GeneNegativeScore < 0) %>%
        select(Gene, Score = GeneNegativeScore, P = P_negative)

    data_pos <- data %>%
        filter(GenePositiveScore > 0) %>%
        select(Gene, Score = GenePositiveScore, P = P_positive)

    plot_data <- bind_rows(data_neg, data_pos)
    plot_data <- plot_data[is.finite(plot_data$Score) & is.finite(plot_data$P), ]
    if (!nrow(plot_data)) stop("No finite scores and P-values to plot.")
    if (any(plot_data$P < 0 | plot_data$P > 1)) stop("P-values must be between 0 and 1.")

    plot_data$group <- "NS"
    plot_data$group[plot_data$Score > score_cutoff & plot_data$P < p_value_cutoff] <- "Up"
    plot_data$group[plot_data$Score < -score_cutoff & plot_data$P < p_value_cutoff] <- "Down"

    positive_pvalues <- plot_data$P[plot_data$P > 0]
    p_floor <- if (length(positive_pvalues)) max(min(positive_pvalues) / 10, .Machine$double.xmin) else .Machine$double.xmin
    plot_data$P[plot_data$P == 0] <- p_floor
    y_limit <- max(1, -log10(plot_data$P)) * 1.2

    if (!is.null(label_genes) && length(label_genes) > 0) {
        label_data <- plot_data[plot_data$Gene %in% label_genes, ]
    } else {
        label_data <- plot_data[FALSE, ]
    }

    p <- ggplot(data = plot_data, aes(x = Score, y = -log10(P), color = group)) +
        geom_point(alpha = 1, size = 1.2) +
        scale_color_manual(values = c("Down" = col_down, "NS" = col_ns, "Up" = col_up), name = "Group") +
        scale_y_continuous(expand = expansion(add = c(0, 0)), limits = c(0, y_limit)) +
        geom_hline(yintercept = -log10(p_value_cutoff), lty = 4, linewidth = 0.6, alpha = 0.8, color = "black", show.legend = FALSE) +
        labs(x = "Gene Score", y = "-log10(P-value)", title = plot_title) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme_minimal(base_size = base_font_size) +
        theme(
            panel.border = element_rect(fill = NA, color = "black", linetype = "solid"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "top",
            axis.line = element_blank()
        ) +
        annotate("text",
            x = max(plot_data$Score, na.rm = TRUE) * 0.8, y = y_limit * 0.9,
            label = paste0("Up = ", sum(plot_data$group == "Up")), size = base_font_size / ggplot2::.pt, color = col_up
        ) +
        annotate("text",
            x = min(plot_data$Score, na.rm = TRUE) * 0.8, y = y_limit * 0.9,
            label = paste0("Down = ", sum(plot_data$group == "Down")), size = base_font_size / ggplot2::.pt, color = col_down
        )

    if (score_cutoff > 0) {
        p <- p + geom_vline(xintercept = c(-score_cutoff, score_cutoff), lty = 4, linewidth = 0.6, alpha = 0.8, color = "black", show.legend = FALSE)
    }

    if (nrow(label_data) > 0) {
        p <- p + geom_text_repel(
            data = label_data, aes(x = Score, y = -log10(P), label = Gene, color = group),
            seed = 123, size = base_font_size / ggplot2::.pt, min.segment.length = 0, show.legend = FALSE,
            box.padding = 0.5,
            max.overlaps = Inf, segment.linetype = 1, segment.alpha = 0.8,
            force = 2
        )
    }

    return(style_publication_plot(p, base_font_size))
}

generate_ranking_plot <- function(data,
                                  type = "positive",
                                  top_n = 10,
                                  gradient_low = "#f2dada",
                                  gradient_high = "#b22222",
                                  base_font_size = 8,
                                  plot_title = NULL) {
    if (type == "positive") {
        required_cols <- c("Gene", "GenePositiveRank", "GenePositiveScore")
    } else {
        required_cols <- c("Gene", "GeneNegativeRank", "GeneNegativeScore")
    }

    if (!all(required_cols %in% colnames(data))) {
        stop(paste("Data missing required columns for", type, "ranking:", paste(setdiff(required_cols, colnames(data)), collapse = ", ")))
    }

    data <- data %>% filter(!grepl("^hsa-mir-", Gene, ignore.case = TRUE))

    if (type == "positive") {
        plot_data <- data %>%
            arrange(GenePositiveRank) %>%
            slice_head(n = top_n) %>%
            mutate(rank = GenePositiveRank, gene = Gene, score = GenePositiveScore)

        default_title <- "Positive gene scores ranking"
    } else {
        plot_data <- data %>%
            arrange(GeneNegativeRank) %>%
            slice_head(n = top_n) %>%
            mutate(rank = GeneNegativeRank, gene = Gene, score = GeneNegativeScore)

        default_title <- "Negative gene scores ranking"
    }

    final_title <- if (is.null(plot_title)) default_title else plot_title

    plot_data$gene <- factor(plot_data$gene, levels = rev(plot_data$gene))

    p <- ggplot(plot_data, aes(x = score, y = gene, fill = score)) +
        geom_col(width = 0.8) +
        scale_fill_gradient(low = gradient_low, high = gradient_high, name = "Score") +
        geom_text(
            aes(x = score + (max(score) * 0.02), y = gene, label = rank),
            color = "black", size = base_font_size / ggplot2::.pt, hjust = 0
        ) +
        labs(title = final_title, subtitle = paste0("Top-", top_n, " only"), x = "Score", y = "Gene") +
        theme_minimal(base_size = base_font_size) +
        theme(
            plot.title = element_text(face = "bold"),
            legend.position = "bottom",
            panel.grid.major.y = element_line(color = "gray90", linewidth = 0.8),
            panel.grid.minor = element_blank()
        ) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.15)))

    if (type == "negative") {
        p <- p + scale_y_discrete(position = "right")
    }

    return(style_publication_plot(p, base_font_size))
}
