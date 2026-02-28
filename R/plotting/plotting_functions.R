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


# plotting_functions.R
# Plotting Related Functions Collection

# --- GSEA Lollipop Plot Function ---
generate_gsea_lollipop_plot_shiny <- function(
    gsea_data_input_df, # Data frame from GSEA results (e.g., gsea_results_positive_df)
    analysis_type_label, # e.g. "Positive" or "Negative" for titling
    num_pathways_to_plot,
    plot_title_main,
    bar_fill_low_hex = "#E0F3F8", # Light blue / cyan
    bar_fill_high_hex = "#045275", # Dark blue
    circle_fill_low_hex = "#FFF2CC", # Light yellow / gold
    circle_fill_high_hex = "#D66000", # Dark orange / brick
    base_font_size = 11,
    bar_width_ratio = 0.1) {
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

    # Filter and sort pathways based on analysis_type_label (NES sign)
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
        # Default to using all, but this should be driven by selection usually
        gsea_data_filtered <- gsea_data_input_df %>%
            arrange(desc(abs(NES)), p.adjust) %>%
            head(num_pathways_to_plot)
    }

    if (is.null(gsea_data_filtered) || nrow(gsea_data_filtered) == 0) {
        stop(paste("No pathways for", analysis_type_label, "after filtering (NES sign & top N). No plot."))
    }

    # Prepare data for plotting
    gsea_data_filtered <- gsea_data_filtered %>%
        mutate(
            Description = factor(Description, levels = rev(unique(.$Description))), # Important for y-axis order
            negLog10PAdjust = -log10(p.adjust + .Machine$double.xmin),
            negLog10PValue = -log10(pvalue + .Machine$double.xmin), # Added for circle color
            absNES = abs(NES),
            DescriptionWrapped = stringr::str_wrap(as.character(Description), width = 50) # Fixed width to 50
        ) %>%
        mutate(DescriptionWrapped = factor(DescriptionWrapped, levels = stringr::str_wrap(levels(Description), width = 50))) # Fixed width to 50

    plot_subtitle <- paste("Top", nrow(gsea_data_filtered), analysis_type_label, "Enriched Pathways (by |NES|)")

    # Main plot object construction
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
            axis.text.y = element_text(size = rel(0.95), hjust = 1), # Relative to base_font_size
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

    return(plot_obj)
}

# --- sgRNA Paired Plot Function ---
generate_sgrna_paired_plot <- function(
    normalized_counts_df,
    target_gene_id,
    sgrna_id_col, # Name of the sgRNA ID column
    gene_col, # Name of the Gene ID column
    condition_1_col,
    condition_2_col,
    condition_1_label = NULL,
    condition_2_label = NULL,
    color_condition_1 = "black",
    color_condition_2 = "red",
    connecting_line_color = "grey50",
    connecting_line_size = 0.5,
    connecting_line_alpha = 0.7,
    point_size = 2.5,
    point_alpha = 0.8,
    plot_title = NULL,
    y_axis_label = "sgRNA Normalized Read Count",
    base_font_size = 12) {
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

    # Ensure columns for plotting are numeric
    gene_sgrna_data_subset[[condition_1_col]] <- as.numeric(gene_sgrna_data_subset[[condition_1_col]])
    gene_sgrna_data_subset[[condition_2_col]] <- as.numeric(gene_sgrna_data_subset[[condition_2_col]])

    # Remove rows with NAs in EITHER of the plotting columns for this gene
    gene_sgrna_data_subset <- gene_sgrna_data_subset[
        complete.cases(gene_sgrna_data_subset[[condition_1_col]]) &
            complete.cases(gene_sgrna_data_subset[[condition_2_col]]),
    ]

    if (nrow(gene_sgrna_data_subset) == 0) {
        stop(paste("No valid (non-NA in both conditions) sgRNA data for plotting gene:", target_gene_id, "after NA removal."))
    }

    # Use sgrna_id_col directly for unique IDs. If it's not unique per gene, this might be an issue.
    # The original script used make.unique, implying potential duplicates. For CRISPR screens, sgRNA IDs within a gene context should be unique.
    long_df_for_plot <- gene_sgrna_data_subset %>%
        rename(unique_sgrna_identifier_for_plot = !!sym(sgrna_id_col)) %>% # Prefer direct use if sgRNA IDs are truly unique
        pivot_longer(
            cols = c(all_of(condition_1_col), all_of(condition_2_col)),
            names_to = "Group",
            values_to = "Value"
        )

    # Explicitly set the factor levels for 'Group' to match the input order
    # This ensures condition_1_col is the first on the x-axis, condition_2_col is second.
    long_df_for_plot$Group <- factor(long_df_for_plot$Group, levels = c(condition_1_col, condition_2_col))

    actual_label_1 <- condition_1_label %||% condition_1_col
    actual_label_2 <- condition_2_label %||% condition_2_col

    # Define colors and labels for the plot aesthetics
    # The order in scale_color_manual and scale_x_discrete depends on the factor levels of 'Group'
    # which pivot_longer creates based on the order of columns provided to 'cols'.
    # So, condition_1_col will be the first level, condition_2_col the second.
    group_colors_map <- setNames(c(color_condition_1, color_condition_2), c(condition_1_col, condition_2_col))
    group_labels_map <- setNames(c(actual_label_1, actual_label_2), c(condition_1_col, condition_2_col))

    plot_final_title <- plot_title %||% paste("sgRNA Paired Plot for Gene:", target_gene_id)

    p <- ggplot(long_df_for_plot, aes(x = Group, y = Value, color = Group)) +
        geom_line(aes(group = unique_sgrna_identifier_for_plot),
            color = connecting_line_color,
            linewidth = connecting_line_size,
            alpha = connecting_line_alpha
        ) +
        geom_point(size = point_size, alpha = point_alpha) +
        scale_color_manual(values = group_colors_map, labels = group_labels_map, name = "Condition") + # Added legend name
        scale_x_discrete(labels = group_labels_map) + # Ensures x-axis ticks match condition labels
        labs(
            y = y_axis_label,
            x = NULL, # Typically no x-axis label for paired plots like this
            title = plot_final_title
        ) +
        ggpubr::theme_pubr(base_size = base_font_size) + # Using theme_pubr for consistency
        theme(
            plot.title = element_text(hjust = 0.5, size = rel(1.1)), # Relative sizing
            legend.position = "top" # Or "bottom", "none"
            # axis.text.x = element_text(angle = 45, hjust = 1) # If labels are long
        )

    return(p)
}

# --- Internal Function for Single GSEA Pathway Plot (based on GSEA.R) ---
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
            # Call gseaNb with parameters from the unification logic.
            # IMPORTANT: curveCol is INTENTIONALLY OMITTED here to avoid transparency issues in single plot.
            # We rely on GseaVis default colors for the single pathway case.
            p <- GseaVis::gseaNb(
                object = gsea_s4_object,
                geneSetID = pathway_id,
                subPlot = subplot_type,
                addPval = add_pval, # Restore user choice for default GseaVis labeling
                # pvalX and pvalY REMOVED to use GseaVis defaults
                addGene = highlighted_genes,
                htCol = highlight_colors
            )

            p <- p + ggplot2::theme(
                text = ggplot2::element_text(size = base_font_size),
                axis.text = ggplot2::element_text(size = base_font_size * 0.9),
                axis.title = ggplot2::element_text(size = base_font_size * 1.1),
                plot.title = ggplot2::element_text(size = base_font_size * 1.2),
                legend.text = ggplot2::element_text(size = base_font_size * 0.9),
                legend.title = ggplot2::element_text(size = base_font_size)
            )
            return(p)
        },
        error = function(e) {
            stop(paste("Failed to generate GSEA plot for pathway '", pathway_id, "': ", e$message))
        }
    )
}

# --- Unified GSEA Visualization Function ---
generate_gsea_multi_pathway_plot <- function(
    gsea_s4_object, # GSEA result object
    selected_pathway_ids, # List of selected pathway IDs
    pathway_colors = NULL, # Pathway colors, default used if NULL
    highlighted_genes = NULL, # List of genes to highlight
    highlight_colors = c("#0E6DB3", "#BB1E38"), # Fixed highlight colors: Blue and Red
    term_width = 20, # Pathway name display width
    legend_position = c(0.85, 0.8), # Legend position
    subplot_type = 2, # Subplot type
    add_pval = TRUE, # Whether to add p-value
    pval_x = 0.02, # p-value x position
    pval_y = 0.04, # p-value y position
    base_font_size = 10 # Base font size
    ) {
    # Basic input validation
    if (is.null(gsea_s4_object) || !inherits(gsea_s4_object, "gseaResult")) {
        stop("Input is not a valid gseaResult S4 object.")
    }
    if (is.null(selected_pathway_ids) || length(selected_pathway_ids) == 0) {
        stop("Please select at least one pathway.")
    }
    if (!requireNamespace("GseaVis", quietly = TRUE)) {
        stop("Package 'GseaVis' is required. Please install it.")
    }

    # --- Logic Dispatch based on pathway count ---
    if (length(selected_pathway_ids) == 1) {
        # --- Single Pathway Logic ---
        # Delegate to internal function.
        # Note: We do NOT pass pathway_colors (curveCol) to avoid transparency issues in single plot.
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
        # --- Multi-Pathway Logic ---
        final_colors <- pathway_colors
        if (is.null(final_colors)) {
            final_colors <- if (length(selected_pathway_ids) <= 8) {
                RColorBrewer::brewer.pal(max(3, length(selected_pathway_ids)), "Set2")[1:length(selected_pathway_ids)]
            } else {
                rainbow(length(selected_pathway_ids))
            }
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
                    curveCol = final_colors, # Pass colors for multi-pathway
                    subPlot = subplot_type,
                    addPval = add_pval,
                    addGene = highlighted_genes,
                    htCol = highlight_colors, # Use passed highlight colors and not hardcoded
                    pvalX = pval_x,
                    pvalY = pval_y
                )
                p <- p + ggplot2::theme(
                    text = ggplot2::element_text(size = base_font_size),
                    axis.text = ggplot2::element_text(size = base_font_size * 0.9),
                    axis.title = ggplot2::element_text(size = base_font_size * 1.1),
                    plot.title = ggplot2::element_text(size = base_font_size * 1.2),
                    legend.text = ggplot2::element_text(size = base_font_size * 0.9),
                    legend.title = ggplot2::element_text(size = base_font_size)
                )
                return(p)
            },
            error = function(e) {
                stop(paste("Failed to generate GSEA plot using GseaVis:", e$message))
            }
        )
    }
}

# --- Refactored Functions from 03_visualization_functions.R ---

# --- Volcano Plot Function ---
generate_volcano_plot <- function(data,
                                  p_value_cutoff = 0.05,
                                  score_cutoff = 0.2,
                                  label_genes = NULL,
                                  col_up = "#E31A1C",
                                  col_down = "#1F78B4",
                                  col_ns = "#ADB6B6",
                                  base_font_size = 14,
                                  plot_title = "Volcano Plot") {
    # Data validation
    required_cols <- c("Gene", "GeneNegativeScore", "P_negative", "GenePositiveScore", "P_positive")
    if (!all(required_cols %in% colnames(data))) {
        stop(paste("Data missing required columns:", paste(setdiff(required_cols, colnames(data)), collapse = ", ")))
    }

    # Data processing
    data_neg <- data %>%
        filter(GeneNegativeScore < 0) %>%
        select(Gene, Score = GeneNegativeScore, P = P_negative)

    data_pos <- data %>%
        filter(GenePositiveScore > 0) %>%
        select(Gene, Score = GenePositiveScore, P = P_positive)

    plot_data <- bind_rows(data_neg, data_pos)

    # Define groups
    plot_data$group <- "NS"
    plot_data$group[plot_data$Score > score_cutoff & plot_data$P < p_value_cutoff] <- "Up"
    plot_data$group[plot_data$Score < -score_cutoff & plot_data$P < p_value_cutoff] <- "Down"

    # Handle P = 0
    min_nonzero_p <- min(plot_data$P[plot_data$P > 0], na.rm = TRUE)
    if (any(plot_data$P == 0)) {
        plot_data$P[plot_data$P == 0] <- min_nonzero_p / 10
    }

    # Calculate Y-axis limit
    max_log_p <- max(-log10(plot_data$P), na.rm = TRUE)
    y_limit <- max_log_p * 1.2

    # Prepare label data
    if (!is.null(label_genes) && length(label_genes) > 0) {
        label_data <- plot_data[plot_data$Gene %in% label_genes, ]
    } else {
        label_data <- plot_data[FALSE, ]
    }

    # Plotting
    p <- ggplot(data = plot_data, aes(x = Score, y = -log10(P), color = group)) +
        geom_point(alpha = 1, size = 1.2) +
        scale_color_manual(values = c("Down" = col_down, "NS" = col_ns, "Up" = col_up), name = "Group") +
        scale_y_continuous(expand = expansion(add = c(0, 0)), limits = c(0, y_limit)) +
        geom_hline(yintercept = -log10(p_value_cutoff), lty = 4, lwd = 0.6, alpha = 0.8, color = "black", show.legend = FALSE) +
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
            x = max(plot_data$Score, na.rm = T) * 0.8, y = y_limit * 0.9,
            label = paste0("Up = ", sum(plot_data$group == "Up")), size = base_font_size / 3, color = col_up
        ) +
        annotate("text",
            x = min(plot_data$Score, na.rm = T) * 0.8, y = y_limit * 0.9,
            label = paste0("Down = ", sum(plot_data$group == "Down")), size = base_font_size / 3, color = col_down
        )

    # Vertical lines
    if (score_cutoff > 0) {
        p <- p + geom_vline(xintercept = c(-score_cutoff, score_cutoff), lty = 4, lwd = 0.6, alpha = 0.8, color = "black", show.legend = FALSE)
    }

    # Labels
    if (nrow(label_data) > 0) {
        p <- p + geom_text_repel(
            data = label_data, aes(x = Score, y = -log10(P), label = Gene, color = group),
            seed = 123, size = base_font_size / 3, min.segment.length = 0, show.legend = F,
            box.padding = 0.5,
            max.overlaps = Inf, segment.linetype = 1, segment.alpha = 0.8,
            force = 2
        )
    }

    return(p)
}

# --- Gene Score Ranking Plot Function ---
generate_ranking_plot <- function(data,
                                  type = "positive", # "positive" or "negative"
                                  top_n = 10,
                                  gradient_low = "#fee0d2",
                                  gradient_high = "#a50f15",
                                  base_font_size = 12,
                                  plot_title = NULL) {
    # Data validation
    if (type == "positive") {
        required_cols <- c("Gene", "GenePositiveRank", "GenePositiveScore")
    } else {
        required_cols <- c("Gene", "GeneNegativeRank", "GeneNegativeScore")
    }

    if (!all(required_cols %in% colnames(data))) {
        stop(paste("Data missing required columns for", type, "ranking:", paste(setdiff(required_cols, colnames(data)), collapse = ", ")))
    }

    # Filter out microRNAs if needed (keeping original logic)
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
            aes(x = score + (max(plot_data$score) * 0.02), y = gene, label = rank),
            color = "black", size = base_font_size / 3.5, hjust = 0
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

    return(p)
}
