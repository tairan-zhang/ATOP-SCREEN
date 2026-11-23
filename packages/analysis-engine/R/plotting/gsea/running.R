# Public orchestration for ATOP GSEA enrichment curves.

generate_gsea_enrichment_plot_atop <- function(gsea_object,
                                                pathway_ids,
                                                pathway_colours = NULL,
                                                highlighted_genes = NULL,
                                                highlight_colours = c("#00609D", "#C02635"),
                                                term_width = 20,
                                                legend_position = c(0.85, 0.8),
                                                subplot_type = 3,
                                                base_font_size = 11,
                                                font_family = "Helvetica",
                                                direction = "positive",
                                                show_nes = TRUE,
                                                show_pvalue = TRUE,
                                                pvalue_type = "p.adjust",
                                                show_rank_explanation = TRUE,
                                                rank_positive_label = "Top-ranked genes (positively correlated)",
                                                rank_negative_label = "Bottom-ranked genes (negatively correlated)",
                                                plot_title = NULL) {
    if (is.null(gsea_object) || !inherits(gsea_object, "gseaResult")) {
        stop("Input is not a valid gseaResult S4 object.", call. = FALSE)
    }
    pathway_ids <- unique(as.character(pathway_ids))
    pathway_ids <- pathway_ids[!is.na(pathway_ids) & nzchar(pathway_ids)]
    if (length(pathway_ids) == 0L) stop("Select at least one pathway.", call. = FALSE)

    available_pathways <- names(gsea_object@geneSets)
    missing_pathways <- setdiff(pathway_ids, available_pathways)
    if (length(missing_pathways) > 0L) {
        stop(
            "Selected pathway(s) were not found in the GSEA object: ",
            paste(missing_pathways, collapse = ", "),
            call. = FALSE
        )
    }

    highlighted_genes <- resolve_gsea_highlighted_genes_atop(
        gsea_object,
        pathway_ids,
        highlighted_genes
    )
    normalized <- normalize_gsea_pathway_labels_atop(gsea_object)
    mapped_pathway_ids <- unname(normalized$id_map[pathway_ids])
    mapped_pathway_ids[is.na(mapped_pathway_ids)] <- pathway_ids[is.na(mapped_pathway_ids)]
    subplot_type <- max(1L, min(3L, as.integer(subplot_type)))
    direction <- if (tolower(as.character(direction[[1]])) == "negative") "negative" else "positive"
    pvalue_type <- if (identical(as.character(pvalue_type[[1]]), "pvalue")) "pvalue" else "p.adjust"

    raw_plot <- tryCatch(
        generate_gsea_running_raw_atop(
            gsea_object = normalized$object,
            pathway_ids = mapped_pathway_ids,
            pathway_colours = pathway_colours,
            highlighted_genes = highlighted_genes,
            highlight_colours = highlight_colours,
            term_width = term_width,
            legend_position = legend_position,
            subplot_type = subplot_type
        ),
        error = function(error) {
            stop("Failed to generate the GseaVis running plot: ", conditionMessage(error), call. = FALSE)
        }
    )

    style_gsea_running_plot_atop(
        plot_object = raw_plot,
        plot_title = plot_title,
        direction = direction,
        pathway_rows = gsea_pathway_rows_atop(normalized$object, mapped_pathway_ids),
        subplot_type = subplot_type,
        pathway_count = length(mapped_pathway_ids),
        legend_position = legend_position,
        base_font_size = base_font_size,
        font_family = font_family,
        show_nes = isTRUE(show_nes),
        show_pvalue = isTRUE(show_pvalue),
        pvalue_type = pvalue_type,
        ranked_gene_list = normalized$object@geneList,
        show_rank_explanation = isTRUE(show_rank_explanation),
        rank_positive_label = as.character(rank_positive_label[[1]]),
        rank_negative_label = as.character(rank_negative_label[[1]])
    )
}
