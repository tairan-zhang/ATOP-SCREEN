# GseaVis calls are isolated from ATOP validation and styling.

resolve_gsea_pathway_colours_atop <- function(pathway_ids, pathway_colours = NULL) {
    colours <- as.character(pathway_colours %||% character())
    colours <- colours[grepl("^#[0-9A-Fa-f]{6}$", colours)]
    if (length(colours) == 0L) {
        colours <- grDevices::colorRampPalette(c("#C02635", "#EFC8CC", "#BFD7E6", "#00609D"))(length(pathway_ids))
    }
    rep(colours, length.out = length(pathway_ids))
}

generate_gsea_running_raw_atop <- function(gsea_object,
                                            pathway_ids,
                                            pathway_colours = NULL,
                                            highlighted_genes = NULL,
                                            highlight_colours = c("#00609D", "#C02635"),
                                            term_width = 20,
                                            legend_position = c(0.85, 0.8),
                                            subplot_type = 3) {
    if (!requireNamespace("GseaVis", quietly = TRUE)) {
        stop("Package 'GseaVis' is required for GSEA enrichment curves.", call. = FALSE)
    }

    # The original ATOP client deliberately uses GseaVis defaults for one pathway.
    if (length(pathway_ids) == 1L) {
        return(GseaVis::gseaNb(
            object = gsea_object,
            geneSetID = pathway_ids[[1]],
            subPlot = subplot_type,
            addPval = FALSE,
            addGene = highlighted_genes,
            rankCol = c("#00609D", "white", "#C02635"),
            newHtCol = c("#00609D", "white", "#C02635"),
            htCol = highlight_colours
        ))
    }

    final_colours <- resolve_gsea_pathway_colours_atop(pathway_ids, pathway_colours)
    plot_object <- GseaVis::gseaNb(
        object = gsea_object,
        termWidth = term_width,
        legend.position = legend_position,
        geneSetID = pathway_ids,
        curveCol = final_colours,
        subPlot = subplot_type,
        addPval = FALSE,
        addGene = highlighted_genes,
        rankCol = c("#00609D", "white", "#C02635"),
            newHtCol = c("#00609D", "white", "#C02635"),
            htCol = highlight_colours,
        rmPrefix = FALSE
    )

    pathway_rows <- gsea_pathway_rows_atop(gsea_object, pathway_ids)
    descriptions <- if ("Description" %in% names(pathway_rows)) {
        as.character(pathway_rows$Description)
    } else {
        pathway_ids
    }
    descriptions[is.na(descriptions)] <- pathway_ids[is.na(descriptions)]
    pathway_scale <- ggplot2::scale_color_manual(
        values = stats::setNames(final_colours, descriptions),
        breaks = descriptions,
        labels = stats::setNames(format_gsea_pathway_name_atop(descriptions, width = term_width), descriptions),
        name = "Pathway"
    )
    if (inherits(plot_object, "patchwork")) {
        plot_object[[1]] <- suppressMessages(plot_object[[1]] + pathway_scale)
    } else {
        plot_object <- suppressMessages(plot_object + pathway_scale)
    }
    plot_object
}
