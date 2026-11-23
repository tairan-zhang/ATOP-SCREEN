visualization_defaults <- function() {
    list(
        plotType = "gene_volcano",
        direction = "positive",
        topN = 10L,
        pValueCutoff = 0.05,
        scoreCutoff = 0.2,
        targetGene = NULL,
        pathwayIds = character(),
        pathwayColors = character(),
        labelGenes = character(),
        knownGenes = character(),
        highlightedGenes = character(),
        colUp = "#C02635",
        colDown = "#00609D",
        colNs = "#CECECE",
        gradientLow = "#EFC8CC",
        gradientHigh = "#C02635",
        positiveColor = "#EFC8CC",
        negativeColor = "#BFD7E6",
        condition1Color = "#00609D",
        condition2Color = "#C02635",
        connectingLineColor = "#7F7F7F",
        title = NULL,
        widthCm = 18,
        heightCm = 14,
        fontSizePt = 8,
        pngDpi = 300,
        pathwayWrapWidth = 36L,
        gseaTermWidth = 20L,
        legendX = 0.85,
        legendY = 0.8,
        subplotType = 3L,
        showNes = TRUE,
        showPvalue = TRUE,
        pvalueType = "p.adjust",
        showRankExplanation = TRUE,
        rankPositiveLabel = "Top-ranked genes (positively correlated)",
        rankNegativeLabel = "Bottom-ranked genes (negatively correlated)",
        condition1Column = NULL,
        condition2Column = NULL
    )
}

visualization_dimension_defaults <- function(plot_type) {
    switch(
        as.character(plot_type[[1]] %||% "gene_volcano"),
        gene_volcano = list(width = 6.5, height = 6.5),
        gene_ranking = list(width = 11, height = 8.5),
        sgrna_paired = list(width = 6.5, height = 6.5),
        gsea_lollipop = list(width = 13, height = 11.5),
        gsea_running = list(width = 10.5, height = 7),
        list(width = 10.5, height = 7)
    )
}

normalize_visualization_color <- function(value, fallback) {
    candidate <- as.character(value[[1]] %||% fallback)
    valid <- tryCatch(
        {
            grDevices::col2rgb(candidate)
            TRUE
        },
        error = function(error) FALSE
    )
    if (valid) candidate else fallback
}

resolve_visualization_parameters <- function(value = NULL) {
    params <- visualization_defaults()
    value <- value %||% list()
    for (key in intersect(names(params), names(value))) {
        if (!is.null(value[[key]])) params[[key]] <- value[[key]]
    }
    params$plotType <- as.character(params$plotType[[1]] %||% "gene_volcano")
    dimension_defaults <- visualization_dimension_defaults(params$plotType)
    params$direction <- if (identical(as.character(params$direction[[1]]), "negative")) "negative" else "positive"
    params$topN <- max(1L, as.integer(params$topN[[1]] %||% 10L))
    params$pValueCutoff <- max(.Machine$double.xmin, min(1, as.numeric(params$pValueCutoff[[1]] %||% 0.05)))
    params$scoreCutoff <- max(0, as.numeric(params$scoreCutoff[[1]] %||% 0.2))
    params$widthCm <- max(3, as.numeric((value$widthCm %||% dimension_defaults$width)[[1]]))
    params$heightCm <- max(3, as.numeric((value$heightCm %||% dimension_defaults$height)[[1]]))
    params$fontSizePt <- max(4, min(36, as.numeric(params$fontSizePt[[1]] %||% 8)))
    params$pngDpi <- max(72L, min(1200L, as.integer(params$pngDpi[[1]] %||% 300L)))
    params$width <- params$widthCm / 2.54
    params$height <- params$heightCm / 2.54
    params$pathwayWrapWidth <- max(8L, as.integer(params$pathwayWrapWidth[[1]] %||% 36L))
    params$gseaTermWidth <- max(8L, as.integer(params$gseaTermWidth[[1]] %||% 20L))
    params$legendX <- max(0, min(1, as.numeric(params$legendX[[1]] %||% 0.85)))
    params$legendY <- max(0, min(1, as.numeric(params$legendY[[1]] %||% 0.8)))
    params$subplotType <- max(1L, min(3L, as.integer(params$subplotType[[1]] %||% 3L)))
    params$showNes <- isTRUE(params$showNes)
    params$showPvalue <- isTRUE(params$showPvalue)
    params$pvalueType <- if (identical(as.character(params$pvalueType[[1]]), "pvalue")) "pvalue" else "p.adjust"
    params$showRankExplanation <- isTRUE(params$showRankExplanation)
    params$rankPositiveLabel <- as.character(params$rankPositiveLabel[[1]] %||% "Top-ranked genes (positively correlated)")
    params$rankNegativeLabel <- as.character(params$rankNegativeLabel[[1]] %||% "Bottom-ranked genes (negatively correlated)")
    params$targetGene <- if (is.null(params$targetGene)) NULL else as.character(params$targetGene[[1]])
    params$pathwayIds <- as.character(unlist(params$pathwayIds %||% character(), use.names = FALSE))
    params$pathwayColors <- as.character(unlist(params$pathwayColors %||% character(), use.names = FALSE))
    if (length(params$pathwayColors) > 0L) {
        params$pathwayColors <- unname(vapply(
            params$pathwayColors,
            normalize_visualization_color,
            character(1),
            fallback = "#C02635"
        ))
    }
    params$labelGenes <- as.character(unlist(params$labelGenes %||% character(), use.names = FALSE))
    params$knownGenes <- as.character(unlist(params$knownGenes %||% character(), use.names = FALSE))
    params$highlightedGenes <- as.character(unlist(params$highlightedGenes %||% character(), use.names = FALSE))
    color_defaults <- visualization_defaults()
    if (params$direction == "negative") {
        color_defaults$gradientLow <- "#BFD7E6"
        color_defaults$gradientHigh <- "#00609D"
        if (is.null(value$gradientLow)) params$gradientLow <- color_defaults$gradientLow
        if (is.null(value$gradientHigh)) params$gradientHigh <- color_defaults$gradientHigh
    }
    for (key in c(
        "colUp", "colDown", "colNs", "gradientLow", "gradientHigh",
        "positiveColor", "negativeColor", "condition1Color", "condition2Color", "connectingLineColor"
    )) {
        params[[key]] <- normalize_visualization_color(params[[key]], color_defaults[[key]])
    }
    params$condition1Column <- if (is.null(params$condition1Column)) NULL else as.character(params$condition1Column[[1]])
    params$condition2Column <- if (is.null(params$condition2Column)) NULL else as.character(params$condition2Column[[1]])
    params
}

save_plot_assets <- function(plot_obj, output_dir, stem, width, height, png_dpi) {
    if (is.null(output_dir) || !nzchar(as.character(output_dir))) stop("Visualization output directory is missing")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
    paths <- list(
        svg = file.path(output_dir, paste0(stem, ".svg")),
        png = file.path(output_dir, paste0(stem, ".png")),
        pdf = file.path(output_dir, paste0(stem, ".pdf"))
    )

    render_with_device <- function(open_device) {
        previous_device <- grDevices::dev.cur()
        tryCatch(
            {
                open_device()
                print(plot_obj)
            },
            finally = {
                if (grDevices::dev.cur() != previous_device) grDevices::dev.off()
            }
        )
    }

    render_svg_in_clean_process <- function() {
        renderer <- file.path(analysis_root, "inst", "worker", "renderers", "svg_renderer.R")
        if (!file.exists(renderer)) stop("SVG renderer is missing: ", renderer, call. = FALSE)
        plot_snapshot <- tempfile(pattern = "atop-svg-", tmpdir = output_dir, fileext = ".rds")
        on.exit(unlink(plot_snapshot, force = TRUE), add = TRUE)
        saveRDS(plot_obj, plot_snapshot, compress = FALSE)
        output <- system2(
            file.path(R.home("bin"), "Rscript"),
            args = c(
                "--vanilla",
                shQuote(renderer),
                shQuote(plot_snapshot),
                shQuote(paths$svg),
                as.character(width),
                as.character(height)
            ),
            stdout = TRUE,
            stderr = TRUE
        )
        status <- attr(output, "status") %||% 0L
        if (!identical(as.integer(status), 0L) || !file.exists(paths$svg)) {
            stop("Failed to render SVG: ", paste(output, collapse = "\n"), call. = FALSE)
        }
    }

    render_svg_in_clean_process()

    render_with_device(function() {
        if (requireNamespace("ragg", quietly = TRUE)) {
            ragg::agg_png(paths$png, width = width, height = height, units = "in", res = png_dpi)
        } else {
            grDevices::png(paths$png, width = round(width * png_dpi), height = round(height * png_dpi), res = png_dpi)
        }
    })

    render_with_device(function() {
        if (identical(Sys.info()[["sysname"]], "Darwin") && capabilities("aqua")) {
            grDevices::quartz(type = "pdf", file = paths$pdf, width = width, height = height, family = "Helvetica")
        } else {
            grDevices::pdf(paths$pdf, width = width, height = height, onefile = FALSE, family = "Helvetica")
        }
    })

    lapply(paths, normalizePath, mustWork = TRUE)
}

visualization_title <- function(params) {
    params$title %||% switch(
        params$plotType,
        gene_volcano = "Volcano",
        gene_ranking = paste("Top", params$direction, "gene scores"),
        gsea_lollipop = paste(if (params$direction == "negative") "Negative" else "Positive", "GSEA pathways"),
        gsea_running = paste(if (params$direction == "negative") "Negative" else "Positive", "GSEA running plot"),
        sgrna_paired = paste("sgRNA counts:", params$targetGene %||% "selected gene")
    )
}

create_gene_plot <- function(plot_type, params, title) {
    gene_data <- state$data_processing$gene_summary_data
    gene_column <- state$data_processing$params$gene_col %||% "Gene"
    atop_require_columns(gene_data, gene_column, "Gene summary")
    gene_data$Gene <- gene_data[[gene_column]]
    switch(
        plot_type,
        gene_volcano = generate_volcano_plot(
            data = gene_data,
            p_value_cutoff = params$pValueCutoff,
            score_cutoff = params$scoreCutoff,
            label_genes = params$labelGenes,
            col_up = params$colUp,
            col_down = params$colDown,
            col_ns = params$colNs,
            plot_title = title,
            font_family = "Helvetica",
            base_font_size = params$fontSizePt
        ),
        gene_ranking = generate_ranking_plot(
            data = gene_data,
            type = params$direction,
            top_n = params$topN,
            gradient_low = params$gradientLow,
            gradient_high = params$gradientHigh,
            known_genes = params$knownGenes,
            plot_title = title,
            font_family = "Helvetica",
            base_font_size = params$fontSizePt
        )
    )
}

create_gsea_plot <- function(plot_type, params, title) {
    if (is.null(state$gsea)) stop("Run GSEA before creating a pathway visualization")
    negative <- identical(params$direction, "negative")
    result_table <- if (negative) state$gsea$gsea_results_negative_df else state$gsea$gsea_results_positive_df
    if (plot_type == "gsea_lollipop") {
        return(generate_gsea_lollipop_plot(
            gsea_data_input_df = result_table,
            analysis_type_label = if (negative) "Negative" else "Positive",
            num_pathways_to_plot = params$topN,
            plot_title_main = title,
            pathway_wrap_width = params$pathwayWrapWidth,
            positive_color = params$positiveColor,
            negative_color = params$negativeColor,
            font_family = "Helvetica",
            base_font_size = params$fontSizePt
        ))
    }
    if (!requireNamespace("GseaVis", quietly = TRUE)) stop("Package 'GseaVis' is required for running plots")
    gsea_object <- if (negative) state$gsea$gsea_object_negative else state$gsea$gsea_object_positive
    if (is.null(gsea_object)) stop("No GSEA object is available for the selected direction")
    pathway_ids <- params$pathwayIds
    if (length(pathway_ids) == 0L && !is.null(result_table) && nrow(result_table) > 0L) {
        pathway_ids <- head(as.character(result_table$ID), 1L)
    }
    generate_gsea_enrichment_plot_atop(
        gsea_object = gsea_object,
        plot_title = params$title,
        pathway_ids = pathway_ids,
        pathway_colours = params$pathwayColors,
        highlighted_genes = params$highlightedGenes,
        term_width = params$gseaTermWidth,
        legend_position = c(params$legendX, params$legendY),
        subplot_type = params$subplotType,
        base_font_size = params$fontSizePt,
        font_family = "Helvetica",
        direction = params$direction,
        show_nes = params$showNes,
        show_pvalue = params$showPvalue,
        pvalue_type = params$pvalueType,
        show_rank_explanation = params$showRankExplanation,
        rank_positive_label = params$rankPositiveLabel,
        rank_negative_label = params$rankNegativeLabel
    )
}

create_sgrna_plot <- function(params, title) {
    normalized <- state$data_processing$normalized_counts
    analysis_params <- state$data_processing$params
    if (is.null(normalized)) stop("Normalized counts are not available for paired plotting")
    if (is.null(params$targetGene) || !nzchar(params$targetGene)) stop("Select a target gene for the sgRNA paired plot")
    select_condition_column <- function(explicit_column, configured_columns, label) {
        configured_columns <- as.character(unlist(configured_columns %||% character(), use.names = FALSE))
        fallback_column <- if (length(configured_columns) > 0L) configured_columns[[1]] else NULL
        selected_column <- explicit_column %||% fallback_column
        if (is.null(selected_column) || !nzchar(selected_column)) {
            stop("Select a ", label, " count column for the sgRNA paired plot", call. = FALSE)
        }
        if (!selected_column %in% names(normalized)) {
            stop("The selected ", label, " count column is not available: ", selected_column, call. = FALSE)
        }
        selected_column
    }
    condition_one <- select_condition_column(params$condition1Column, analysis_params$condition1_replicate_cols, "first group")
    condition_two <- select_condition_column(params$condition2Column, analysis_params$condition2_replicate_cols, "second group")
    if (identical(condition_one, condition_two)) stop("Select two different count columns for the sgRNA paired plot", call. = FALSE)
    generate_sgrna_paired_plot(
        normalized_counts_df = normalized,
        target_gene_id = params$targetGene,
        sgrna_id_col = as.character(analysis_params$gRNA_col),
        gene_col = as.character(analysis_params$gene_col),
        condition_1_col = condition_one,
        condition_2_col = condition_two,
        color_condition_1 = params$condition1Color,
        color_condition_2 = params$condition2Color,
        connecting_line_color = params$connectingLineColor,
        plot_title = title,
        font_family = "Helvetica",
        base_font_size = params$fontSizePt
    )
}

run_visualization <- function(request) {
    source_modules(c("analysis", "gsea", "plotting"))
    if (is.null(state$data_processing)) stop("Complete data processing before creating a visualization")
    params <- resolve_visualization_parameters(request$parameters)
    supported <- c("gene_volcano", "gene_ranking", "gsea_lollipop", "gsea_running", "sgrna_paired")
    if (!params$plotType %in% supported) stop("Unsupported visualization type: ", params$plotType)
    title <- visualization_title(params)
    plot_obj <- if (params$plotType %in% c("gene_volcano", "gene_ranking")) {
        create_gene_plot(params$plotType, params, title)
    } else if (params$plotType %in% c("gsea_lollipop", "gsea_running")) {
        create_gsea_plot(params$plotType, params, title)
    } else {
        create_sgrna_plot(params, title)
    }
    if (params$plotType == "sgrna_paired" && params$widthCm < 5) {
        plot_obj <- plot_obj + ggplot2::theme(legend.position = "none")
    }
    stem <- paste0(
        "atop-", gsub("[^A-Za-z0-9_-]", "-", params$plotType), "-",
        gsub("[^A-Za-z0-9_-]", "-", request$id %||% "plot")
    )
    assets <- save_plot_assets(plot_obj, request$artifactDir, stem, params$width, params$height, params$pngDpi)
    list(status = "completed", plotType = params$plotType, title = title, assets = assets)
}
