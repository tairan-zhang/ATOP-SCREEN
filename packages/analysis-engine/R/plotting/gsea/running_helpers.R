# Helpers shared by the ATOP GSEA running-plot pipeline.

format_gsea_pathway_name_atop <- function(x, width = 55) {
    value <- gsub("_", " ", toupper(gsub("\n", " ", as.character(x))))
    if (!is.null(width) && is.finite(width) && width > 0) {
        value <- stringr::str_wrap(value, width = width)
    }
    value
}

format_gsea_p_atop <- function(x) {
    if (!is.finite(x)) return("NA")
    if (x < 0.001) return("<0.001")
    formatC(x, digits = 3, format = "fg")
}

normalize_gsea_pathway_labels_atop <- function(gsea_object) {
    if (is.null(gsea_object) || !inherits(gsea_object, "gseaResult")) {
        stop("Input is not a valid gseaResult S4 object.", call. = FALSE)
    }

    original_ids <- as.character(gsea_object@result$ID)
    display_ids <- toupper(gsub("_", " ", original_ids, fixed = TRUE))
    id_map <- stats::setNames(display_ids, original_ids)

    gsea_object@result$ID <- unname(id_map[original_ids])
    if ("Description" %in% names(gsea_object@result)) {
        gsea_object@result$Description <- unname(id_map[original_ids])
    }
    rownames(gsea_object@result) <- unname(id_map[original_ids])

    remap_ids <- function(ids) {
        mapped <- unname(id_map[ids])
        mapped[is.na(mapped)] <- ids[is.na(mapped)]
        mapped
    }

    if (!is.null(gsea_object@geneSets) && length(gsea_object@geneSets) > 0L) {
        names(gsea_object@geneSets) <- remap_ids(names(gsea_object@geneSets))
    }
    if (!is.null(gsea_object@permScores)) {
        if (!is.null(rownames(gsea_object@permScores))) {
            rownames(gsea_object@permScores) <- remap_ids(rownames(gsea_object@permScores))
        }
        if (!is.null(colnames(gsea_object@permScores))) {
            colnames(gsea_object@permScores) <- remap_ids(colnames(gsea_object@permScores))
        }
    }
    if (!is.null(gsea_object@termsim)) {
        if (!is.null(rownames(gsea_object@termsim))) {
            rownames(gsea_object@termsim) <- remap_ids(rownames(gsea_object@termsim))
        }
        if (!is.null(colnames(gsea_object@termsim))) {
            colnames(gsea_object@termsim) <- remap_ids(colnames(gsea_object@termsim))
        }
    }

    list(object = gsea_object, id_map = id_map)
}

gsea_pathway_rows_atop <- function(gsea_object, pathway_ids) {
    rows <- as.data.frame(gsea_object@result, stringsAsFactors = FALSE)
    rows <- rows[rows$ID %in% pathway_ids, , drop = FALSE]
    rows[match(pathway_ids, rows$ID), , drop = FALSE]
}

resolve_gsea_highlighted_genes_atop <- function(gsea_object, pathway_ids, highlighted_genes = NULL) {
    requested <- unique(trimws(as.character(highlighted_genes %||% character())))
    requested <- requested[!is.na(requested) & nzchar(requested)]
    if (length(requested) == 0L) return(NULL)

    result_rows <- as.data.frame(gsea_object@result, stringsAsFactors = FALSE)
    result_rows <- result_rows[result_rows$ID %in% pathway_ids, , drop = FALSE]
    pathway_genes <- if ("core_enrichment" %in% names(result_rows)) {
        unique(unlist(strsplit(as.character(result_rows$core_enrichment), "/", fixed = TRUE), use.names = FALSE))
    } else {
        character()
    }
    pathway_genes <- pathway_genes[!is.na(pathway_genes) & nzchar(pathway_genes)]
    missing_genes <- setdiff(requested, pathway_genes)
    if (length(missing_genes) == length(requested)) {
        stop(
            "None of the specified highlight genes are present in the selected pathway(s): ",
            paste(requested, collapse = ", "), ".",
            call. = FALSE
        )
    }
    if (length(missing_genes) > 0L) {
        warning(
            "Ignored highlight genes not present in the selected pathway(s): ",
            paste(missing_genes, collapse = ", "), ".",
            call. = FALSE
        )
    }
    intersect(requested, pathway_genes)
}
