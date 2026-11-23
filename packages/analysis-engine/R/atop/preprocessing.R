# Data preparation and guide-level scoring for the original ATOP method.

atop_progress_reporter <- function(progress_callback = NULL) {
    function(message) {
        if (is.function(progress_callback)) {
            progress_callback(message)
        } else {
            cat(message, "\n")
        }
        invisible(NULL)
    }
}

atop_prepare_input <- function(raw_data_df, params, progress) {
    id_columns <- c(params$gRNA_col, params$gene_col, params$sequence_col)
    id_columns <- id_columns[!is.null(id_columns) & nzchar(id_columns)]
    numeric_columns <- if (params$skip_normalization) c(params$diff_score_col1, params$diff_score_col2) else c(params$condition1_replicate_cols, params$condition2_replicate_cols)
    data <- raw_data_df[, c(id_columns, numeric_columns), drop = FALSE]
    for (column in numeric_columns) data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
    if (!params$skip_normalization) {
        had_missing_counts <- anyNA(data[numeric_columns])
        data[numeric_columns] <- lapply(data[numeric_columns], function(value) { value[is.na(value)] <- 0; value + 1 })
        progress(if (had_missing_counts) "Missing counts replaced with zero before normalization." else "Count columns prepared for normalization.")
    }
    guide_counts <- dplyr::summarise(
        dplyr::group_by(dplyr::distinct(data, .data[[params$gene_col]], .data[[params$gRNA_col]]), .data[[params$gene_col]]),
        actual_sgrna_count = dplyr::n(), .groups = "drop"
    )
    low_guides <- dplyr::filter(guide_counts, .data$actual_sgrna_count < params$min_sgrna_threshold)
    removed <- NULL
    if (nrow(low_guides) > 0L) {
        removed_genes <- low_guides[[params$gene_col]]
        removed <- dplyr::left_join(dplyr::filter(data, .data[[params$gene_col]] %in% removed_genes), low_guides, by = params$gene_col)
        data <- dplyr::filter(data, !.data[[params$gene_col]] %in% removed_genes)
    }
    if (nrow(data) == 0L) stop("All guides were removed by the minimum sgRNA threshold.", call. = FALSE)
    list(data = data, id_columns = id_columns, guide_counts = guide_counts, removed = removed)
}

atop_normalize_counts <- function(data, count_columns, id_columns) {
    matrix_counts <- as.matrix(data[count_columns])
    guide_geomeans <- apply(matrix_counts, 1, geomean_custom_for_gsea)
    guide_geomeans[is.na(guide_geomeans) | guide_geomeans == 0] <- 1
    factors <- apply(sweep(matrix_counts, 1, guide_geomeans, "/"), 2, stats::median, na.rm = TRUE)
    factors[is.na(factors) | factors == 0] <- 1
    normalized_matrix <- round(sweep(matrix_counts, 2, factors, "/"))
    normalized <- data[, id_columns, drop = FALSE]
    normalized[colnames(normalized_matrix)] <- as.data.frame(normalized_matrix)
    normalized
}

atop_score_guides <- function(prepared, params, progress) {
    data <- prepared$data
    if (params$skip_normalization) {
        guide_results <- data[, prepared$id_columns, drop = FALSE]
        guide_results$diff_score <- data[[params$diff_score_col1]] - data[[params$diff_score_col2]]
        normalized <- NULL
    } else {
        count_columns <- c(params$condition1_replicate_cols, params$condition2_replicate_cols)
        normalized <- atop_normalize_counts(data, count_columns, prepared$id_columns)
        guide_results <- data[, prepared$id_columns, drop = FALSE]
        endpoints <- params$condition1_replicate_cols
        baselines <- params$condition2_replicate_cols
        lfc <- function(endpoint, baseline) {
            delta <- log2((endpoint + params$pseudo_count_lfc) / (baseline + params$pseudo_count_lfc))
            delta[!is.finite(delta)] <- NA_real_
            delta
        }
        if (identical(params$comparison_mode, "group_means")) {
            guide_results$diff_score <- lfc(rowMeans(normalized[endpoints]), rowMeans(normalized[baselines]))
        } else {
            if (identical(params$comparison_mode, "shared_baseline")) baselines <- rep(baselines, length(endpoints))
            differences <- Map(function(a, b) lfc(normalized[[a]], normalized[[b]]), endpoints, baselines)
            lfc_matrix <- do.call(cbind, differences)
            colnames(lfc_matrix) <- paste0("deltaLFC_Rep", seq_along(endpoints))
            guide_results$diff_score <- rowMeans(lfc_matrix, na.rm = TRUE)
        }
    }
    valid_guides <- !is.na(guide_results$diff_score)
    guide_results <- guide_results[valid_guides, , drop = FALSE]
    if (nrow(guide_results) == 0L) stop("No valid sgRNA differential scores remain.", call. = FALSE)
    if (!params$skip_normalization && !identical(params$comparison_mode, "group_means")) {
        lfc <- lfc_matrix[valid_guides, , drop = FALSE]
        components <- atop_ncs(lfc, params)
        guide_results <- cbind(guide_results, as.data.frame(t(components)))
    }
    ranked <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(guide_results, .data[[params$gene_col]]), sgPositiveRank = rank(-.data$diff_score, ties.method = "first"), sgNegativeRank = rank(.data$diff_score, ties.method = "first")))
    if (!params$skip_normalization && !identical(params$comparison_mode, "group_means")) ranked <- cbind(ranked, as.data.frame(lfc_matrix)[valid_guides, , drop = FALSE])
    progress("sgRNA normalization, LFC calculation, and ranking completed.")
    list(ranked = ranked, normalized = normalized)
}

# Keep the R engine usable when the optional native library is unavailable.
atop_ncs <- function(lfc, params) {
    if (exists("calculate_ncs_cpp", mode = "function")) {
        return(calculate_ncs_cpp(lfc, params$ncs_delta0, params$ncs_delta1, params$ncs_delta2))
    }
    vapply(seq_len(nrow(lfc)), function(i) {
        x <- lfc[i, ]
        x <- x[is.finite(x)]
        if (!length(x)) return(c(dir_score = NA_real_, mag_score = NA_real_, stab_score = NA_real_, NCS = NA_real_))
        up <- sum(x > params$ncs_delta0)
        down <- sum(x < -params$ncs_delta0)
        direction <- max(up, down, length(x) - up - down) / length(x)
        magnitude <- min(1, median(abs(x)) / params$ncs_delta1)
        stability <- 1 / (1 + median(abs(x - median(x))) / params$ncs_delta2)
        c(dir_score = direction, mag_score = magnitude, stab_score = stability,
          NCS = direction * magnitude * stability)
    }, numeric(4))
}

atop_filter_ncs <- function(ranked, params) {
    if (params$ncs_threshold == 0 || is.null(ranked$NCS)) return(list(ranked = ranked, qc = NULL))
    keep <- is.na(ranked$NCS) | ranked$NCS >= params$ncs_threshold
    qc <- ranked
    qc$.retained <- keep
    qc <- dplyr::summarise(dplyr::group_by(qc, .data[[params$gene_col]]),
        total_sgrnas = dplyr::n(), retained_sgrnas = sum(.data$.retained),
        removed_sgrnas = total_sgrnas - retained_sgrnas,
        retention_rate = retained_sgrnas / total_sgrnas, .groups = "drop")
    list(ranked = ranked[keep, , drop = FALSE], excluded = ranked[!keep, , drop = FALSE], qc = dplyr::arrange(qc, dplyr::desc(.data$removed_sgrnas)))
}
