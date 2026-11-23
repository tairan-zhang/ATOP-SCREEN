# Adaptive Top-N gene scoring is kept independent from input and permutation IO.

atop_score_genes <- function(ranked_guides, guide_counts, gene_col, min_sgrna_threshold) {
    counted <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(ranked_guides, .data[[gene_col]]), n_sgrnas = dplyr::n()))
    score_direction <- function(data, descending, score_name) {
        dplyr::summarise(
            dplyr::filter(
                dplyr::mutate(
                    dplyr::arrange(dplyr::group_by(dplyr::filter(data, .data$n_sgrnas >= min_sgrna_threshold), .data[[gene_col]]), if (descending) dplyr::desc(.data$diff_score) else .data$diff_score),
                    k = ceiling(2 * .data$n_sgrnas / 3)
                ),
                dplyr::row_number() <= .data$k
            ),
            !!score_name := mean(.data$diff_score, na.rm = TRUE), n_sgrnas_used = dplyr::first(.data$k), .groups = "drop"
        )
    }
    positive <- score_direction(counted, TRUE, "GenePositiveScore")
    negative <- score_direction(counted, FALSE, "GeneNegativeScore")
    summary <- dplyr::full_join(dplyr::rename(positive, n_sgrnas_used_pos = n_sgrnas_used), dplyr::rename(negative, n_sgrnas_used_neg = n_sgrnas_used), by = gene_col)
    summary <- dplyr::left_join(summary, guide_counts, by = gene_col)
    dplyr::mutate(
        dplyr::arrange(summary, .data$GeneNegativeScore),
        GenePositiveRank = rank(-.data$GenePositiveScore, ties.method = "first", na.last = "keep"),
        GeneNegativeRank = rank(.data$GeneNegativeScore, ties.method = "first", na.last = "keep")
    )
}

atop_apply_bh_correction <- function(summary) {
    for (direction in c("positive", "negative")) {
        source <- paste0("P_", direction)
        target <- paste0(source, "_adj_bh")
        summary[[target]] <- rep(NA_real_, nrow(summary))
        valid <- !is.na(summary[[source]])
        if (any(valid)) summary[[target]][valid] <- stats::p.adjust(summary[[source]][valid], method = "BH")
    }
    summary
}
