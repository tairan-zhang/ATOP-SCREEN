






















geomean_custom_for_gsea <- function(x, na.rm = TRUE) {
    if (all(is.na(x))) {
        return(NA_real_)
    }
    valid_x <- x[x > 0 & !is.na(x) & is.finite(x)]
    if (length(valid_x) == 0) {
        return(NA_real_)
    }
    exp(mean(log(valid_x), na.rm = na.rm))
}

#' Run the canonical ATOP analysis workflow.
#'
#' This is intentionally an orchestration function: domain logic lives in
#' validation.R, preprocessing.R, scoring.R and statistics.R.
perform_crispr_screen_analysis <- function(raw_data_df,
                                           gRNA_col = "gRNA",
                                           gene_col = "Gene",
                                           sequence_col = NULL,
                                           condition1_replicate_cols,
                                           condition2_replicate_cols,
                                           N_perm = 1000,
                                           pseudo_count_lfc = 1.0,
                                           min_sgrna_threshold = 3,
                                           user_engine_choice = "cpp",
                                           skip_normalization = FALSE,
                                           diff_score_col1 = NULL,
                                           diff_score_col2 = NULL,
                                           progress_callback = NULL,
                                           initial_progress_value_abs = 0,
                                           ncs_delta0 = 0.1, ncs_delta1 = 0.3,
                                           ncs_delta2 = 0.2, ncs_threshold = 0,
                                           comparison_mode = "auto", permutation_seed = 42L) {
    params <- atop_validate_request(
        raw_data_df, gRNA_col, gene_col, sequence_col, condition1_replicate_cols,
        condition2_replicate_cols, N_perm, pseudo_count_lfc, min_sgrna_threshold,
        user_engine_choice, skip_normalization, diff_score_col1, diff_score_col2,
        initial_progress_value_abs, comparison_mode
    )
    params$ncs_delta0 <- atop_scalar(ncs_delta0, "NCS direction dead zone", 0)
    params$ncs_delta1 <- atop_scalar(ncs_delta1, "NCS magnitude threshold", 0)
    params$ncs_delta2 <- atop_scalar(ncs_delta2, "NCS stability scale", 0)
    params$ncs_threshold <- atop_scalar(ncs_threshold, "NCS cutoff", 0)
    if (params$ncs_delta1 == 0 || params$ncs_delta2 == 0) stop("NCS magnitude and stability scales must be positive.")
    if (params$ncs_threshold > 1) stop("NCS cutoff must be between zero and one.")
    params$ncs_enabled <- params$comparison_mode %in% c("paired", "shared_baseline")
    if (!params$ncs_enabled) params$ncs_threshold <- 0
    if (length(permutation_seed) != 1 || !is.numeric(permutation_seed) || !is.finite(permutation_seed) ||
        permutation_seed > .Machine$integer.max) stop("Seed must be an integer from 0 to 2147483647.")
    params$permutation_seed <- atop_scalar(permutation_seed, "Permutation seed", 0, TRUE)
    params$top_n_fraction <- 2/3
    params$top_n_rule <- "ceiling(2 * retained_guides / 3)"
    progress <- atop_progress_reporter(progress_callback)
    progress("Preparing ATOP input data.")
    prepared <- atop_prepare_input(raw_data_df, params, progress)
    guide_scores <- atop_score_guides(prepared, params, progress)
    filtered <- atop_filter_ncs(guide_scores$ranked, params)
    progress("Calculating adaptive Top-N gene scores.")
    gene_summary <- atop_score_genes(filtered$ranked, prepared$guide_counts, params$gene_col, params$min_sgrna_threshold)
    eligible <- filtered$ranked[[params$gene_col]] %in% gene_summary[[params$gene_col]]
    gene_summary <- atop_attach_permutation_pvalues(filtered$ranked[eligible, , drop = FALSE], gene_summary, params, progress)
    params$permutation_engine_used <- attr(gene_summary, "permutation_engine") %||% "skipped"
    gene_summary <- atop_apply_bh_correction(gene_summary)
    progress("ATOP analysis completed.")
    list(
        processed_sg_data = filtered$ranked,
        ncs_excluded = filtered$excluded,
        sg_ncs_filter_qc = filtered$qc,
        gene_summary_data = gene_summary,
        normalized_counts = guide_scores$normalized,
        filtered_genes_sgrna_data = prepared$removed,
        params = params
    )
}
