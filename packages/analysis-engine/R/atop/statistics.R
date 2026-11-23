# Permutation dispatch and multiple-testing correction for the ATOP workflow.

atop_available_cores <- function() {
    detected <- suppressWarnings(parallel::detectCores())
    if (length(detected) != 1L || !is.finite(detected) || detected < 1L) return(1L)
    max(1L, min(as.integer(detected) - 1L, 8L))
}

atop_attach_permutation_pvalues <- function(ranked_guides, gene_summary, params, progress) {
    empty_pvalues <- function(summary) {
        summary$P_positive <- rep(NA_real_, nrow(summary))
        summary$P_negative <- rep(NA_real_, nrow(summary))
        summary
    }
    if (params$N_perm < 10L) {
        progress("Fewer than 10 permutations requested; P-value calculation skipped.")
        return(empty_pvalues(gene_summary))
    }
    permutation_data <- dplyr::filter(dplyr::select(ranked_guides, dplyr::all_of(c(params$gene_col, "diff_score"))), !is.na(.data$diff_score))
    if (nrow(permutation_data) == 0L || nrow(gene_summary) == 0L) return(empty_pvalues(gene_summary))
    core_count <- atop_available_cores()
    progress(paste("Running", params$N_perm, "permutations with", core_count, "CPU cores using the", params$user_engine_choice, "engine."))
    started <- Sys.time()
    result <- engine_selector(
        sgrna_info_dt = permutation_data,
        gene_summary_dt = gene_summary,
        gene_col = params$gene_col,
        N_perm = params$N_perm,
        n_cores = core_count,
        progress_callback = function(message) progress(message),
        user_engine_choice = params$user_engine_choice,
        min_sgrna_threshold = params$min_sgrna_threshold, seed = params$permutation_seed
    )
    attr(gene_summary, "permutation_engine") <- result$engine
    gene_summary$P_positive <- result$P_positive
    gene_summary$P_negative <- result$P_negative
    progress(paste("Permutation test completed in", round(as.numeric(Sys.time() - started, units = "secs"), 2), "seconds."))
    gene_summary
}
