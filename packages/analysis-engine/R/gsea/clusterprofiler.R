























prepare_ranked_gene_list_gsea <- function(df, score_column, gene_id_col) {

    df_filtered <- df[!is.na(df[[score_column]]) & !is.na(df[[gene_id_col]]), ]
    if (nrow(df_filtered) == 0) {
        return(NULL)
    }


    if (any(duplicated(df_filtered[[gene_id_col]]))) {
        df_filtered <- df_filtered %>%
            group_by(!!sym(gene_id_col)) %>%
            slice_max(order_by = abs(!!sym(score_column)), n = 1, with_ties = FALSE) %>%
            ungroup()
    }


    gene_list_vector <- df_filtered[[score_column]]
    names(gene_list_vector) <- df_filtered[[gene_id_col]]


    gene_list_vector <- gene_list_vector[!is.na(names(gene_list_vector)) & !is.na(gene_list_vector)]

    if (length(gene_list_vector) == 0) {
        return(NULL)
    }



    sorted_gene_list <- sort(gene_list_vector, decreasing = TRUE)

    return(sorted_gene_list)
}

perform_gsea_analysis_clusterprofiler <- function(gene_summary_df,
                                                  gmt_file_path,
                                                  gene_id_column,
                                                  clp_n_permutations = 1000,
                                                  clp_pvalue_cutoff = 1,
                                                  clp_min_gene_set_size = 15,
                                                  max_gene_set_size = 500,
                                                  clp_seed_logical = TRUE,
                                                  clp_by = "fgsea",
                                                  progress_callback = NULL, permutation_seed = 123L) {
    report_progress <- function(message) {
        if (is.function(progress_callback)) {
            progress_callback(message)
        } else {
            cat(message, "\n")
        }
        invisible(NULL)
    }

    report_progress("Starting GSEA analysis...")
    report_progress("Reading GMT file...")
    term2gene_df <- tryCatch(
        {
            clusterProfiler::read.gmt(gmt_file_path)
        },
        error = function(e) stop(paste("Error reading GMT file:", e$message))
    )


    if (!gene_id_column %in% names(gene_summary_df)) stop(paste("Gene ID column '", gene_id_column, "' not found in gene summary."))
    if (!"GenePositiveScore" %in% names(gene_summary_df)) stop("'GenePositiveScore' not found in gene summary.")
    if (!"GeneNegativeScore" %in% names(gene_summary_df)) stop("'GeneNegativeScore' not found in gene summary.")


    pos_scores <- gene_summary_df$GenePositiveScore[!is.na(gene_summary_df$GenePositiveScore)]
    neg_scores <- gene_summary_df$GeneNegativeScore[!is.na(gene_summary_df$GeneNegativeScore)]

    report_progress("Preparing gene lists for GSEA analysis...")



    gene_list_positive <- prepare_ranked_gene_list_gsea(gene_summary_df, "GenePositiveScore", gene_id_column)



    gene_list_negative <- prepare_ranked_gene_list_gsea(gene_summary_df, "GeneNegativeScore", gene_id_column)

    run_single_gsea <- function(gene_list, term2gene, score_type_label) {
        report_progress(paste("Performing GSEA for:", score_type_label))
        gsea_s4_result <- NULL
        if (!is.null(gene_list) && length(gene_list) > 0) {
            if (length(permutation_seed) != 1 || !is.finite(permutation_seed) || permutation_seed < 0 || permutation_seed > .Machine$integer.max || permutation_seed != floor(permutation_seed)) stop("Invalid GSEA permutation seed.")
            if (isTRUE(clp_seed_logical)) {
                set.seed(as.integer(permutation_seed))
            }
            gsea_s4_result <- tryCatch(
                {
                    clusterProfiler::GSEA(
                        geneList = gene_list,
                        TERM2GENE = term2gene,
                        minGSSize = clp_min_gene_set_size,
                        maxGSSize = max_gene_set_size,
                        pvalueCutoff = clp_pvalue_cutoff,
                        nPerm = clp_n_permutations,
                        method = if (identical(clp_by, "permute")) "permute" else "multilevel",
                        verbose = FALSE
                    )
                },
                error = function(e) {
                    stop("GSEA failed: ", conditionMessage(e), call. = FALSE)
                }
            )
        } else {
            return(list(table = data.frame(), object = NULL))
        }
        report_progress(paste("GSEA for", score_type_label, "completed."))
        if (is.null(gsea_s4_result) || nrow(as.data.frame(gsea_s4_result)) == 0) {
            return(list(table = data.frame(), object = NULL))
        } else {
            return(list(table = as.data.frame(gsea_s4_result), object = gsea_s4_result))
        }
    }


    results_positive <- run_single_gsea(gene_list_positive, term2gene_df, "GenePositiveScore")


    results_negative <- run_single_gsea(gene_list_negative, term2gene_df, "GeneNegativeScore")

    report_progress("GSEA analysis function completed.")

    return(list(
        gsea_results_positive_df = results_positive$table,
        gsea_results_negative_df = results_negative$table,
        gsea_object_positive = results_positive$object,
        gsea_object_negative = results_negative$object,
        gene_list_positive = gene_list_positive,
        gene_list_negative = gene_list_negative
    ))
}
