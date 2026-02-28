# SPDX-License-Identifier: GPL-3.0-or-later
# gsea_functions.R
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


# gsea_functions.R
# GSEA Analysis Functions Collection

# --- GSEA Analysis Function ---
prepare_ranked_gene_list_gsea <- function(df, score_column, gene_id_col, progress_callback) {
    # Filter out NA values
    df_filtered <- df[!is.na(df[[score_column]]) & !is.na(df[[gene_id_col]]), ]
    if (nrow(df_filtered) == 0) {
        return(NULL)
    }

    # Handle duplicate gene IDs: select the score with the largest absolute value
    if (any(duplicated(df_filtered[[gene_id_col]]))) {
        df_filtered <- df_filtered %>%
            group_by(!!sym(gene_id_col)) %>%
            slice_max(order_by = abs(!!sym(score_column)), n = 1, with_ties = FALSE) %>%
            ungroup()
    }

    # Create named gene score vector
    gene_list_vector <- df_filtered[[score_column]]
    names(gene_list_vector) <- df_filtered[[gene_id_col]]

    # Remove any remaining NA values
    gene_list_vector <- gene_list_vector[!is.na(names(gene_list_vector)) & !is.na(gene_list_vector)]

    if (length(gene_list_vector) == 0) {
        return(NULL)
    }


    # Sort by score in descending order directly
    sorted_gene_list <- sort(gene_list_vector, decreasing = TRUE)

    return(sorted_gene_list)
}

perform_gsea_analysis_clusterProfiler_shiny <- function(gene_summary_df,
                                                        gmt_file_path,
                                                        gene_id_column,
                                                        clp_n_permutations = 1000,
                                                        clp_pvalue_cutoff = 1,
                                                        clp_min_gene_set_size = 15,
                                                        max_gene_set_size = 500,
                                                        clp_seed_logical = TRUE,
                                                        clp_by = "fgsea",
                                                        shiny_session = NULL,
                                                        initial_progress_value_abs = 0.5) {
    func_current_progress_abs <- initial_progress_value_abs
    gsea_progress_msg <- function(msg, val_increment_abs = 0) {
        if (!is.null(shiny_session)) {
            if (val_increment_abs > 0) {
                func_current_progress_abs <<- func_current_progress_abs + val_increment_abs
                shiny::setProgress(value = min(func_current_progress_abs, 0.99), message = msg, detail = paste0(round(min(func_current_progress_abs, 0.99) * 100), "%"), session = shiny_session)
            }
        } else {
            cat(paste(msg, "\n"))
        }
    }

    total_span_for_gsea <- 0.99 - initial_progress_value_abs
    num_gsea_major_steps <- 6
    gsea_step_unit <- total_span_for_gsea / num_gsea_major_steps

    gsea_progress_msg("Starting GSEA analysis...")

    gsea_progress_msg("Reading GMT file...", val_increment_abs = gsea_step_unit)
    term2gene_df <- tryCatch(
        {
            clusterProfiler::read.gmt(gmt_file_path)
        },
        error = function(e) stop(paste("Error reading GMT file:", e$message))
    )

    # Validate required columns exist
    if (!gene_id_column %in% names(gene_summary_df)) stop(paste("Gene ID column '", gene_id_column, "' not found in gene summary."))
    if (!"GenePositiveScore" %in% names(gene_summary_df)) stop("'GenePositiveScore' not found in gene summary.")
    if (!"GeneNegativeScore" %in% names(gene_summary_df)) stop("'GeneNegativeScore' not found in gene summary.")

    # Silently validate data (do not show detail)
    pos_scores <- gene_summary_df$GenePositiveScore[!is.na(gene_summary_df$GenePositiveScore)]
    neg_scores <- gene_summary_df$GeneNegativeScore[!is.na(gene_summary_df$GeneNegativeScore)]

    gsea_progress_msg("Preparing gene lists for GSEA analysis...", val_increment_abs = gsea_step_unit)

    # Use GenePositiveScore for positive GSEA analysis
    # High scores will be at the top of the ranked list, indicating promoting effects
    gene_list_positive <- prepare_ranked_gene_list_gsea(gene_summary_df, "GenePositiveScore", gene_id_column, gsea_progress_msg)

    # Use GeneNegativeScore for negative GSEA analysis
    # Raw scores are used directly for GSEA, negative values indicate inhibitory effects
    gene_list_negative <- prepare_ranked_gene_list_gsea(gene_summary_df, "GeneNegativeScore", gene_id_column, gsea_progress_msg)

    run_single_gsea_shiny <- function(gene_list, term2gene, score_type_label) {
        gsea_progress_msg(paste("Performing GSEA for:", score_type_label), val_increment_abs = gsea_step_unit * 0.5)
        gsea_s4_result <- NULL
        if (!is.null(gene_list) && length(gene_list) > 0) {
            seed_to_use <- FALSE
            if (isTRUE(clp_seed_logical)) {
                set.seed(123)
                seed_to_use <- TRUE
            }
            gsea_s4_result <- tryCatch(
                {
                    GSEA(
                        geneList = gene_list,
                        TERM2GENE = term2gene,
                        minGSSize = clp_min_gene_set_size,
                        maxGSSize = max_gene_set_size,
                        pvalueCutoff = clp_pvalue_cutoff,
                        nPermSimple = clp_n_permutations, # Changed from nPerm to nPermSimple based on common use with fgsea
                        verbose = FALSE,
                        seed = seed_to_use,
                        by = clp_by
                    )
                },
                error = function(e) {
                    return(list(table = data.frame(), object = NULL)) # Ensure list structure on error
                }
            )
        } else {
            return(list(table = data.frame(), object = NULL)) # Ensure list structure if skipped
        }
        gsea_progress_msg(paste("GSEA for", score_type_label, "completed."), val_increment_abs = gsea_step_unit * 0.5)
        if (is.null(gsea_s4_result) || nrow(as.data.frame(gsea_s4_result)) == 0) {
            return(list(table = data.frame(), object = NULL))
        } else {
            return(list(table = as.data.frame(gsea_s4_result), object = gsea_s4_result))
        }
    }

    # Run GSEA analysis based on GenePositiveScore
    results_positive <- run_single_gsea_shiny(gene_list_positive, term2gene_df, "GenePositiveScore")

    # Run GSEA analysis based on GeneNegativeScore
    results_negative <- run_single_gsea_shiny(gene_list_negative, term2gene_df, "GeneNegativeScore")

    gsea_progress_msg("GSEA analysis function completed.")

    return(list(
        gsea_results_positive_df = results_positive$table,
        gsea_results_negative_df = results_negative$table,
        gsea_object_positive = results_positive$object,
        gsea_object_negative = results_negative$object,
        gene_list_positive = gene_list_positive, # Return raw gene list
        gene_list_negative = gene_list_negative # Return raw gene list
    ))
}
