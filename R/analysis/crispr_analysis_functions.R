# SPDX-License-Identifier: GPL-3.0-or-later
# crispr_analysis_functions.R
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
                                           shiny_session = NULL,
                                           initial_progress_value_abs = 0) {
    N_perm <- validate_permutation_integer(N_perm, "N_perm", 0L)
    min_sgrna_threshold <- validate_permutation_integer(min_sgrna_threshold, "min_sgrna_threshold")
    if (length(pseudo_count_lfc) != 1L || !is.finite(pseudo_count_lfc) || pseudo_count_lfc < 0) {
        stop("pseudo_count_lfc must be a finite nonnegative number.")
    }
    func_current_progress_abs <- initial_progress_value_abs
    progress_msg <- function(msg, val_increment = 0) {
        if (!is.null(shiny_session)) {
            if (val_increment > 0) {
                func_current_progress_abs <<- func_current_progress_abs + val_increment
                shiny::setProgress(value = min(func_current_progress_abs, 0.49), message = msg, detail = paste0(round(min(func_current_progress_abs, 0.49) * 100), "%"), session = shiny_session)
            } else {
                cat(paste("[PROGRESS]", msg, "\n"))
            }
        } else {
            cat(paste(msg, "\n"))
        }
    }
    total_span_for_this_func <- 0.49 - initial_progress_value_abs
    num_major_steps_approx <- 18
    progress_step_increment <- total_span_for_this_func / num_major_steps_approx

    if (skip_normalization) {
        if (is.null(diff_score_col1) || is.null(diff_score_col2)) {
            stop("Skip normalization mode requires both diff_score_col1 and diff_score_col2 to be specified.")
        }
        if (!diff_score_col1 %in% colnames(raw_data_df) || !diff_score_col2 %in% colnames(raw_data_df)) {
            stop(paste0("Specified diff_score columns not found in data. Looking for: ", diff_score_col1, ", ", diff_score_col2))
        }
        progress_msg(paste("Skip normalization mode: Using", diff_score_col1, "-", diff_score_col2), val_increment = progress_step_increment / 2)
    } else {
        progress_msg("Validating replicate column inputs...")
        if (missing(condition1_replicate_cols) || length(condition1_replicate_cols) == 0) stop("Condition (numerator) replicate cols missing.")
        if (missing(condition2_replicate_cols) || length(condition2_replicate_cols) == 0) stop("Control (denominator) replicate cols missing.")

        n_reps_condition1 <- length(condition1_replicate_cols)
        n_reps_condition2 <- length(condition2_replicate_cols)
        if (n_reps_condition1 != n_reps_condition2) {
            stop(paste(
                "Unequal number of replicates: Condition (", n_reps_condition1,
                "), Control (", n_reps_condition2, ")"
            ))
        }
        n_reps <- n_reps_condition1
        progress_msg(paste("Using", n_reps, "replicates."), val_increment = progress_step_increment / 2)
        all_condition_replicate_cols <- c(condition1_replicate_cols, condition2_replicate_cols)
    }

    if (skip_normalization) {
        all_condition_replicate_cols <- character(0)
    }

    id_cols <- c(gRNA_col, gene_col)
    if (!is.null(sequence_col) && nzchar(sequence_col)) {
        id_cols <- c(id_cols, sequence_col)
    } else {
        sequence_col <- NULL
    }

    progress_msg("Step I: Data Initial Preparation", val_increment = progress_step_increment)
    data_df <- as.data.frame(raw_data_df)

    if (skip_normalization) {
        all_needed_cols_from_file <- c(id_cols, diff_score_col1, diff_score_col2)
    } else {
        all_needed_cols_from_file <- c(id_cols, all_condition_replicate_cols)
    }

    missing_cols <- setdiff(all_needed_cols_from_file, names(data_df))
    if (length(missing_cols) > 0) stop(paste("Missing columns in input:", paste(missing_cols, collapse = ", ")))
    data_df <- data_df[, all_needed_cols_from_file, drop = FALSE]
    for (column in c(gRNA_col, gene_col)) {
        labels <- as.character(data_df[[column]])
        if (any(is.na(labels) | !nzchar(trimws(labels)))) stop("Missing identifiers in ", column)
    }
    if (anyDuplicated(data_df[c(gene_col, gRNA_col)])) stop("Duplicate sgRNA identifiers within a gene.")

    if (!skip_normalization) {
        for (col in all_condition_replicate_cols) {
            data_df[[col]] <- suppressWarnings(as.numeric(data_df[[col]]))
            if (any(is.infinite(data_df[[col]]) | data_df[[col]] < 0, na.rm = TRUE)) {
                stop("Count columns must contain finite nonnegative values.")
            }
        }

        na_present_initial <- FALSE

        for (col_name_idx in seq_along(all_condition_replicate_cols)) {
            if (any(is.na(data_df[[all_condition_replicate_cols[col_name_idx]]]))) {
                na_present_initial <- TRUE
                break
            }
        }

        if (na_present_initial) {
            for (col in all_condition_replicate_cols) {
                data_df[[col]][is.na(data_df[[col]])] <- 0
            }
        }

        data_df[all_condition_replicate_cols] <- data_df[all_condition_replicate_cols] + 1
    } else {
        data_df[[diff_score_col1]] <- suppressWarnings(as.numeric(data_df[[diff_score_col1]]))
        data_df[[diff_score_col2]] <- suppressWarnings(as.numeric(data_df[[diff_score_col2]]))
    }

    progress_msg("Step I completed.", val_increment = progress_step_increment)

    progress_msg("Step Ia: Analyzing sgRNA counts per gene...", val_increment = progress_step_increment)
    sgrna_counts_per_gene <- data_df %>%
        distinct(!!sym(gene_col), !!sym(gRNA_col)) %>%
        group_by(!!sym(gene_col)) %>%
        summarise(actual_sgrna_count = n(), .groups = "drop")

    low_sgrna_genes_df <- sgrna_counts_per_gene %>% filter(actual_sgrna_count < min_sgrna_threshold)

    filtered_sgrna_data <- NULL
    if (nrow(low_sgrna_genes_df) > 0) {
        genes_to_remove <- low_sgrna_genes_df[[gene_col]]
        progress_msg(sprintf("Removing %d genes with < %d sgRNAs", length(genes_to_remove), min_sgrna_threshold))

        filtered_sgrna_data <- data_df %>%
            filter(!!sym(gene_col) %in% genes_to_remove) %>%
            left_join(low_sgrna_genes_df, by = gene_col)

        data_df <- data_df %>% filter(!(!!sym(gene_col) %in% genes_to_remove))
        if (nrow(data_df) == 0) stop(sprintf("All data removed after filtering genes with < %d sgRNAs.", min_sgrna_threshold))
    }

    if (skip_normalization) {
        progress_msg("Skip Mode: Calculating differential score from provided columns...", val_increment = progress_step_increment * 4)

        data_sg_level_results <- data_df[, id_cols, drop = FALSE]
        data_sg_level_results$diff_score <- data_df[[diff_score_col1]] - data_df[[diff_score_col2]]

        data_sg_processed <- data_sg_level_results[is.finite(data_sg_level_results$diff_score), ]
        if (nrow(data_sg_processed) < nrow(data_sg_level_results)) {
            progress_msg("Removed sgRNAs with NA differential score.")
        }
        if (nrow(data_sg_processed) == 0) stop("No sgRNAs remaining after calculating differential score.")

        progress_msg("Step III: sgRNA Ranking (skip mode)...", val_increment = progress_step_increment)
        data_sg_ranked <- data_sg_processed %>%
            group_by(!!sym(gene_col)) %>%
            mutate(
                sgPositiveRank = rank(-diff_score, ties.method = "first"),
                sgNegativeRank = rank(diff_score, ties.method = "first")
            ) %>%
            ungroup()

        data_normalized <- NULL
        progress_msg("Skip mode: Steps Ib-III skipped. Using differential score directly.", val_increment = progress_step_increment)
    } else {
        progress_msg("Step Ib: Normalizing samples...", val_increment = progress_step_increment * 1.5)
        counts_matrix_all_reps <- as.matrix(data_df[all_condition_replicate_cols])
        x_hat_i <- apply(counts_matrix_all_reps, 1, geomean_custom_for_gsea)
        x_hat_i[is.na(x_hat_i) | (x_hat_i == 0 & !is.na(x_hat_i))] <- 1
        r_ij_matrix_all_reps <- sweep(counts_matrix_all_reps, 1, x_hat_i, "/")
        s_j_vec_all_reps <- apply(r_ij_matrix_all_reps, 2, median, na.rm = TRUE)
        s_j_vec_all_reps[is.na(s_j_vec_all_reps) | (s_j_vec_all_reps == 0 & !is.na(s_j_vec_all_reps))] <- 1
        normalized_counts_matrix_all_reps <- round(sweep(counts_matrix_all_reps, 2, s_j_vec_all_reps, "/"))
        data_normalized <- data_df[, id_cols, drop = FALSE]
        data_normalized[colnames(normalized_counts_matrix_all_reps)] <- as.data.frame(normalized_counts_matrix_all_reps)
        progress_msg("Step Ib: Normalization completed.", val_increment = progress_step_increment)

        progress_msg("Step II: Calculating LFCs...", val_increment = progress_step_increment * 1.5)
        delta_lfc_replicates_list <- vector("list", n_reps)
        for (k in 1:n_reps) {
            xi_cond1_k_norm <- data_normalized[[condition1_replicate_cols[k]]]
            xi_cond2_k_norm <- data_normalized[[condition2_replicate_cols[k]]]
            cond1_pseudo <- xi_cond1_k_norm + pseudo_count_lfc
            cond2_pseudo <- xi_cond2_k_norm + pseudo_count_lfc

            delta_LFC_i_k <- log2(cond1_pseudo / cond2_pseudo)
            delta_LFC_i_k[!is.finite(delta_LFC_i_k)] <- NA
            delta_lfc_replicates_list[[k]] <- delta_LFC_i_k
        }
        delta_lfc_matrix <- do.call(cbind, delta_lfc_replicates_list)
        colnames(delta_lfc_matrix) <- paste0("deltaLFC_Rep", 1:n_reps)

        data_sg_level_results <- data_df[, id_cols, drop = FALSE]
        data_sg_level_results$diff_score <- rowMeans(delta_lfc_matrix, na.rm = TRUE)
        data_sg_processed <- data_sg_level_results[is.finite(data_sg_level_results$diff_score), ]
        if (nrow(data_sg_processed) < nrow(data_sg_level_results)) progress_msg("Removed sgRNAs with NA diff_score.")
        if (nrow(data_sg_processed) == 0) stop("No sgRNAs remaining after LFC calculation.")
        progress_msg("Step II: LFC calculation completed.", val_increment = progress_step_increment)

        progress_msg("Step III: sgRNA Ranking...", val_increment = progress_step_increment)
        data_sg_ranked <- data_sg_processed %>%
            group_by(!!sym(gene_col)) %>%
            mutate(sgPositiveRank = rank(-diff_score, ties.method = "first"), sgNegativeRank = rank(diff_score, ties.method = "first")) %>%
            ungroup()

        data_sg_ranked <- cbind(data_sg_ranked, as.data.frame(delta_lfc_matrix)[is.finite(data_sg_level_results$diff_score), , drop = FALSE])

        progress_msg("Step III: sgRNA ranking completed.", val_increment = progress_step_increment)
    }

    progress_msg("Step IV: Gene Scoring & Ranking (Adaptive Top-N)...", val_increment = progress_step_increment)

    data_with_counts <- data_sg_ranked %>%
        group_by(!!sym(gene_col)) %>%
        mutate(n_sgrnas = n()) %>%
        ungroup()

    gene_positive_scores <- data_with_counts %>%
        filter(n_sgrnas >= min_sgrna_threshold) %>%
        group_by(!!sym(gene_col)) %>%
        arrange(sgPositiveRank) %>%
        mutate(k = ceiling(2 * n_sgrnas / 3)) %>%
        filter(row_number() <= k) %>%
        summarise(
            GenePositiveScore = mean(diff_score, na.rm = TRUE),
            n_sgrnas_used = first(k), .groups = "drop"
        )

    gene_negative_scores <- data_with_counts %>%
        filter(n_sgrnas >= min_sgrna_threshold) %>%
        group_by(!!sym(gene_col)) %>%
        arrange(diff_score) %>%
        mutate(k = ceiling(2 * n_sgrnas / 3)) %>%
        filter(row_number() <= k) %>%
        summarise(
            GeneNegativeScore = mean(diff_score, na.rm = TRUE),
            n_sgrnas_used = first(k), .groups = "drop"
        )

    gene_summary_df <- full_join(
        gene_positive_scores %>% rename(n_sgrnas_used_pos = n_sgrnas_used),
        gene_negative_scores %>% rename(n_sgrnas_used_neg = n_sgrnas_used),
        by = gene_col
    )

    gene_summary_df <- gene_summary_df %>%
        left_join(sgrna_counts_per_gene, by = gene_col)

    gene_summary_df <- gene_summary_df %>%
        arrange(desc(GenePositiveScore)) %>%
        mutate(GenePositiveRank = rank(-GenePositiveScore, ties.method = "first", na.last = "keep")) %>%
        arrange(GeneNegativeScore) %>%
        mutate(GeneNegativeRank = rank(GeneNegativeScore, ties.method = "first", na.last = "keep"))
    progress_msg("Step IV: Gene scoring completed.", val_increment = progress_step_increment)

    gene_summary_df$P_positive <- rep(NA_real_, nrow(gene_summary_df))
    gene_summary_df$P_negative <- rep(NA_real_, nrow(gene_summary_df))
    if (N_perm >= 10L && nrow(gene_summary_df) > 0L) {
        progress_msg("Calculating permutation P-values.", val_increment = progress_step_increment)
        permutation_start <- func_current_progress_abs
        report_permutation <- function(message) {
            progress_msg(message)
            if (!is.null(shiny_session)) {
                percentage <- as.numeric(sub(".*\\(([0-9.]+)%\\).*", "\\1", message))
                shiny::setProgress(
                    value = permutation_start + (0.49 - permutation_start) * percentage / 100,
                    message = message, session = shiny_session
                )
            }
        }
        permutation_results <- engine_selector(
            sgrna_info_dt = data_sg_ranked,
            gene_summary_dt = gene_summary_df,
            gene_col = gene_col,
            N_perm = N_perm,
            progress_callback = report_permutation,
            user_engine_choice = user_engine_choice,
            min_sgrna_threshold = min_sgrna_threshold,
            shiny_session = NULL
        )
        gene_summary_df$P_positive <- permutation_results$P_positive
        gene_summary_df$P_negative <- permutation_results$P_negative
    } else {
        progress_msg("Fewer than 10 permutations or no eligible genes; P-values skipped.")
    }
    for (direction in c("positive", "negative")) {
        p_column <- paste0("P_", direction)
        gene_summary_df[[paste0(p_column, "_adj_bh")]] <- p.adjust(gene_summary_df[[p_column]], method = "BH")
    }

    return(list(
        processed_sg_data = data_sg_ranked,
        gene_summary_data = gene_summary_df,
        normalized_counts = data_normalized,
        filtered_genes_sgrna_data = filtered_sgrna_data,
        params = list(
            gRNA_col = gRNA_col,
            gene_col = gene_col,
            sequence_col = sequence_col,
            analysis_tool = "ATOP",
            software_version = trimws(readLines("VERSION", warn = FALSE)[1])
        )
    ))
}
