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


# CRISPR Screen Analysis Core Functions

# Geometric Mean Helper
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

# Main CRISPR Screen Analysis Function
perform_crispr_screen_analysis <- function(raw_data_df,
                                           gRNA_col = "gRNA",
                                           gene_col = "Gene",
                                           sequence_col = NULL,
                                           condition1_replicate_cols, # Treatment
                                           condition2_replicate_cols, # Control
                                           N_perm = 1000,
                                           pseudo_count_lfc = 1.0,
                                           min_sgrna_threshold = 3,
                                           user_engine_choice = "cpp",
                                           skip_normalization = FALSE,
                                           diff_score_col1 = NULL,
                                           diff_score_col2 = NULL,
                                           shiny_session = NULL,
                                           initial_progress_value_abs = 0) {
    func_current_progress_abs <- initial_progress_value_abs
    progress_msg <- function(msg, val_increment = 0) {
        if (!is.null(shiny_session)) {
            if (val_increment > 0) {
                func_current_progress_abs <<- func_current_progress_abs + val_increment
                shiny::setProgress(value = min(func_current_progress_abs, 0.49), message = msg, detail = paste0(round(min(func_current_progress_abs, 0.49) * 100), "%"), session = shiny_session)
            } else {
                # Simple notification: console only
                cat(paste("[PROGRESS]", msg, "\n"))
            }
        } else {
            cat(paste(msg, "\n"))
        }
    }
    total_span_for_this_func <- 0.49 - initial_progress_value_abs
    num_major_steps_approx <- 18
    progress_step_increment <- total_span_for_this_func / num_major_steps_approx


    # Skip normalization mode validation
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

    # Set empty replicate array for skip mode
    if (skip_normalization) {
        all_condition_replicate_cols <- character(0)
    }

    id_cols <- c(gRNA_col, gene_col)
    if (!is.null(sequence_col) && nzchar(sequence_col)) {
        id_cols <- c(id_cols, sequence_col)
    } else {
        sequence_col <- NULL # Ensure it's truly NULL if not provided or empty
    }

    progress_msg("Step I: Data Initial Preparation", val_increment = progress_step_increment)
    data_df <- raw_data_df

    # Determine which columns to check based on mode
    if (skip_normalization) {
        # Check for diff_score columns in skip mode
        all_needed_cols_from_file <- c(id_cols, diff_score_col1, diff_score_col2)
    } else {
        # Check for replicate columns in full pipeline mode
        all_needed_cols_from_file <- c(id_cols, all_condition_replicate_cols)
    }

    missing_cols <- setdiff(all_needed_cols_from_file, names(data_df))
    if (length(missing_cols) > 0) stop(paste("Missing columns in input:", paste(missing_cols, collapse = ", ")))
    data_df <- data_df[, all_needed_cols_from_file, drop = FALSE]

    # Only process replicate columns in full pipeline mode
    if (!skip_normalization) {
        for (col in all_condition_replicate_cols) {
            data_df[[col]] <- suppressWarnings(as.numeric(data_df[[col]]))
        }

        # Data Cleaning: Handle NAs and pseudo-counts
        na_present_initial <- FALSE
        # Check for NAs efficiently across all relevant columns
        for (col_name_idx in seq_along(all_condition_replicate_cols)) {
            if (any(is.na(data_df[[all_condition_replicate_cols[col_name_idx]]]))) {
                na_present_initial <- TRUE
                break
            }
        }

        if (na_present_initial) {
            detail_msg_na <- "NAs found in count columns. Replacing with 0."
            for (col in all_condition_replicate_cols) {
                data_df[[col]][is.na(data_df[[col]])] <- 0 # NAs become 0
            }
        }

        # Add 1 to all count columns (original values and former NAs that are now 0)
        data_df[all_condition_replicate_cols] <- data_df[all_condition_replicate_cols] + 1

        detail_msg_increment <- if (na_present_initial) {
            "Original NAs are now 1 (0+1). Other counts also incremented by 1."
        } else {
            "All count columns incremented by 1 (no NAs found initially)."
        }
    } else {
        # Skip mode: validate diff_score columns are numeric
        data_df[[diff_score_col1]] <- suppressWarnings(as.numeric(data_df[[diff_score_col1]]))
        data_df[[diff_score_col2]] <- suppressWarnings(as.numeric(data_df[[diff_score_col2]]))
        detail_msg_increment <- "Skip mode: diff_score columns loaded."
    }
    # if (!is.null(shiny_session)) showNotification(detail_msg_increment, type = "message", session = shiny_session) else cat(paste0(detail_msg_increment, "\\n"))

    sgrna_count_report_df <- NULL # This is for Step Ia
    progress_msg("Step I completed.", val_increment = progress_step_increment)

    # Step Ia: Filter genes based on sgRNA count threshold
    progress_msg("Step Ia: Analyzing sgRNA counts per gene...", val_increment = progress_step_increment)
    sgrna_counts_per_gene <- data_df %>%
        distinct(!!sym(gene_col), !!sym(gRNA_col)) %>%
        group_by(!!sym(gene_col)) %>%
        summarise(actual_sgrna_count = n(), .groups = "drop")


    low_sgrna_genes_df <- sgrna_counts_per_gene %>% filter(actual_sgrna_count < min_sgrna_threshold)

    # Store sgRNA data for filtered genes
    filtered_sgrna_data <- NULL
    if (nrow(low_sgrna_genes_df) > 0) {
        genes_to_remove <- low_sgrna_genes_df[[gene_col]]
        progress_msg(sprintf("Removing %d genes with < %d sgRNAs", length(genes_to_remove), min_sgrna_threshold))

        # Get all sgRNA data for filtered genes before removing
        filtered_sgrna_data <- data_df %>%
            filter(!!sym(gene_col) %in% genes_to_remove) %>%
            left_join(low_sgrna_genes_df, by = gene_col)

        # Remove these genes from analysis
        data_df <- data_df %>% filter(!(!!sym(gene_col) %in% genes_to_remove))
        if (nrow(data_df) == 0) stop(sprintf("All data removed after filtering genes with < %d sgRNAs.", min_sgrna_threshold))
    }

    # Store dropped genes for reporting
    sgrna_count_report_df <- low_sgrna_genes_df

    # Branch: Skip normalization mode vs. Full pipeline
    if (skip_normalization) {
        # Skip Mode: Use pre-calculated diff_score columns directly
        progress_msg("Skip Mode: Calculating ddiff_score from provided columns...", val_increment = progress_step_increment * 4)

        # Calculate ddiff_score = diff_score_col1 - diff_score_col2
        data_sg_level_results <- data_df[, id_cols, drop = FALSE]
        data_sg_level_results$diff_score <- data_df[[diff_score_col1]] - data_df[[diff_score_col2]]

        # Filter out NA values
        data_sg_processed <- data_sg_level_results[!is.na(data_sg_level_results$diff_score), ]
        if (nrow(data_sg_processed) < nrow(data_sg_level_results)) {
            progress_msg("Removed sgRNAs with NA ddiff_score.")
        }
        if (nrow(data_sg_processed) == 0) stop("No sgRNAs remaining after calculating ddiff_score.")

        # Step III equivalent: sgRNA Ranking
        progress_msg("Step III: sgRNA Ranking (skip mode)...", val_increment = progress_step_increment)
        data_sg_ranked <- data_sg_processed %>%
            group_by(!!sym(gene_col)) %>%
            mutate(
                sgPositiveRank = rank(-diff_score, ties.method = "first"),
                sgNegativeRank = rank(diff_score, ties.method = "first")
            ) %>%
            ungroup()

        # No normalized counts or LFC data in skip mode
        data_normalized <- NULL
        progress_msg("Skip mode: Steps Ib-III skipped. Using ddiff_score directly.", val_increment = progress_step_increment)
    } else {
        # Full Pipeline: Normalization + LFC calculation
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
        data_sg_processed <- data_sg_level_results[!is.na(data_sg_level_results$diff_score), ]
        if (nrow(data_sg_processed) < nrow(data_sg_level_results)) progress_msg("Removed sgRNAs with NA diff_score.")
        if (nrow(data_sg_processed) == 0) stop("No sgRNAs remaining after LFC calculation.")
        progress_msg("Step II: LFC calculation completed.", val_increment = progress_step_increment)

        progress_msg("Step III: sgRNA Ranking...", val_increment = progress_step_increment)
        data_sg_ranked <- data_sg_processed %>%
            group_by(!!sym(gene_col)) %>%
            mutate(sgPositiveRank = rank(-diff_score, ties.method = "first"), sgNegativeRank = rank(diff_score, ties.method = "first")) %>%
            ungroup()

        # Add LFC details (deltaLFC per replicate) to the sgRNA data
        data_sg_ranked <- cbind(data_sg_ranked, as.data.frame(delta_lfc_matrix)[!is.na(data_sg_level_results$diff_score), , drop = FALSE])

        progress_msg("Step III: sgRNA ranking completed.", val_increment = progress_step_increment)
    }

    progress_msg("Step IV: Gene Scoring & Ranking (Adaptive Top-N)...", val_increment = progress_step_increment)

    # Adaptive Scoring: Top-k (k = ceil(2n/3))

    data_with_counts <- data_sg_ranked %>%
        group_by(!!sym(gene_col)) %>%
        mutate(n_sgrnas = n()) %>%
        ungroup()

    # Process Positive Scores
    gene_positive_scores <- data_with_counts %>%
        filter(n_sgrnas >= min_sgrna_threshold) %>%
        group_by(!!sym(gene_col)) %>%
        arrange(sgPositiveRank) %>%
        mutate(k = ceiling(2 * n_sgrnas / 3)) %>%
        filter(row_number() <= k) %>% # take top k (filter handles varying k per group correctly)
        summarise(
            GenePositiveScore = mean(diff_score, na.rm = TRUE),
            n_sgrnas_used = first(k), .groups = "drop"
        )

    # Process Negative Scores
    gene_negative_scores <- data_with_counts %>%
        filter(n_sgrnas >= min_sgrna_threshold) %>%
        group_by(!!sym(gene_col)) %>%
        arrange(diff_score) %>% # Select most negative
        mutate(k = ceiling(2 * n_sgrnas / 3)) %>%
        filter(row_number() <= k) %>% # take top k most negative
        summarise(
            GeneNegativeScore = mean(diff_score, na.rm = TRUE),
            n_sgrnas_used = first(k), .groups = "drop"
        )
    # Merge Positive and Negative scores
    gene_summary_df <- full_join(
        gene_positive_scores %>% rename(n_sgrnas_used_pos = n_sgrnas_used),
        gene_negative_scores %>% rename(n_sgrnas_used_neg = n_sgrnas_used),
        by = gene_col
    )

    # Add sgRNA count information
    gene_summary_df <- gene_summary_df %>%
        left_join(sgrna_counts_per_gene, by = gene_col)

    gene_summary_df <- gene_summary_df %>%
        arrange(desc(GenePositiveScore)) %>%
        mutate(GenePositiveRank = rank(-GenePositiveScore, ties.method = "first", na.last = "keep")) %>%
        arrange(GeneNegativeScore) %>%
        mutate(GeneNegativeRank = rank(GeneNegativeScore, ties.method = "first", na.last = "keep"))
    progress_msg("Step IV: Gene scoring completed.", val_increment = progress_step_increment)

    if (N_perm > 0 && N_perm >= 10) {
        progress_msg(paste("Step V: P-value Calculation (Permutations: ", N_perm, ")..."), val_increment = progress_step_increment)
        perm_loop_total_increment <- progress_step_increment * 4

        sgrna_info_for_perm <- data_sg_ranked %>%
            select(!!sym(gene_col), "diff_score") %>%
            filter(!is.na(diff_score))
        if (nrow(sgrna_info_for_perm) == 0) { # Added check for empty sgrna_info_for_perm
            progress_msg("No sgRNAs with valid diff_scores for P-value permutations. Skipping.", val_increment = perm_loop_total_increment) # Use up the allocated progress
            gene_summary_df$P_positive <- NA_real_
            gene_summary_df$P_negative <- NA_real_
        } else {
            unique_gene_names_from_summary <- unique(gene_summary_df[[gene_col]][!is.na(gene_summary_df[[gene_col]])])
            if (length(unique_gene_names_from_summary) == 0) {
                progress_msg("No unique genes in summary for P-value permutations. Skipping.", val_increment = perm_loop_total_increment)
                gene_summary_df$P_positive <- NA_real_
                gene_summary_df$P_negative <- NA_real_
            } else {
                # Select Engine
                data_size <- nrow(sgrna_info_for_perm)
                n_genes <- length(unique_gene_names_from_summary)

                algorithm_choice <- user_engine_choice
                engine_name <- switch(algorithm_choice,
                    "cpp" = "C++ Ultra-Fast Engine",
                    "r_parallel" = "R Parallel Engine",
                    "Engine"
                )
                progress_msg(paste("Using", engine_name), val_increment = progress_step_increment * 0.05)

                # Dataset size warning
                if (n_genes > 10000) {
                    progress_msg(paste("Large dataset (", n_genes, " genes) - Using R Parallel Engine"), val_increment = progress_step_increment * 0.05)
                }

                progress_msg(paste("Data size:", data_size, "sgRNA,", n_genes, "genes;", N_perm, "permutations"), val_increment = progress_step_increment * 0.1)

                use_high_performance <- TRUE

                if (use_high_performance) {
                    # Algorithm selection message
                    algorithm_msg <- "R Parallel Algorithm Enabled (Estimated 3-15x speedup)"
                    progress_msg(algorithm_msg, val_increment = progress_step_increment * 0.1)

                    # System Resource Check
                    n_cores_available <- parallel::detectCores()
                    n_cores_to_use <- min(n_cores_available - 1, 8)
                    n_cores_to_use <- max(n_cores_to_use, 1)

                    progress_msg(paste("System detected", n_cores_available, "CPU cores, using", n_cores_to_use, "cores"), val_increment = progress_step_increment * 0.1)

                    # Progress Tracking Callback
                    high_perf_progress <- function(msg) {
                        # Parse progress message, extract info
                        if (!is.null(shiny_session)) {
                            # Check if it is a progress update message
                            if (grepl("Permutation Test:", msg) && grepl("\\(.*%\\)", msg)) {
                                # Extract percentage
                                percent_match <- regmatches(msg, regexpr("\\d+\\.?\\d*%", msg))
                                if (length(percent_match) > 0) {
                                    percent_value <- as.numeric(gsub("%", "", percent_match[1])) / 100

                                    # Calculate absolute progress
                                    current_base_progress <- func_current_progress_abs
                                    total_perm_progress <- progress_step_increment * 3.6
                                    adjusted_progress <- current_base_progress + (percent_value * total_perm_progress)

                                    # Update Shiny Progress
                                    shiny::setProgress(
                                        value = min(adjusted_progress, 0.95),
                                        message = "Permutation test in progress...",
                                        detail = msg,
                                        session = shiny_session
                                    )
                                }
                            } else {
                                # Console only
                                cat(paste("[INFO]", msg, "\n"))
                            }
                        } else {
                            cat(paste(msg, "\n"))
                        }
                    }

                    start_time <- Sys.time()
                    permutation_results <- tryCatch(
                        {
                            # Use R parallel engine

                            # Expose session for C++ logic
                            if (!is.null(shiny_session)) {
                                assign(".shiny_session", shiny_session, envir = .GlobalEnv)
                            }

                            # Use engine selector
                            engine_selector(
                                sgrna_info_dt = sgrna_info_for_perm,
                                gene_summary_dt = gene_summary_df,
                                gene_col = gene_col,
                                N_perm = N_perm,
                                n_cores = n_cores_to_use,
                                progress_callback = high_perf_progress,
                                user_engine_choice = algorithm_choice, # Pass user selection
                                min_sgrna_threshold = min_sgrna_threshold # Pass threshold
                            )
                        },
                        error = function(e) {
                            progress_msg(paste("Algorithm failed, falling back to standard method:", e$message))
                            NULL
                        }
                    )

                    if (!is.null(permutation_results)) {
                        elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
                        gene_summary_df$P_positive <- permutation_results$P_positive
                        gene_summary_df$P_negative <- permutation_results$P_negative

                        gene_summary_df$P_positive <- permutation_results$P_positive
                        gene_summary_df$P_negative <- permutation_results$P_negative

                        # Performance report
                        speed_msg <- paste("Permutation test completed! Total time:", round(elapsed_time, 2), "seconds")
                        if (N_perm >= 1000) {
                            estimated_standard_time <- N_perm * data_size * 0.0001 # Estimate standard method time
                            if (estimated_standard_time > elapsed_time) {
                                speedup <- round(estimated_standard_time / elapsed_time, 1)
                                speed_msg <- paste(speed_msg, "- Estimated speedup", speedup, "x")
                            }
                        }
                        progress_msg(speed_msg, val_increment = progress_step_increment * 3.6)
                    } else {
                        use_high_performance <- FALSE # Fallback
                    }
                }

                if (!use_high_performance) {
                    # Legacy Method
                    progress_msg("Using standard permutation test method...", val_increment = progress_step_increment * 0.2)

                    # Initialize Tracker
                    standard_global_tracker <- create_global_permutation_tracker(
                        total_permutations = N_perm,
                        shiny_session = shiny_session
                    )
                    standard_global_tracker(0, "Starting standard permutation test...", force_update = TRUE)

                    # Throttling updates
                    update_every_n_perms <- if (N_perm <= 100) {
                        5 # <100: Update every 5
                    } else if (N_perm <= 1000) {
                        max(20, round(N_perm / 20)) # 100-1000: Update 20 times
                    } else {
                        max(50, round(N_perm / 30)) # >1000: Update 30 times
                    }

                    perm_single_iter_increment <- (perm_loop_total_increment * 0.8) / N_perm

                    perm_pos_matrix <- matrix(NA_real_, nrow = length(unique_gene_names_from_summary), ncol = N_perm, dimnames = list(unique_gene_names_from_summary, NULL))
                    perm_neg_matrix <- matrix(NA_real_, nrow = length(unique_gene_names_from_summary), ncol = N_perm, dimnames = list(unique_gene_names_from_summary, NULL))

                    start_standard_time <- Sys.time()

                    for (p_idx in 1:N_perm) {
                        perm_sgrna_data <- sgrna_info_for_perm %>% mutate(shuffled_Gene_col_values = sample(sgrna_info_for_perm[[gene_col]]))

                        # Adaptive Top-N Permutation Test: Positive Scores
                        perm_pos_iter <- perm_sgrna_data %>%
                            group_by(shuffled_Gene_col_values) %>%
                            filter(n() >= min_sgrna_threshold) %>% # Only consider genes with >=min_sgrna_threshold sgRNAs
                            arrange(desc(diff_score)) %>%
                            mutate(
                                k = ceiling(2 * n() / 3)
                            ) %>%
                            filter(row_number() <= k) %>%
                            summarise(GPS = mean(diff_score, na.rm = TRUE), .groups = "drop") %>%
                            rename(PermGene = shuffled_Gene_col_values)

                        # Adaptive Top-N Permutation Test: Negative Scores
                        perm_neg_iter <- perm_sgrna_data %>%
                            group_by(shuffled_Gene_col_values) %>%
                            filter(n() >= min_sgrna_threshold) %>% # Only consider genes with >=min_sgrna_threshold sgRNAs
                            arrange(diff_score) %>%
                            mutate(
                                k = ceiling(2 * n() / 3)
                            ) %>%
                            filter(row_number() <= k) %>%
                            summarise(GNS = mean(diff_score, na.rm = TRUE), .groups = "drop") %>%
                            rename(PermGene = shuffled_Gene_col_values)
                        common_pos <- intersect(rownames(perm_pos_matrix), perm_pos_iter$PermGene)
                        if (length(common_pos) > 0) perm_pos_matrix[common_pos, p_idx] <- perm_pos_iter$GPS[match(common_pos, perm_pos_iter$PermGene)]
                        common_neg <- intersect(rownames(perm_neg_matrix), perm_neg_iter$PermGene)
                        if (length(common_neg) > 0) perm_neg_matrix[common_neg, p_idx] <- perm_neg_iter$GNS[match(common_neg, perm_neg_iter$PermGene)]

                        # Update global progress
                        standard_global_tracker(1) # Update for each permutation

                        # Update Shiny progress (less frequent)
                        if (p_idx %% update_every_n_perms == 0 || p_idx == 1 || p_idx == N_perm) {
                            progress_increment <- perm_single_iter_increment * update_every_n_perms
                            if (p_idx == N_perm) progress_increment <- perm_single_iter_increment * (N_perm %% update_every_n_perms)
                            progress_msg("", val_increment = progress_increment) # Silent update
                        }
                    }

                    # Standard method completed
                    standard_elapsed <- as.numeric(Sys.time() - start_standard_time, units = "secs")
                    if (standard_elapsed > 1) {
                        speed_standard <- round(N_perm / standard_elapsed, 1)
                        progress_msg(paste("Standard permutation completed! Time:", round(standard_elapsed, 2), "seconds"))
                    }
                    gene_summary_df$P_positive <- NA_real_
                    gene_summary_df$P_negative <- NA_real_
                    for (i in 1:nrow(gene_summary_df)) {
                        cur_gene <- gene_summary_df[[gene_col]][i]
                        if (is.na(cur_gene) || !cur_gene %in% unique_gene_names_from_summary) next
                        obs_pos_score <- gene_summary_df$GenePositiveScore[i]
                        obs_neg_score <- gene_summary_df$GeneNegativeScore[i]
                        if (!is.na(obs_pos_score)) {
                            null_pos <- perm_pos_matrix[cur_gene, ]
                            valid_perms <- sum(!is.na(null_pos))
                            if (valid_perms > 0) gene_summary_df$P_positive[i] <- (sum(null_pos >= obs_pos_score, na.rm = TRUE) + 1) / (valid_perms + 1)
                        }
                        if (!is.na(obs_neg_score)) {
                            null_neg <- perm_neg_matrix[cur_gene, ]
                            valid_perms <- sum(!is.na(null_neg))
                            if (valid_perms > 0) gene_summary_df$P_negative[i] <- (sum(null_neg <= obs_neg_score, na.rm = TRUE) + 1) / (valid_perms + 1)
                        }
                    }
                }
            } # End else for unique_gene_names_from_summary check
        } # End else for sgrna_info_for_perm check

        progress_msg("Step V: P-value calculation completed.", val_increment = progress_step_increment)
        progress_msg("Step VI: Applying BH correction...", val_increment = progress_step_increment)
        gene_summary_df$P_positive_adj_bh <- NA_real_
        gene_summary_df$P_negative_adj_bh <- NA_real_
        if ("P_positive" %in% names(gene_summary_df) && any(!is.na(gene_summary_df$P_positive))) {
            valid_idx <- !is.na(gene_summary_df$P_positive)
            gene_summary_df$P_positive_adj_bh[valid_idx] <- p.adjust(gene_summary_df$P_positive[valid_idx], method = "BH")
        }
        if ("P_negative" %in% names(gene_summary_df) && any(!is.na(gene_summary_df$P_negative))) {
            valid_idx <- !is.na(gene_summary_df$P_negative)
            gene_summary_df$P_negative_adj_bh[valid_idx] <- p.adjust(gene_summary_df$P_negative[valid_idx], method = "BH")
        }
        progress_msg("Step VI: BH correction completed.", val_increment = progress_step_increment)
    } else {
        progress_msg("N_perm < 10, P-value calculation and BH correction skipped.", val_increment = progress_step_increment * 2)
        gene_summary_df$P_positive <- NA_real_
        gene_summary_df$P_negative <- NA_real_
        gene_summary_df$P_positive_adj_bh <- NA_real_
        gene_summary_df$P_negative_adj_bh <- NA_real_
    }


    return(list(
        processed_sg_data = data_sg_ranked,
        gene_summary_data = gene_summary_df,
        normalized_counts = data_normalized,
        filtered_genes_sgrna_data = filtered_sgrna_data,
        params = list( # Store used parameters
            gRNA_col = gRNA_col,
            gene_col = gene_col,
            sequence_col = sequence_col
        )
    ))
}
