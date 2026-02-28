# SPDX-License-Identifier: GPL-3.0-or-later
# cpp_progress_wrapper.R
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


# C++ Progress Tracker Wrapper
# Converts C++ progress output to Shiny progress updates

#' Create C++ Progress Tracker Wrapper
create_cpp_progress_wrapper <- function(n_permutations, progress_callback = NULL, shiny_session = NULL) {
    if (is.null(progress_callback) && is.null(shiny_session)) {
        return(NULL) # No progress tracking needed
    }

    # Create global progress tracker
    if (!is.null(shiny_session)) {
        global_tracker <- create_global_permutation_tracker(n_permutations, shiny_session)
    } else {
        global_tracker <- NULL
    }

    # Function to monitor C++ output
    monitor_cpp_progress <- function() {
        if (!is.null(global_tracker)) {
            # Simulate progress update
            for (i in 1:10) {
                Sys.sleep(0.1)
                global_tracker(n_permutations / 10, paste("C++ Permutation test in progress..."))
            }
        }

        if (!is.null(progress_callback)) {
            progress_callback("C++ Permutation test in progress...")
        }
    }

    return(monitor_cpp_progress)
}

#' C++ Permutation Test with Progress Tracking
perform_cpp_permutation_with_progress <- function(
    sgrna_data,
    gene_summary_data = NULL,
    gene_col = "Gene",
    score_col = "diff_score",
    n_permutations = 1000,
    min_sgrna_threshold = 3,
    n_threads = 0,
    progress_callback = NULL,
    shiny_session = NULL) {
    if (!is.null(progress_callback)) {
        progress_callback("Starting C++ Permutation Engine...")
    }

    # Data preprocessing (reuse existing logic)
    if (!is.null(progress_callback)) {
        progress_callback("Preprocessing data...")
    }

    # Extract necessary columns
    if (!gene_col %in% colnames(sgrna_data)) {
        stop("Gene column not found: ", gene_col)
    }
    if (!score_col %in% colnames(sgrna_data)) {
        stop("Score column not found: ", score_col)
    }

    # Filter valid data
    valid_rows <- !is.na(sgrna_data[[gene_col]]) &
        !is.na(sgrna_data[[score_col]]) &
        is.finite(sgrna_data[[score_col]])

    if (sum(valid_rows) == 0) {
        stop("No valid data rows")
    }

    clean_data <- sgrna_data[valid_rows, ]

    # Filter genes with sgRNA count >= min_sgrna_threshold
    gene_counts <- table(clean_data[[gene_col]])
    valid_genes <- names(gene_counts)[gene_counts >= min_sgrna_threshold]

    if (length(valid_genes) == 0) {
        stop(paste0("No genes meet the condition (sgRNA count >= ", min_sgrna_threshold, ")"))
    }

    final_data <- clean_data[clean_data[[gene_col]] %in% valid_genes, ]

    if (!is.null(progress_callback)) {
        progress_callback(paste("Data preparation complete:", nrow(final_data), "sgRNAs,", length(valid_genes), "genes"))
    }

    # Prepare C++ input data
    diff_scores <- as.numeric(final_data[[score_col]])
    gene_labels <- as.character(final_data[[gene_col]])
    unique_genes <- sort(unique(gene_labels))

    # Create global progress tracker
    if (!is.null(shiny_session)) {
        global_tracker <- create_global_permutation_tracker(n_permutations, shiny_session)

        # Simplify notifications to avoid redundant popups
        if (!is.null(progress_callback)) {
            progress_callback("Executing permutation test...")
        }

        # Simulate batched execution of C++ calculation to provide progress updates
        batch_size <- max(100, n_permutations %/% 20) # Split into max 20 batches
        n_batches <- ceiling(n_permutations / batch_size)

        all_results <- list()
        cumulative_time <- 0

        for (batch_i in 1:n_batches) {
            start_perm <- (batch_i - 1) * batch_size + 1
            end_perm <- min(batch_i * batch_size, n_permutations)
            current_batch_size <- end_perm - start_perm + 1

            batch_start_time <- Sys.time()

            # Execute current batch
            seed <- sample(.Machine$integer.max, 1)
            batch_result <- tryCatch(
                {
                    perform_cpp_permutation_test(
                        diff_scores = diff_scores,
                        gene_labels = gene_labels,
                        unique_genes = unique_genes,
                        n_permutations = current_batch_size,
                        min_sgrna_threshold = min_sgrna_threshold,
                        seed = seed,
                        show_progress = FALSE # Do not show progress within batch
                    )
                },
                error = function(e) {
                    stop("C++ Permutation test failed: ", e$message)
                }
            )

            batch_end_time <- Sys.time()
            batch_elapsed <- as.numeric(difftime(batch_end_time, batch_start_time, units = "secs"))
            cumulative_time <- cumulative_time + batch_elapsed

            # Update global progress
            completed_perms <- end_perm
            progress_percent <- round((completed_perms / n_permutations) * 100)

            # Use global progress tracker only, no extra progress_callback
            global_tracker(
                current_batch_size,
                paste0("Permutation Test: ", batch_i, "/", n_batches, " Completed ", progress_percent, "%")
            )

            all_results[[batch_i]] <- batch_result
        }

        # Merge results (Aggregate P-values from all batches)
        first_res <- all_results[[1]]
        n_genes <- length(first_res$P_positive)

        sum_pos_counts <- numeric(n_genes)
        sum_neg_counts <- numeric(n_genes)
        valid_total_perms <- 0

        for (i in 1:n_batches) {
            start_perm <- (i - 1) * batch_size + 1
            end_perm <- min(i * batch_size, n_permutations)
            current_batch_size <- end_perm - start_perm + 1

            res <- all_results[[i]]

            # Reconstruct counts from P-values (P * N)
            # Use current_batch_size as weight
            sum_pos_counts <- sum_pos_counts + (res$P_positive * current_batch_size)
            sum_neg_counts <- sum_neg_counts + (res$P_negative * current_batch_size)
            valid_total_perms <- valid_total_perms + current_batch_size
        }

        final_result <- first_res
        final_result$P_positive <- sum_pos_counts / valid_total_perms
        final_result$P_negative <- sum_neg_counts / valid_total_perms
        final_result$n_permutations <- valid_total_perms
        final_result$elapsed_time <- cumulative_time
    } else {
        # Execute directly when no Shiny session
        if (!is.null(progress_callback)) {
            progress_callback(paste("Executing", n_permutations, "permutations..."))
        }

        start_time <- Sys.time()
        seed <- sample(.Machine$integer.max, 1)

        final_result <- tryCatch(
            {
                perform_cpp_permutation_test(
                    diff_scores = diff_scores,
                    gene_labels = gene_labels,
                    unique_genes = unique_genes,
                    n_permutations = n_permutations,
                    min_sgrna_threshold = min_sgrna_threshold,
                    seed = seed,
                    show_progress = TRUE
                )
            },
            error = function(e) {
                if (!is.null(progress_callback)) {
                    progress_callback(paste("C++ Calculation Failed:", e$message))
                }
                stop("C++ Permutation test failed: ", e$message)
            }
        )

        end_time <- Sys.time()
        final_result$elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    }

    # Final completion notification (Simplified)
    if (!is.null(progress_callback)) {
        progress_callback(paste("Permutation test completed, time:", round(final_result$elapsed_time, 2), "seconds"))
    }

    # Process return results (Reuse existing logic)
    if (!is.null(gene_summary_data)) {
        gene_summary_genes <- gene_summary_data[[gene_col]]
        n_summary_genes <- length(gene_summary_genes)

        matched_p_positive <- rep(NA_real_, n_summary_genes)
        matched_p_negative <- rep(NA_real_, n_summary_genes)

        result_map <- setNames(1:length(final_result$gene_names), final_result$gene_names)

        for (i in 1:n_summary_genes) {
            gene_name <- gene_summary_genes[i]
            if (!is.na(gene_name) && gene_name %in% names(result_map)) {
                result_idx <- result_map[gene_name]
                matched_p_positive[i] <- final_result$P_positive[result_idx]
                matched_p_negative[i] <- final_result$P_negative[result_idx]
            }
        }

        return(list(
            P_positive = matched_p_positive,
            P_negative = matched_p_negative,
            elapsed_time = final_result$elapsed_time,
            n_threads = 1,
            engine = "C++",
            performance_stats = list(
                n_valid_genes = final_result$n_valid_genes,
                n_permutations = final_result$n_permutations,
                speed = paste(round(final_result$n_permutations / final_result$elapsed_time), "perms/sec")
            )
        ))
    } else {
        return(list(
            P_positive = final_result$P_positive,
            P_negative = final_result$P_negative,
            gene_names = final_result$gene_names,
            elapsed_time = final_result$elapsed_time,
            n_threads = 1,
            engine = "C++",
            performance_stats = list(
                n_valid_genes = final_result$n_valid_genes,
                n_permutations = final_result$n_permutations,
                speed = paste(round(final_result$n_permutations / final_result$elapsed_time), "perms/sec")
            )
        ))
    }
}
