# SPDX-License-Identifier: GPL-3.0-or-later
# permutation_functions.R
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


# permutation_functions.R
# Permutation Test Logic

# --- Permutation Test Selector ---
#' Permutation Test Selector
select_permutation_engine <- function(
    sgrna_info_dt,
    gene_summary_dt,
    gene_col,
    N_perm = 1000,
    n_cores = NULL,
    progress_callback = NULL,
    user_engine_choice = "cpp", # Default to C++ engine
    min_sgrna_threshold = 3 # Add threshold parameter
    ) {
    # Engine Selection
    engine_choice <- user_engine_choice

    # Validate Engine
    if (!engine_choice %in% c("cpp", "r_parallel")) {
        engine_choice <- "cpp" # Default fallback to C++
        if (!is.null(progress_callback)) {
            progress_callback("Unknown engine choice, using C++ engine")
        }
    }

    engine_name <- switch(engine_choice,
        "cpp" = "C++ Ultra-Fast Engine",
        "r_parallel" = "R Parallel Engine",
        "Engine"
    )

    if (!is.null(progress_callback)) {
        progress_callback(paste("User Selection:", engine_name))
    }

    # Execute Engine
    tryCatch(
        {
            perform_permutation(
                sgrna_info_dt = sgrna_info_dt,
                gene_summary_dt = gene_summary_dt,
                gene_col = gene_col,
                N_perm = N_perm,
                n_cores = n_cores,
                use_data_table = TRUE,
                progress_callback = progress_callback,
                min_sgrna_threshold = min_sgrna_threshold # Pass threshold
            )
        },
        error = function(e) {
            if (!is.null(progress_callback)) {
                progress_callback(paste("Engine", engine_choice, "failed, falling back to R Parallel version:", e$message))
            }
            perform_permutation(
                sgrna_info_dt = sgrna_info_dt,
                gene_summary_dt = gene_summary_dt,
                gene_col = gene_col,
                N_perm = N_perm,
                n_cores = n_cores,
                use_data_table = TRUE,
                progress_callback = progress_callback,
                min_sgrna_threshold = min_sgrna_threshold # Pass threshold
            )
        }
    )
}

# --- Permutation Test Function (New) ---
#' Permutation Test - Parallel Version
#'
#' @param sgrna_info_dt data.table containing gene and diff_score columns
#' @param gene_summary_dt data.table containing gene summary data
#' @param gene_col column name for genes
#' @param N_perm number of permutations
#' @param n_cores number of CPU cores to use (NULL for auto-detect)
#' @param use_data_table logical, whether to use data.table for speed
#' @param progress_callback function to report progress
#' @param min_sgrna_threshold minimum number of sgRNAs per gene for scoring
perform_permutation <- function(
    sgrna_info_dt,
    gene_summary_dt,
    gene_col,
    N_perm = 1000,
    n_cores = NULL,
    use_data_table = TRUE,
    progress_callback = NULL,
    min_sgrna_threshold = 3) {
    # Auto-detect CPU Cores
    if (is.null(n_cores)) {
        n_cores <- min(parallel::detectCores() - 1, 8) # Keep 1 core, use max 8
        n_cores <- max(n_cores, 1) # Use at least 1 core
    }

    if (!is.null(progress_callback)) {
        progress_callback(paste("Using", n_cores, "CPU cores for parallel calculation"))
    }

    # Convert to data.table
    if (use_data_table) {
        if (!inherits(sgrna_info_dt, "data.table")) {
            sgrna_info_dt <- data.table::as.data.table(sgrna_info_dt)
        }
        if (!inherits(gene_summary_dt, "data.table")) {
            gene_summary_dt <- data.table::as.data.table(gene_summary_dt)
        }
        data.table::setkeyv(sgrna_info_dt, gene_col)
    }

    unique_genes <- unique(gene_summary_dt[[gene_col]][!is.na(gene_summary_dt[[gene_col]])])
    n_genes <- length(unique_genes)

    if (n_genes == 0) {
        if (!is.null(progress_callback)) {
            progress_callback("No valid genes for permutation test")
        }
        return(list(
            P_positive = rep(NA_real_, nrow(gene_summary_dt)),
            P_negative = rep(NA_real_, nrow(gene_summary_dt))
        ))
    }

    # Setup Parallel Backend
    cl <- NULL
    if (n_cores > 1) {
        if (!is.null(progress_callback)) {
            progress_callback("Setting up parallel calculation environment...")
        }

        # Estimate Memory for Large Data
        data_size_mb <- object.size(sgrna_info_dt) / 1024^2 # Convert to MB
        estimated_need_gb <- max(0.5, min(8, data_size_mb * n_cores / 500)) # Estimate needed GB

        # Lightweight Parallel Strategy
        if (n_genes > 10000) {
            if (!is.null(progress_callback)) {
                progress_callback("Large dataset detected: Enabling ultra-lightweight parallel mode")
            }
            # Use parallel package instead of future (avoid complex global object transfer)
            cl <- parallel::makeCluster(n_cores)
            on.exit(parallel::stopCluster(cl), add = TRUE)

            # Only transfer basic packages to worker nodes
            parallel::clusterEvalQ(cl, {
                library(data.table)
            })
        } else {
            # Standard datasets also use parallel package
            cl <- parallel::makeCluster(n_cores)
            on.exit(parallel::stopCluster(cl), add = TRUE)

            # Transfer necessary packages to worker nodes
            parallel::clusterEvalQ(cl, {
                library(data.table)
            })
        }
    }

    # Dynamic Batch Sizing
    if (n_genes > 15000) {
        batch_size <- min(25, N_perm) # Huge dataset: small batch
    } else if (n_genes > 5000) {
        batch_size <- min(50, N_perm) # Large dataset: medium batch
    } else {
        batch_size <- min(100, N_perm) # Normal size: standard batch
    }
    n_batches <- ceiling(N_perm / batch_size)

    if (!is.null(progress_callback)) {
        progress_callback(paste("Using batch size:", batch_size, "- Total", n_batches, "batches"))
    }

    # Initialize result matrix
    perm_pos_matrix <- matrix(NA_real_,
        nrow = n_genes, ncol = N_perm,
        dimnames = list(unique_genes, NULL)
    )
    perm_neg_matrix <- matrix(NA_real_,
        nrow = n_genes, ncol = N_perm,
        dimnames = list(unique_genes, NULL)
    )

    # Init Global Tracker
    global_tracker <- NULL
    if (!is.null(progress_callback)) {
        progress_callback(paste("Initializing permutation matrix -", n_genes, "genes ×", N_perm, "permutations"))
        global_tracker <- create_global_permutation_tracker(
            total_permutations = N_perm,
            shiny_session = if (!is.null(progress_callback)) getDefaultReactiveDomain() else NULL
        )
        # Initialize display
        global_tracker(0, "Preparing to start permutation test...", force_update = TRUE)
    }

    if (use_data_table) {
        # data.table Version
        gene_sgrna_dt <- sgrna_info_dt[, c(gene_col, "diff_score"), with = FALSE]
        data.table::setnames(gene_sgrna_dt, gene_col, "gene") # Rename for consistency
        gene_sgrna_dt <- gene_sgrna_dt[!is.na(diff_score)]

        for (batch_idx in 1:n_batches) {
            start_perm <- (batch_idx - 1) * batch_size + 1
            end_perm <- min(batch_idx * batch_size, N_perm)
            current_batch_size <- end_perm - start_perm + 1

            # Parallel calculation - choose strategy based on dataset size
            if (n_cores > 1) {
                gene_data <- gene_sgrna_dt$gene
                score_data <- gene_sgrna_dt$diff_score
                unique_genes_vec <- unique_genes

                if (n_genes > 10000) {
                    # Large dataset: use parallel package (ultra-lightweight)
                    # Pre-transfer data to cluster nodes
                    parallel::clusterExport(cl, c("gene_data", "score_data", "unique_genes_vec"),
                        envir = environment()
                    )

                    batch_results <- parallel::parLapply(cl, 1:current_batch_size, function(i) {
                        # Lightweight single permutation function
                        shuffled_genes <- sample(gene_data)

                        temp_dt <- data.table::data.table(
                            gene = gene_data,
                            diff_score = score_data,
                            shuffled_gene = shuffled_genes
                        )
                        data.table::setkey(temp_dt, shuffled_gene)

                        # Adaptive Top-N calculation: Large dataset version
                        pos_scores_dt <- temp_dt[,
                            {
                                n_sgrnas <- .N
                                n_sgrnas <- .N
                                k <- ceiling(2 * n_sgrnas / 3)
                                .(pos_score = mean(head(sort(diff_score, decreasing = TRUE), k), na.rm = TRUE))
                            },
                            by = shuffled_gene
                        ][shuffled_gene %in% unique_genes_vec]

                        neg_scores_dt <- temp_dt[,
                            {
                                n_sgrnas <- .N
                                n_sgrnas <- .N
                                k <- ceiling(2 * n_sgrnas / 3)
                                .(neg_score = mean(head(sort(diff_score, decreasing = FALSE), k), na.rm = TRUE))
                            },
                            by = shuffled_gene
                        ][shuffled_gene %in% unique_genes_vec]

                        pos_scores <- pos_scores_dt$pos_score
                        names(pos_scores) <- pos_scores_dt$shuffled_gene

                        neg_scores <- neg_scores_dt$neg_score
                        names(neg_scores) <- neg_scores_dt$shuffled_gene

                        list(pos_scores = pos_scores, neg_scores = neg_scores)
                    })
                } else {
                    # Standard dataset: use parallel package
                    # Pre-transfer data to cluster nodes
                    parallel::clusterExport(cl, c("gene_data", "score_data", "unique_genes_vec"),
                        envir = environment()
                    )

                    batch_results <- parallel::parLapply(cl, 1:current_batch_size, function(i) {
                        shuffled_genes <- sample(gene_data)

                        temp_dt <- data.table::data.table(
                            gene = gene_data,
                            diff_score = score_data,
                            shuffled_gene = shuffled_genes
                        )
                        data.table::setkey(temp_dt, shuffled_gene)

                        # Adaptive Top-N calculation: Standard dataset version
                        pos_scores_dt <- temp_dt[,
                            {
                                n_sgrnas <- .N
                                n_sgrnas <- .N
                                k <- ceiling(2 * n_sgrnas / 3)
                                .(pos_score = mean(head(sort(diff_score, decreasing = TRUE), k), na.rm = TRUE))
                            },
                            by = shuffled_gene
                        ][shuffled_gene %in% unique_genes_vec]

                        neg_scores_dt <- temp_dt[,
                            {
                                n_sgrnas <- .N
                                n_sgrnas <- .N
                                k <- ceiling(2 * n_sgrnas / 3)
                                .(neg_score = mean(head(sort(diff_score, decreasing = FALSE), k), na.rm = TRUE))
                            },
                            by = shuffled_gene
                        ][shuffled_gene %in% unique_genes_vec]

                        pos_scores <- pos_scores_dt$pos_score
                        names(pos_scores) <- pos_scores_dt$shuffled_gene

                        neg_scores <- neg_scores_dt$neg_score
                        names(neg_scores) <- neg_scores_dt$shuffled_gene

                        list(pos_scores = pos_scores, neg_scores = neg_scores)
                    })
                }
            } else {
                batch_results <- lapply(1:current_batch_size, function(i) {
                    perform_single_permutation_dt(gene_sgrna_dt, unique_genes)
                })
            }

            # Store results and update global progress
            for (i in 1:current_batch_size) {
                perm_idx <- start_perm + i - 1
                result <- batch_results[[i]]
                if (!is.null(result$pos_scores)) {
                    common_genes <- intersect(rownames(perm_pos_matrix), names(result$pos_scores))
                    if (length(common_genes) > 0) {
                        perm_pos_matrix[common_genes, perm_idx] <- result$pos_scores[common_genes]
                    }
                }
                if (!is.null(result$neg_scores)) {
                    common_genes <- intersect(rownames(perm_neg_matrix), names(result$neg_scores))
                    if (length(common_genes) > 0) {
                        perm_neg_matrix[common_genes, perm_idx] <- result$neg_scores[common_genes]
                    }
                }

                # 更新全局进度
                if (!is.null(global_tracker)) {
                    global_tracker(1) # 增加1个置换
                }
            }
        }
    } else {
        # Standard dplyr version (compatibility)
        for (batch_idx in 1:n_batches) {
            start_perm <- (batch_idx - 1) * batch_size + 1
            end_perm <- min(batch_idx * batch_size, N_perm)
            current_batch_size <- end_perm - start_perm + 1

            # Parallel calculation - dplyr version
            if (n_cores > 1) {
                gene_data <- sgrna_info_dt[[gene_col]]
                score_data <- sgrna_info_dt$diff_score
                gene_col_name <- gene_col

                if (n_genes > 10000) {
                    # Large dataset: use parallel package (ultra-lightweight)
                    parallel::clusterExport(cl, c("gene_data", "score_data", "gene_col_name"),
                        envir = environment()
                    )

                    # Load dplyr to worker nodes
                    parallel::clusterEvalQ(cl, {
                        library(dplyr)
                    })

                    batch_results <- parallel::parLapply(cl, 1:current_batch_size, function(i) {
                        shuffled_genes <- sample(gene_data)

                        temp_df <- data.frame(
                            gene = gene_data,
                            diff_score = score_data,
                            shuffled_gene = shuffled_genes,
                            stringsAsFactors = FALSE
                        )

                        perm_pos_iter <- temp_df %>%
                            dplyr::group_by(shuffled_gene) %>%
                            dplyr::filter(dplyr::n() >= 1) %>%
                            dplyr::arrange(dplyr::desc(diff_score)) %>%
                            dplyr::mutate(
                                n_sgrnas = dplyr::n(),
                                k = ceiling(2 * n_sgrnas / 3)
                            ) %>%
                            dplyr::slice_head(n = k[1]) %>%
                            dplyr::summarise(GPS = mean(diff_score, na.rm = TRUE), .groups = "drop") %>%
                            dplyr::rename(PermGene = shuffled_gene)

                        perm_neg_iter <- temp_df %>%
                            dplyr::group_by(shuffled_gene) %>%
                            dplyr::filter(dplyr::n() >= 1) %>%
                            dplyr::arrange(diff_score) %>%
                            dplyr::mutate(
                                n_sgrnas = dplyr::n(),
                                k = ceiling(2 * n_sgrnas / 3)
                            ) %>%
                            dplyr::slice_head(n = k[1]) %>%
                            dplyr::summarise(GNS = mean(diff_score, na.rm = TRUE), .groups = "drop") %>%
                            dplyr::rename(PermGene = shuffled_gene)

                        list(pos_scores = perm_pos_iter, neg_scores = perm_neg_iter)
                    })
                } else {
                    # Standard dataset: use parallel package
                    # Pre-transfer data to cluster nodes
                    parallel::clusterExport(cl, c("gene_data", "score_data", "gene_col_name"),
                        envir = environment()
                    )

                    # Load dplyr to worker nodes
                    parallel::clusterEvalQ(cl, {
                        library(dplyr)
                    })

                    batch_results <- parallel::parLapply(cl, 1:current_batch_size, function(i) {
                        shuffled_genes <- sample(gene_data)

                        temp_df <- data.frame(
                            gene = gene_data,
                            diff_score = score_data,
                            shuffled_gene = shuffled_genes,
                            stringsAsFactors = FALSE
                        )

                        perm_pos_iter <- temp_df %>%
                            dplyr::group_by(shuffled_gene) %>%
                            dplyr::filter(dplyr::n() >= 1) %>%
                            dplyr::arrange(dplyr::desc(diff_score)) %>%
                            dplyr::mutate(
                                n_sgrnas = dplyr::n(),
                                k = ceiling(2 * n_sgrnas / 3)
                            ) %>%
                            dplyr::slice_head(n = k[1]) %>%
                            dplyr::summarise(GPS = mean(diff_score, na.rm = TRUE), .groups = "drop") %>%
                            dplyr::rename(PermGene = shuffled_gene)

                        perm_neg_iter <- temp_df %>%
                            dplyr::group_by(shuffled_gene) %>%
                            dplyr::filter(dplyr::n() >= 1) %>%
                            dplyr::arrange(diff_score) %>%
                            dplyr::mutate(
                                n_sgrnas = dplyr::n(),
                                k = ceiling(2 * n_sgrnas / 3)
                            ) %>%
                            dplyr::slice_head(n = k[1]) %>%
                            dplyr::summarise(GNS = mean(diff_score, na.rm = TRUE), .groups = "drop") %>%
                            dplyr::rename(PermGene = shuffled_gene)

                        list(pos_scores = perm_pos_iter, neg_scores = perm_neg_iter)
                    })
                }
            } else {
                batch_results <- lapply(1:current_batch_size, function(i) {
                    perform_single_permutation_dplyr(sgrna_info_dt, gene_col)
                })
            }

            # Store results and update global progress
            for (i in 1:current_batch_size) {
                perm_idx <- start_perm + i - 1
                result <- batch_results[[i]]
                if (!is.null(result$pos_scores)) {
                    common_genes <- intersect(rownames(perm_pos_matrix), result$pos_scores$PermGene)
                    if (length(common_genes) > 0) {
                        perm_pos_matrix[common_genes, perm_idx] <- result$pos_scores$GPS[match(common_genes, result$pos_scores$PermGene)]
                    }
                }
                if (!is.null(result$neg_scores)) {
                    common_genes <- intersect(rownames(perm_neg_matrix), result$neg_scores$PermGene)
                    if (length(common_genes) > 0) {
                        perm_neg_matrix[common_genes, perm_idx] <- result$neg_scores$GNS[match(common_genes, result$neg_scores$PermGene)]
                    }
                }

                # 更新全局进度
                if (!is.null(global_tracker)) {
                    global_tracker(1) # 增加1个置换
                }
            }
        }
    }



    if (!is.null(progress_callback)) {
        progress_callback("Starting p-value calculation...")
        if (!is.null(global_tracker)) {
            global_tracker(0, "Calculating p-values...", force_update = TRUE)
        }
    }

    # Vectorized p-value calculation
    P_positive <- rep(NA_real_, nrow(gene_summary_dt))
    P_negative <- rep(NA_real_, nrow(gene_summary_dt))

    for (i in 1:nrow(gene_summary_dt)) {
        cur_gene <- gene_summary_dt[[gene_col]][i]
        if (is.na(cur_gene) || !cur_gene %in% unique_genes) next

        obs_pos_score <- gene_summary_dt$GenePositiveScore[i]
        obs_neg_score <- gene_summary_dt$GeneNegativeScore[i]

        if (!is.na(obs_pos_score)) {
            null_pos <- perm_pos_matrix[cur_gene, ]
            valid_perms <- sum(!is.na(null_pos))
            if (valid_perms > 0) {
                P_positive[i] <- (sum(null_pos >= obs_pos_score, na.rm = TRUE) + 1) / (valid_perms + 1)
            }
        }

        if (!is.na(obs_neg_score)) {
            null_neg <- perm_neg_matrix[cur_gene, ]
            valid_perms <- sum(!is.na(null_neg))
            if (valid_perms > 0) {
                P_negative[i] <- (sum(null_neg <= obs_neg_score, na.rm = TRUE) + 1) / (valid_perms + 1)
            }
        }
    }

    return(list(P_positive = P_positive, P_negative = P_negative))
}

#' Single Permutation - data.table Version
perform_single_permutation_dt <- function(gene_sgrna_dt, unique_genes) {
    # Efficiently shuffle gene labels
    shuffled_genes <- sample(gene_sgrna_dt$gene)

    # Use data.table for efficient grouping operations
    perm_dt <- data.table::copy(gene_sgrna_dt)
    perm_dt[, shuffled_gene := shuffled_genes]
    data.table::setkey(perm_dt, shuffled_gene)

    # Adaptive Top-N calculation for positive and negative scores
    pos_scores_dt <- perm_dt[,
        {
            n_sgrnas <- .N
            k <- ceiling(2 * n_sgrnas / 3)
            .(pos_score = mean(head(sort(diff_score, decreasing = TRUE), k), na.rm = TRUE))
        },
        by = shuffled_gene
    ][shuffled_gene %in% unique_genes]

    neg_scores_dt <- perm_dt[,
        {
            n_sgrnas <- .N
            k <- ceiling(2 * n_sgrnas / 3)
            .(neg_score = mean(head(sort(diff_score, decreasing = FALSE), k), na.rm = TRUE))
        },
        by = shuffled_gene
    ][shuffled_gene %in% unique_genes]

    # Convert to named vector
    pos_scores <- pos_scores_dt$pos_score
    names(pos_scores) <- pos_scores_dt$shuffled_gene

    neg_scores <- neg_scores_dt$neg_score
    names(neg_scores) <- neg_scores_dt$shuffled_gene

    return(list(pos_scores = pos_scores, neg_scores = neg_scores))
}

#' Single Permutation - dplyr Version (Compatibility)
perform_single_permutation_dplyr <- function(sgrna_info_df, gene_col) {
    perm_sgrna_data <- sgrna_info_df %>%
        mutate(shuffled_Gene_col_values = sample(.data[[gene_col]]))

    perm_pos_iter <- perm_sgrna_data %>%
        group_by(shuffled_Gene_col_values) %>%
        filter(n() >= 1) %>%
        arrange(desc(diff_score)) %>%
        mutate(
            n_sgrnas = n(),
            k = ceiling(2 * n_sgrnas / 3)
        ) %>%
        slice_head(n = k[1]) %>%
        summarise(GPS = mean(diff_score, na.rm = TRUE), .groups = "drop") %>%
        rename(PermGene = shuffled_Gene_col_values)

    perm_neg_iter <- perm_sgrna_data %>%
        group_by(shuffled_Gene_col_values) %>%
        filter(n() >= 1) %>%
        arrange(diff_score) %>%
        mutate(
            n_sgrnas = n(),
            k = ceiling(2 * n_sgrnas / 3)
        ) %>%
        slice_head(n = k[1]) %>%
        summarise(GNS = mean(diff_score, na.rm = TRUE), .groups = "drop") %>%
        rename(PermGene = shuffled_Gene_col_values)

    return(list(pos_scores = perm_pos_iter, neg_scores = perm_neg_iter))
}
