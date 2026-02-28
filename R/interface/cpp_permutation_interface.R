# SPDX-License-Identifier: GPL-3.0-or-later
# cpp_permutation_interface.R
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


# C++ Permutation Engine - R Interface
# Rcpp implementation

#' Check and install necessary packages
setup_cpp_environment <- function() {
    required_packages <- c("Rcpp", "RcppParallel")

    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            cat(paste("Installing", pkg, "package...\n"))
            tryCatch(
                {
                    install.packages(pkg, dependencies = TRUE)
                    cat(paste(pkg, "installed successfully\n"))
                },
                error = function(e) {
                    stop(paste("Failed to install", pkg, ":", e$message))
                }
            )
        }
    }

    # Load packages
    library(Rcpp)
    library(RcppParallel)

    cat("C++ environment setup complete\n")
    return(TRUE)
}

#' Compile C++ Engine
compile_cpp_engine <- function() {
    # Check if C++ source file exists
    cpp_file <- "src/cpp_permutation_engine.cpp"
    if (!file.exists(cpp_file)) {
        stop("C++ source file not found: ", cpp_file)
    }

    cat("Compiling C++ Permutation Engine...\n")

    tryCatch(
        {
            # Compile C++ code using Rcpp (suppress verbose output)
            suppressMessages({
                capture.output(
                    {
                        Rcpp::sourceCpp(cpp_file, verbose = FALSE, rebuild = TRUE)
                    },
                    type = "output"
                )
            })

            # Wait for compilation to complete
            Sys.sleep(1)

            # Verify if function is available
            if (exists("perform_cpp_permutation_test")) {
                cat("C++ Engine compiled successfully\n")

                return(TRUE)
            } else {
                stop("Function not available after compilation")
            }
        },
        error = function(e) {
            cat("C++ Engine compilation failed:", e$message, "\n")
            return(FALSE)
        }
    )
}

#' System Configuration Detection
detect_system_config <- function() {
    # Detect CPU cores
    n_cores <- parallel::detectCores()
    recommended_threads <- max(1, n_cores - 1)

    # Detect Memory
    memory_info <- tryCatch(
        {
            if (Sys.info()["sysname"] == "Darwin") {
                # macOS
                memory_bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
                memory_gb <- round(memory_bytes / 1024^3, 1)
            } else {
                # Default assumption: 8GB
                memory_gb <- 8
            }
            memory_gb
        },
        error = function(e) 8
    )

    cat("System Configuration:\n")
    cat(paste("  • CPU Cores:", n_cores, "\n"))
    cat(paste("  • Recommended Threads:", recommended_threads, "\n"))
    cat(paste("  • Memory:", memory_info, "GB\n"))

    return(list(
        n_cores = n_cores,
        optimal_threads = recommended_threads,
        memory_gb = memory_info
    ))
}

#' C++ Permutation Test Main Function
#' Identical to R algorithm but with significant performance improvement
perform_cpp_permutation <- function(
    sgrna_data,
    gene_summary_data = NULL,
    gene_col = "Gene",
    score_col = "diff_score",
    n_permutations = 1000,
    min_sgrna_threshold = 3,
    n_threads = 0,
    progress_callback = NULL) {
    if (!is.null(progress_callback)) {
        progress_callback("Starting C++ Permutation Engine...")
    }

    # Detect system configuration
    sys_config <- detect_system_config()

    if (n_threads <= 0) {
        n_threads <- sys_config$optimal_threads
    }

    if (!is.null(progress_callback)) {
        progress_callback(paste0("Using ", n_threads, " threads for parallel calculation"))
    }

    # Data Preprocessing
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

    if (!is.null(progress_callback)) {
        progress_callback(paste("Executing", n_permutations, "permutations..."))
    }

    # Execute C++ Permutation Test
    start_time <- Sys.time()

    # Generate random seed
    seed <- sample(.Machine$integer.max, 1)

    cpp_results <- tryCatch(
        {
            perform_cpp_permutation_test(
                diff_scores = diff_scores,
                gene_labels = gene_labels,
                unique_genes = unique_genes,
                n_permutations = n_permutations,
                min_sgrna_threshold = min_sgrna_threshold,
                seed = seed,
                show_progress = !is.null(progress_callback) # Show progress if callback provided
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
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    if (!is.null(progress_callback)) {
        progress_callback(paste("Permutation test completed, time:", round(elapsed_time, 2), "seconds"))
        progress_callback(paste("Performance:", cpp_results$n_valid_genes, "genes,", cpp_results$n_permutations, "permutations"))
    }

    # If gene_summary_data is provided, match p-values by gene name
    if (!is.null(gene_summary_data)) {
        gene_summary_genes <- gene_summary_data[[gene_col]]
        n_summary_genes <- length(gene_summary_genes)

        # Initialize p-value vectors
        matched_p_positive <- rep(NA_real_, n_summary_genes)
        matched_p_negative <- rep(NA_real_, n_summary_genes)

        # Create result mapping
        result_map <- setNames(1:length(cpp_results$gene_names), cpp_results$gene_names)

        # Match by gene name
        for (i in 1:n_summary_genes) {
            gene_name <- gene_summary_genes[i]
            if (!is.na(gene_name) && gene_name %in% names(result_map)) {
                result_idx <- result_map[gene_name]
                matched_p_positive[i] <- cpp_results$P_positive[result_idx]
                matched_p_negative[i] <- cpp_results$P_negative[result_idx]
            }
        }

        return(list(
            P_positive = matched_p_positive,
            P_negative = matched_p_negative,
            elapsed_time = elapsed_time,
            n_threads = 1, # Simplified version does not support multi-threading display
            engine = "C++",
            performance_stats = list(
                n_valid_genes = cpp_results$n_valid_genes,
                n_permutations = cpp_results$n_permutations,
                speed = paste(round(n_permutations / elapsed_time), "perms/sec")
            )
        ))
    } else {
        return(list(
            P_positive = cpp_results$P_positive,
            P_negative = cpp_results$P_negative,
            gene_names = cpp_results$gene_names,
            elapsed_time = elapsed_time,
            n_threads = 1, # Simplified version does not support multi-threading display
            engine = "C++",
            performance_stats = list(
                n_valid_genes = cpp_results$n_valid_genes,
                n_permutations = cpp_results$n_permutations,
                speed = paste(round(n_permutations / elapsed_time), "perms/sec")
            )
        ))
    }
}

#' Automatically Initialize C++ Engine
initialize_cpp_engine <- function() {
    cat("Initializing C++ Engine...\n")

    # First check if compiled
    if (exists("perform_cpp_permutation_test")) {
        cat("C++ Engine available, skipping compilation\n")
        return(TRUE)
    }

    # Step 1: Check and install packages
    env_success <- tryCatch(
        {
            setup_cpp_environment()
            TRUE
        },
        error = function(e) {
            cat("Environment setup failed:", e$message, "\n")
            FALSE
        }
    )

    if (!env_success) {
        return(FALSE)
    }

    # Step 2: Compile C++ code
    compile_success <- compile_cpp_engine()

    return(compile_success)
}
