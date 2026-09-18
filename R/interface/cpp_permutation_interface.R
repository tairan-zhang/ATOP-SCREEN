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

setup_cpp_environment <- function() {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Install Rcpp using scripts/install_dependencies.R.")
    }
    TRUE
}

compile_cpp_engine <- function() {
    tryCatch({
        setup_cpp_environment()
        Rcpp::sourceCpp(
            "src/cpp_permutation_engine.cpp",
            env = environment(compile_cpp_engine),
            verbose = FALSE
        )
        TRUE
    }, error = function(e) {
        message("C++ compilation failed: ", conditionMessage(e))
        FALSE
    })
}

initialize_cpp_engine <- function() {
    if (exists("perform_cpp_permutation_test", mode = "function")) return(TRUE)
    compile_cpp_engine()
}

validate_permutation_integer <- function(value, name, minimum = 1L) {
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
        value != floor(value) || value < minimum || value > .Machine$integer.max) {
        stop(name, " must be an integer between ", minimum, " and ", .Machine$integer.max)
    }
    as.integer(value)
}

resolve_permutation_threads <- function(n_threads = NULL, n_permutations = 1L) {
    available <- parallel::detectCores()
    if (length(available) != 1L || !is.finite(available)) available <- 1L
    if (is.null(n_threads) || identical(n_threads, 0) || identical(n_threads, 0L)) {
        n_threads <- max(1L, min(available - 1L, 8L))
    }
    n_threads <- validate_permutation_integer(n_threads, "n_threads")
    as.integer(min(n_threads, available, n_permutations))
}

normalize_cpp_seed <- function(seed = NULL) {
    if (is.null(seed)) return(sample.int(.Machine$integer.max, 1L))
    validate_permutation_integer(seed, "seed", 0L)
}

prepare_cpp_permutation_data <- function(
    sgrna_data, gene_col = "Gene", score_col = "diff_score", min_sgrna_threshold = 3L
) {
    min_sgrna_threshold <- validate_permutation_integer(min_sgrna_threshold, "min_sgrna_threshold")
    data <- as.data.frame(sgrna_data)
    if (!all(c(gene_col, score_col) %in% names(data))) stop("Gene or score column is missing.")
    if (!is.numeric(data[[score_col]])) stop("Score column must be numeric.")
    gene_labels <- as.character(data[[gene_col]])
    valid_rows <- !is.na(gene_labels) & nzchar(trimws(gene_labels)) & is.finite(data[[score_col]])
    data <- data[valid_rows, , drop = FALSE]
    gene_labels <- gene_labels[valid_rows]
    gene_counts <- table(gene_labels)
    valid_genes <- names(gene_counts)[gene_counts >= min_sgrna_threshold]
    if (!length(valid_genes)) stop("No genes have enough valid sgRNAs.")
    keep <- gene_labels %in% valid_genes
    list(
        diff_scores = data[[score_col]][keep],
        gene_labels = gene_labels[keep],
        unique_genes = sort(valid_genes)
    )
}

format_cpp_permutation_result <- function(result, gene_summary_data, gene_col, elapsed_time) {
    indices <- seq_along(result$gene_names)
    if (!is.null(gene_summary_data)) {
        if (!gene_col %in% names(gene_summary_data)) stop("Gene column is missing from gene summary.")
        indices <- match(as.character(gene_summary_data[[gene_col]]), result$gene_names)
    }
    list(
        P_positive = result$P_positive[indices],
        P_negative = result$P_negative[indices],
        gene_names = result$gene_names[indices],
        positive_extreme_counts = result$positive_extreme_counts[indices],
        negative_extreme_counts = result$negative_extreme_counts[indices],
        valid_permutation_counts = result$valid_permutation_counts[indices],
        elapsed_time = elapsed_time,
        n_threads = result$n_threads,
        seed = result$seed,
        engine = "C++17",
        performance_stats = list(
            n_valid_genes = result$n_valid_genes,
            n_permutations = result$n_permutations,
            speed = paste(round(result$n_permutations / max(elapsed_time, .Machine$double.eps)), "perms/sec")
        )
    )
}

perform_cpp_permutation <- function(
    sgrna_data, gene_summary_data = NULL, gene_col = "Gene", score_col = "diff_score",
    n_permutations = 1000L, min_sgrna_threshold = 3L, n_threads = NULL,
    progress_callback = NULL, seed = NULL, batch_size = NULL
) {
    n_permutations <- validate_permutation_integer(n_permutations, "n_permutations")
    prepared <- prepare_cpp_permutation_data(sgrna_data, gene_col, score_col, min_sgrna_threshold)
    n_threads <- resolve_permutation_threads(n_threads, n_permutations)
    seed <- normalize_cpp_seed(seed)
    if (is.null(batch_size)) {
        batch_size <- if (is.null(progress_callback)) n_permutations else max(100L, n_permutations %/% 20L)
    }
    batch_size <- validate_permutation_integer(batch_size, "batch_size")
    if (!initialize_cpp_engine()) stop("C++ engine is unavailable.")

    started <- proc.time()[["elapsed"]]
    positive_counts <- negative_counts <- numeric(length(prepared$unique_genes))
    completed <- 0L
    while (completed < n_permutations) {
        current_size <- min(batch_size, n_permutations - completed)
        result <- perform_cpp_permutation_test(
            diff_scores = prepared$diff_scores,
            gene_labels = prepared$gene_labels,
            unique_genes = prepared$unique_genes,
            n_permutations = current_size,
            min_sgrna_threshold = min_sgrna_threshold,
            seed = seed,
            show_progress = FALSE,
            n_threads = n_threads,
            permutation_offset = completed
        )
        positive_counts <- positive_counts + result$positive_extreme_counts
        negative_counts <- negative_counts + result$negative_extreme_counts
        completed <- completed + current_size
        if (!is.null(progress_callback)) {
            progress_callback(sprintf("Permutation Test: %d/%d (%.1f%%)", completed, n_permutations, 100 * completed / n_permutations))
        }
    }

    # Apply the pseudocount once, after all batches have contributed raw counts.
    result$P_positive <- (positive_counts + 1) / (as.double(n_permutations) + 1)
    result$P_negative <- (negative_counts + 1) / (as.double(n_permutations) + 1)
    result$positive_extreme_counts <- positive_counts
    result$negative_extreme_counts <- negative_counts
    result$valid_permutation_counts <- rep(n_permutations, length(positive_counts))
    result$n_permutations <- n_permutations
    result$n_threads <- min(n_threads, batch_size)
    format_cpp_permutation_result(result, gene_summary_data, gene_col, proc.time()[["elapsed"]] - started)
}
