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

calculate_gene_score_pairs <- function(diff_scores, gene_labels, unique_genes) {
    groups <- split(diff_scores, factor(gene_labels, levels = unique_genes))
    vapply(groups, function(scores) {
        scores <- sort(scores)
        k <- ceiling(2 * length(scores) / 3)
        c(positive = mean(rev(scores)[seq_len(k)]), negative = mean(scores[seq_len(k)]))
    }, numeric(2))
}

calculate_permutation_extremes <- function(seed, diff_scores, gene_labels, unique_genes, observed) {
    set.seed(seed)
    shuffled <- gene_labels[sample.int(length(gene_labels))]
    scores <- calculate_gene_score_pairs(diff_scores, shuffled, unique_genes)
    rbind(scores[1, ] >= observed[1, ], scores[2, ] <= observed[2, ])
}

perform_permutation <- function(
    sgrna_info_dt, gene_summary_dt, gene_col, N_perm = 1000L, n_cores = NULL,
    use_data_table = TRUE, progress_callback = NULL, min_sgrna_threshold = 3L,
    seed = NULL
) {
    N_perm <- validate_permutation_integer(N_perm, "N_perm")
    if (!nrow(gene_summary_dt)) {
        return(list(P_positive = numeric(), P_negative = numeric()))
    }
    prepared <- prepare_cpp_permutation_data(sgrna_info_dt, gene_col, "diff_score", min_sgrna_threshold)
    n_cores <- resolve_permutation_threads(n_cores, N_perm)
    observed <- calculate_gene_score_pairs(prepared$diff_scores, prepared$gene_labels, prepared$unique_genes)
    if (!is.null(seed)) set.seed(normalize_cpp_seed(seed))
    seeds <- sample.int(.Machine$integer.max, N_perm, replace = TRUE)
    counts <- matrix(0, nrow = 2L, ncol = length(prepared$unique_genes))
    cluster <- NULL
    if (n_cores > 1L) {
        cluster <- parallel::makeCluster(n_cores)
        on.exit(parallel::stopCluster(cluster), add = TRUE)
        parallel::clusterExport(
            cluster, c("calculate_gene_score_pairs", "calculate_permutation_extremes"),
            envir = environment(perform_permutation)
        )
    }
    completed <- 0L
    while (completed < N_perm) {
        end <- min(as.double(completed) + 100, N_perm)
        batch_seeds <- seeds[seq.int(completed + 1L, end)]
        arguments <- list(
            X = batch_seeds, fun = calculate_permutation_extremes,
            diff_scores = prepared$diff_scores, gene_labels = prepared$gene_labels,
            unique_genes = prepared$unique_genes, observed = observed
        )
        if (is.null(cluster)) {
            arguments$FUN <- arguments$fun
            arguments$fun <- NULL
            results <- do.call(lapply, arguments)
        } else {
            results <- do.call(parallel::parLapply, c(list(cl = cluster), arguments))
        }
        for (result in results) counts <- counts + result
        completed <- end
        if (!is.null(progress_callback)) {
            progress_callback(sprintf("Permutation Test: %d/%d (%.1f%%)", completed, N_perm, 100 * completed / N_perm))
        }
    }
    indices <- match(as.character(gene_summary_dt[[gene_col]]), prepared$unique_genes)
    list(
        P_positive = ((counts[1, ] + 1) / (as.double(N_perm) + 1))[indices],
        P_negative = ((counts[2, ] + 1) / (as.double(N_perm) + 1))[indices]
    )
}
