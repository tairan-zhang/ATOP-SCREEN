suppressPackageStartupMessages(library(dplyr))
source("R/interface/cpp_permutation_interface.R")
source("R/interface/cpp_progress_wrapper.R")
source("R/analysis/permutation_functions.R")
source("R/analysis/engine_selector.R")
source("R/analysis/crispr_analysis_functions.R")

assert_error <- function(expression) {
    stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
}
stopifnot(initialize_cpp_engine())
data <- data.frame(Gene = rep(sprintf("G%02d", 1:10), each = 6), diff_score = seq_len(60))
summary <- data.frame(Gene = c("G10", "G01", "absent", "G05"))
run <- function(batch_size, threads = 1L) {
    perform_cpp_permutation(
        data, summary, n_permutations = 1003L, n_threads = threads,
        seed = 42L, batch_size = batch_size
    )
}
direct <- run(1003L)
for (result in list(run(100L), run(137L, 2L))) {
    stopifnot(identical(direct$P_positive, result$P_positive))
    stopifnot(identical(direct$P_negative, result$P_negative))
    stopifnot(identical(direct$positive_extreme_counts, result$positive_extreme_counts))
}
stopifnot(is.na(direct$P_positive[3]))
stopifnot(all(direct$P_positive == (direct$positive_extreme_counts + 1) / 1004, na.rm = TRUE))
stopifnot(all(direct$P_negative == (direct$negative_extreme_counts + 1) / 1004, na.rm = TRUE))
stopifnot(direct$positive_extreme_counts[1] == 0, direct$P_positive[1] == 1 / 1004)

messages <- character()
wrapped <- perform_cpp_permutation_with_progress(
    data, summary, n_permutations = 1003L, n_threads = 2L, seed = 42L,
    progress_callback = function(message) messages <<- c(messages, message)
)
stopifnot(identical(direct$P_positive, wrapped$P_positive), length(messages) > 1)
stopifnot(grepl("100.0%", tail(messages, 1), fixed = TRUE))

scores <- calculate_observed_scores_cpp(
    c(5, 1, 3, -1, 2, 7, 4), c("A", "A", "B", "B", "B", "B", "B"), c("B", "A"), 2L
)
stopifnot(isTRUE(all.equal(scores$observed_pos_scores, c(4, 3))))
stopifnot(isTRUE(all.equal(scores$observed_neg_scores, c(2, 3))))
assert_error(perform_cpp_permutation_test(c(1, Inf), c("A", "A"), "A"))
assert_error(perform_cpp_permutation_test(1, "A", c("A", "A")))
assert_error(perform_cpp_permutation_test(1, NA_character_, "A"))
assert_error(perform_cpp_permutation(data, n_permutations = 0))
assert_error(perform_cpp_permutation(data, n_permutations = 1.5))

tied <- transform(data, diff_score = 0)
ties <- perform_cpp_permutation(tied, n_permutations = 101L, seed = 1L, batch_size = 30L)
stopifnot(all(ties$P_positive == 1), all(ties$P_negative == 1))
dirty <- rbind(data, data.frame(Gene = c("bad", "bad", NA), diff_score = c(Inf, NA, 1)))
stopifnot(length(prepare_cpp_permutation_data(dirty)$diff_scores) == 60L)
empty <- perform_cpp_permutation(data, data.frame(Gene = character()), n_permutations = 1L, seed = 42L)
stopifnot(length(empty$P_positive) == 0L)

r_one <- perform_permutation(data, summary, "Gene", 211L, 1L, seed = 42L)
r_two <- perform_permutation(data, summary, "Gene", 211L, 2L, seed = 42L, use_data_table = FALSE)
stopifnot(identical(r_one, r_two))
r_tied <- perform_permutation(tied, data.frame(Gene = unique(data$Gene)), "Gene", 21L, 1L, seed = 1L)
stopifnot(all(r_tied$P_positive == 1), all(r_tied$P_negative == 1))
selected <- engine_selector(data, summary, "Gene", 1003L, 1L, seed = 42L)
stopifnot(identical(selected$P_positive, direct$P_positive))

raw <- data.frame(gRNA = paste0("s", 1:12), Gene = rep(c("A", "B", "C"), each = 4),
                  treatment = c(12, 14, 16, 18, 3, 4, 5, 6, 8, 9, 10, 11),
                  control = rep(10, 12))
pipeline <- perform_crispr_screen_analysis(
    raw, condition1_replicate_cols = "treatment", condition2_replicate_cols = "control", N_perm = 21L
)
stopifnot(nrow(pipeline$processed_sg_data) == 12L, nrow(pipeline$gene_summary_data) == 3L)
stopifnot(all(is.finite(pipeline$gene_summary_data$P_positive)))
stopifnot(isTRUE(all.equal(pipeline$gene_summary_data$P_positive_adj_bh,
                          p.adjust(pipeline$gene_summary_data$P_positive, "BH"))))
skipped <- perform_crispr_screen_analysis(
    raw, N_perm = 0L, skip_normalization = TRUE, diff_score_col1 = "treatment", diff_score_col2 = "control"
)
stopifnot(all(is.na(skipped$gene_summary_data$P_positive)))
assert_error(perform_crispr_screen_analysis(
    rbind(raw, raw[1, ]), N_perm = 0L, skip_normalization = TRUE,
    diff_score_col1 = "treatment", diff_score_col2 = "control"
))
cat("PASS: batching, threads, progress, zero extremes, tails, ties, validation, R fallback, and full ATOP pipeline.\n")

small <- data.frame(Gene = rep(c("A", "B"), each = 3), diff_score = 1:6)
draws <- combn(1:6, 3)
positive_null <- apply(draws, 2, function(x) mean(sort(x, decreasing = TRUE)[1:2]))
negative_null <- apply(draws, 2, function(x) mean(sort(x)[1:2]))
expected_positive <- vapply(c(2.5, 5.5), function(x) mean(positive_null >= x), numeric(1))
expected_negative <- vapply(c(1.5, 4.5), function(x) mean(negative_null <= x), numeric(1))
estimated <- perform_cpp_permutation(small, n_permutations = 20000L, seed = 42L, n_threads = 2L)
stopifnot(max(abs(estimated$P_positive - expected_positive)) < 0.02)
stopifnot(max(abs(estimated$P_negative - expected_negative)) < 0.02)

session <- shiny::MockShinySession$new()
shiny::withReactiveDomain(session, {
    shiny::withProgress(
        message = "Regression test", session = session,
        expr = {
            session_result <- perform_cpp_permutation_with_progress(
                data, summary, n_permutations = 1003L, seed = 42L,
                n_threads = 2L, shiny_session = session
            )
        }
    )
})
session$close()
stopifnot(identical(session_result$P_positive, direct$P_positive))
cat("PASS: exhaustive-null calibration and Shiny-session batch path.\n")
