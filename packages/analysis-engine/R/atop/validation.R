# Input contracts for the ATOP workflow. Keeping these checks outside the
# workflow makes new front ends and batch callers share identical semantics.

atop_is_column_name <- function(value) is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)

atop_scalar <- function(value, label, minimum = -Inf, integer = FALSE) {
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value) || value < minimum || (integer && value != floor(value))) {
        stop(label, " must be a ", if (integer) "whole " else "", "number greater than or equal to ", minimum, ".", call. = FALSE)
    }
    if (integer) as.integer(value) else value
}

atop_replicate_columns <- function(columns, label, available_columns) {
    columns <- as.character(columns %||% character())
    columns <- columns[!is.na(columns) & nzchar(columns)]
    if (length(columns) == 0L) stop(label, " replicate columns are required.", call. = FALSE)
    duplicated_columns <- unique(columns[duplicated(columns)])
    if (length(duplicated_columns) > 0L) stop(label, " repeats columns: ", paste(duplicated_columns, collapse = ", "), call. = FALSE)
    missing <- setdiff(columns, available_columns)
    if (length(missing) > 0L) stop(label, " contains unknown columns: ", paste(missing, collapse = ", "), call. = FALSE)
    columns
}

atop_validate_request <- function(raw_data_df, gRNA_col, gene_col, sequence_col,
                                  condition1_replicate_cols, condition2_replicate_cols,
                                  N_perm, pseudo_count_lfc, min_sgrna_threshold,
                                  user_engine_choice, skip_normalization,
                                  diff_score_col1, diff_score_col2, initial_progress_value_abs, comparison_mode = "auto") {
    if (!is.data.frame(raw_data_df)) stop("Input data must be a data frame.", call. = FALSE)
    for (column in c(gRNA_col, gene_col)) {
        if (!atop_is_column_name(column) || !column %in% names(raw_data_df)) stop("A valid sgRNA and gene column must be selected.", call. = FALSE)
    }
    if (!atop_is_column_name(sequence_col)) sequence_col <- NULL
    if (!is.null(sequence_col) && !sequence_col %in% names(raw_data_df)) stop("The selected sequence column is not present in the input data.", call. = FALSE)
    skip_normalization <- isTRUE(skip_normalization)
    mode <- match.arg(comparison_mode, c("auto", "paired", "shared_baseline", "group_means"))
    if (skip_normalization) {
        mode <- "precomputed"
        if (!atop_is_column_name(diff_score_col1) || !atop_is_column_name(diff_score_col2)) stop("Skip-normalization mode requires two differential-score columns.", call. = FALSE)
        unknown <- setdiff(c(diff_score_col1, diff_score_col2), names(raw_data_df))
        if (length(unknown) > 0L) stop("Differential-score columns are missing: ", paste(unknown, collapse = ", "), call. = FALSE)
        condition1_replicate_cols <- character()
        condition2_replicate_cols <- character()
    } else {
        condition1_replicate_cols <- atop_replicate_columns(condition1_replicate_cols, "Condition", names(raw_data_df))
        condition2_replicate_cols <- atop_replicate_columns(condition2_replicate_cols, "Control", names(raw_data_df))
        endpoints <- length(condition1_replicate_cols)
        baselines <- length(condition2_replicate_cols)
        if (mode == "auto") mode <- if (endpoints == baselines) "paired" else if (baselines == 1L) "shared_baseline" else "group_means"
        if (mode == "paired" && endpoints != baselines) stop("Paired LFC requires equally sized groups.")
        if (mode == "shared_baseline" && baselines != 1L) stop("Shared baseline requires one baseline sample.")
        overlap <- intersect(condition1_replicate_cols, condition2_replicate_cols)
        if (length(overlap) > 0L) stop("Condition and control cannot reuse count columns: ", paste(overlap, collapse = ", "), call. = FALSE)
    }
    list(
        comparison_mode = mode,
        gRNA_col = gRNA_col, gene_col = gene_col, sequence_col = sequence_col,
        condition1_replicate_cols = condition1_replicate_cols, condition2_replicate_cols = condition2_replicate_cols,
        N_perm = atop_scalar(N_perm, "Number of permutations", 0, TRUE),
        pseudo_count_lfc = atop_scalar(pseudo_count_lfc, "LFC pseudo-count", 0),
        min_sgrna_threshold = atop_scalar(min_sgrna_threshold, "Minimum sgRNA count per gene", 1, TRUE),
        user_engine_choice = if (identical(user_engine_choice, "r_parallel")) "r_parallel" else "cpp",
        skip_normalization = skip_normalization, diff_score_col1 = diff_score_col1, diff_score_col2 = diff_score_col2,
        initial_progress_value_abs = atop_scalar(initial_progress_value_abs, "Initial progress value", 0)
    )
}
