# mageck_wrappers.R
# Wrappers for running MAGeCK RRA and MLE analyses

# Helper to save minimal count table for MAGeCK
# Helper to locate MAGeCK executable
# Helper to locate MAGeCK executable
# Helper to locate MAGeCK executable with Caching
find_mageck_executable <- function() {
    config_file <- "config/mageck_path.txt"

    # 1. Check Config File First
    if (file.exists(config_file)) {
        stored_path <- trimws(readLines(config_file, warn = FALSE)[1])
        if (nzchar(stored_path) && file.exists(stored_path)) {
            message("[MAGeCK Config] Using cached path from config: ", stored_path)
            return(stored_path)
        }
    }

    # 2. Check System/Conda (Standard Search)
    found_path <- NULL

    # Try system path
    sys_path <- Sys.which("mageck")
    if (nzchar(sys_path)) found_path <- sys_path

    # Try Conda if not in system path
    if (is.null(found_path)) {
        home_dir <- Sys.getenv("HOME")

        # Check for 'conda' binary
        conda_locations <- c(
            # User home paths
            file.path(home_dir, "anaconda3", "bin", "conda"),
            file.path(home_dir, "miniconda3", "bin", "conda"),
            file.path(home_dir, "opt", "anaconda3", "bin", "conda"),
            file.path(home_dir, "opt", "miniconda3", "bin", "conda"),
            file.path(home_dir, "mambaforge", "bin", "conda"),
            file.path(home_dir, "miniforge3", "bin", "conda"),

            # System paths (Root level)
            "/opt/anaconda3/bin/conda",
            "/opt/miniconda3/bin/conda",
            "/usr/local/anaconda3/bin/conda",
            "/usr/local/miniconda3/bin/conda",
            "/usr/local/bin/conda",
            "/opt/homebrew/bin/conda"
        )

        found_conda <- NULL
        for (c_path in conda_locations) {
            if (file.exists(c_path)) {
                found_conda <- c_path
                message("[MAGeCK Config] Found conda binary at: ", found_conda)
                break
            }
        }

        if (!is.null(found_conda)) {
            message("[MAGeCK Config] Attempting to locate mageck via 'conda run -n mageck-vispr'...")
            tryCatch(
                {
                    res <- system2(found_conda, args = c("run", "-n", "mageck-vispr", "which", "mageck"), stdout = TRUE, stderr = FALSE)
                    if (length(res) > 0) {
                        valid_paths <- res[grepl("^/", res)]
                        if (length(valid_paths) > 0 && file.exists(valid_paths[1])) {
                            found_path <- valid_paths[1]
                            message("[MAGeCK Config] Conda reported mageck at: ", found_path)
                        }
                    }
                },
                error = function(e) {
                    message("[MAGeCK Config] 'conda run' detection failed: ", e$message)
                }
            )
        }
    }

    # Try Static Paths if still not found
    if (is.null(found_path)) {
        home_dir <- Sys.getenv("HOME")
        common_roots <- c(
            file.path(home_dir, "anaconda3"),
            file.path(home_dir, "miniconda3"),
            file.path(home_dir, ".conda"),
            # System paths
            "/opt/anaconda3",
            "/opt/miniconda3",
            "/usr/local/anaconda3",
            "/usr/local/miniconda3"
        )
        for (root in common_roots) {
            guess <- file.path(root, "envs", "mageck-vispr", "bin", "mageck")
            if (file.exists(guess)) {
                found_path <- guess
                message("[MAGeCK Config] Found at static path: ", found_path)
                break
            }
        }
    }

    # 3. Save to Config and Return using found path or default to "mageck"
    final_path <- if (!is.null(found_path)) found_path else "mageck"

    if (!is.null(found_path)) {
        # Save validated path to config
        message("[MAGeCK Config] Caching valid path to ", config_file)
        if (!dir.exists("config")) dir.create("config")
        writeLines(final_path, config_file)
    } else {
        message("[MAGeCK Config] WARNING: Could not auto-locate 'mageck'. Using default command 'mageck'.")
    }

    return(final_path)
}

# Helper to execute MAGeCK command with real-time progress updates
execute_mageck_with_progress <- function(cmd_binary, args, shiny_session = NULL, progress_message = "Running MAGeCK...", temp_dir = NULL) {
    # Prepare logs
    if (is.null(temp_dir)) {
        idx <- which(args == "-n")
        if (length(idx) > 0) temp_dir <- dirname(args[idx + 1])
        if (length(temp_dir) == 0 || !dir.exists(temp_dir)) temp_dir <- tempdir()
    }

    stdout_log <- file.path(temp_dir, "stdout.log")
    stderr_log <- file.path(temp_dir, "stderr.log")

    # Clean old logs
    if (file.exists(stdout_log)) file.remove(stdout_log)
    if (file.exists(stderr_log)) file.remove(stderr_log)

    # Configure PATH for processx
    mageck_dir <- dirname(cmd_binary)
    current_path <- Sys.getenv("PATH")
    new_path <- paste(mageck_dir, current_path, sep = ":")
    env_vars <- c("PATH" = new_path)

    # Start process using processx
    bg_process <- processx::process$new(
        command = cmd_binary,
        args = args,
        stdout = stdout_log,
        stderr = stderr_log,
        env = env_vars
    )

    # Monitor loop with progress updates
    start_time <- Sys.time()
    poll_count <- 0

    while (bg_process$is_alive()) {
        Sys.sleep(0.5) # Poll every 0.5 seconds
        poll_count <- poll_count + 1

        if (!is.null(shiny_session)) {
            # Read last line from stderr for status
            detail_msg <- ""
            if (file.exists(stderr_log)) {
                tryCatch(
                    {
                        lines <- readLines(stderr_log, warn = FALSE)
                        if (length(lines) > 0) {
                            # Get last INFO line
                            info_lines <- grep("^INFO", lines, value = TRUE)
                            if (length(info_lines) > 0) {
                                last_info <- tail(info_lines, 1)
                                # Extract message after timestamp
                                detail_msg <- sub("^INFO\\s+@\\s+[^:]+:\\s*", "", last_info)
                                if (nchar(detail_msg) > 80) detail_msg <- substr(detail_msg, 1, 80)
                            }
                        }
                    },
                    error = function(e) {}
                )
            }

            if (detail_msg == "") {
                elapsed <- round(difftime(Sys.time(), start_time, units = "secs"))
                detail_msg <- paste0("Processing... (", elapsed, "s)")
            }

            # Simulate progress (pulse between 0.3 and 0.7)
            progress_val <- 0.3 + 0.4 * (sin(poll_count * 0.5) + 1) / 2

            shiny::setProgress(
                value = progress_val,
                message = progress_message,
                detail = detail_msg,
                session = shiny_session
            )
        }
    }

    # Get exit code
    exit_code <- bg_process$get_exit_status()

    # Final progress update
    if (!is.null(shiny_session)) {
        if (exit_code == 0) {
            shiny::setProgress(value = 0.95, message = progress_message, detail = "Completed successfully", session = shiny_session)
        } else {
            shiny::setProgress(value = 1.0, message = "MAGeCK Failed", detail = "Check error logs", session = shiny_session)
        }
    }

    return(exit_code)
}

# Helper to save minimal count table for MAGeCK
prepare_mageck_input_file <- function(data_df, id_cols, count_cols, filepath) {
    # Keep only ID and Count columns
    # MAGeCK expects Tab-separated: sgRNA  Gene  Sample1  Sample2 ...

    # Select columns
    df_to_save <- data_df[, c(id_cols, count_cols)]

    # Important: MAGeCK expects raw counts.
    # Replace NAs with 0 (as missing data implies 0 counts in raw files often),
    # but DO NOT apply +1 or normalization here.

    # Convert count columns to numeric (just in case) and handle NAs
    for (col in count_cols) {
        if (col %in% names(df_to_save)) {
            df_to_save[[col]] <- as.numeric(df_to_save[[col]])
            df_to_save[[col]][is.na(df_to_save[[col]])] <- 0
        }
    }

    # Write to tab-separated file
    write.table(df_to_save, file = filepath, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Wrapper for MAGeCK RRA
perform_mageck_rra_analysis_wrapper <- function(raw_data_df,
                                                gRNA_col,
                                                gene_col,
                                                condition1_replicate_cols, # Treatment (Test)
                                                condition2_replicate_cols, # Control
                                                output_prefix_base = "mageck_rra_run",
                                                shiny_session = NULL) {
    # 1. Prepare Workspace
    temp_dir <- file.path(tempdir(), paste0("mageck_rra_", format(Sys.time(), "%H%M%S")))
    if (!dir.exists(temp_dir)) dir.create(temp_dir)

    input_file_path <- file.path(temp_dir, "counts.txt")
    output_prefix <- file.path(temp_dir, "output")

    # 2. Write Input File
    prepare_mageck_input_file(raw_data_df, c(gRNA_col, gene_col), c(condition1_replicate_cols, condition2_replicate_cols), input_file_path)

    # 3. Construct Command
    # mageck test -k counts.txt -t Treat1,Treat2 -c Ctrl1,Ctrl2 -n output
    treat_str <- paste(condition1_replicate_cols, collapse = ",")
    ctrl_str <- paste(condition2_replicate_cols, collapse = ",")

    # Check if mageck is accessible (auto-detect)
    mageck_cmd <- find_mageck_executable()

    args <- c(
        "test",
        "-k", input_file_path,
        "-t", treat_str,
        "-c", ctrl_str,
        "-n", output_prefix,
        "--remove-zero", "both", # Recommended
        "--remove-zero-threshold", "0"
    )

    if (!is.null(shiny_session)) {
        msg <- paste0("Executing MAGeCK RRA: ", mageck_cmd)
        shiny::setProgress(value = 0.2, message = "Running MAGeCK RRA...", detail = msg)
    }

    # 4. Run Command with Progress
    exit_code <- execute_mageck_with_progress(
        mageck_cmd,
        args,
        shiny_session,
        progress_message = "Running MAGeCK RRA...",
        temp_dir = temp_dir
    )

    if (exit_code != 0) {
        err_msg <- paste("MAGeCK RRA execution failed with code", exit_code)
        if (file.exists(file.path(temp_dir, "stderr.log"))) {
            err_log <- readLines(file.path(temp_dir, "stderr.log"))
            err_msg <- paste(err_msg, "\nLast stderr:", paste(tail(err_log, 5), collapse = "\n"))
            err_msg <- paste(err_msg, "\nAttempted binary:", mageck_cmd)
        }
        stop(err_msg)
    }

    if (!is.null(shiny_session)) shiny::setProgress(value = 0.7, message = "Parsing MAGeCK Results...", detail = "Importing output...")

    # 5. Parse Output
    gene_summary_file <- paste0(output_prefix, ".gene_summary.txt")
    sgrna_summary_file <- paste0(output_prefix, ".sgrna_summary.txt")

    if (!file.exists(gene_summary_file)) stop("MAGeCK RRA output file not found: gene_summary.txt")

    # Read Gene Summary
    gene_res <- read.delim(gene_summary_file, stringsAsFactors = FALSE)

    # Map to ATOP columns
    # MAGeCK RRA: id, num, neg|score, neg|p-value, neg|fdr, neg|rank, neg|goodsgrna, neg|lfc, pos|score, ...

    gene_summary_df <- data.frame(
        Gene = gene_res$id,
        n_sgrnas = gene_res$num,
        GeneNegativeScore = gene_res$neg.lfc,
        P_negative = gene_res$neg.p.value,
        P_negative_adj_bh = gene_res$neg.fdr, # MAGeCK calculates FDR
        GenePositiveScore = gene_res$pos.lfc, # Usually same as neg|lfc
        P_positive = gene_res$pos.p.value,
        P_positive_adj_bh = gene_res$pos.fdr,
        GenePositiveRank = gene_res$pos.rank,
        GeneNegativeRank = gene_res$neg.rank
    )
    # Ensure column name consistency with parsing code if 'Gene' col was custom
    names(gene_summary_df)[1] <- gene_col

    # Read sgRNA Summary (for completeness)
    sgrna_res <- read.delim(sgrna_summary_file, stringsAsFactors = FALSE)
    # Need to map sgrna_res to processed_sg_data structure roughly
    # ATOP expects: sgrna_id, gene_id, diff_score, sgPositiveRank, sgNegativeRank
    # MAGeCK sgrna_summary: sgrna, Gene, control_count, treatment_count, control_mean, treat_mean, LFC, ...
    # We can use LFC as diff_score.

    # Note: sgrna_res column names depend entirely on MAGeCK version? Usually standard.
    # Columns: sgrna, Gene, control_count, treatment_count, control_mean, treat_mean, LFC, control_var, adjust_p_value, ...

    processed_sg_data <- data.frame(
        gRNA_col = sgrna_res$sgrna,
        gene_col = sgrna_res$Gene,
        diff_score = sgrna_res$LFC
    )
    names(processed_sg_data) <- c(gRNA_col, gene_col, "diff_score")

    # Add rankings to processed_sg_data (per gene)
    processed_sg_data <- processed_sg_data %>%
        dplyr::group_by(!!dplyr::sym(gene_col)) %>%
        dplyr::mutate(
            sgPositiveRank = rank(-diff_score, ties.method = "first"),
            sgNegativeRank = rank(diff_score, ties.method = "first")
        ) %>%
        dplyr::ungroup()

    # Standard Normalization (Median of Ratios) in R for plotting consistency
    # (Since MAGeCK test doesn't output normalized table file easily without 'mageck count')
    # Use existing ATOP normalization logic just for the `normalized_counts` return object
    # We can reuse the ATOP function but we need to source it or copy small logic.
    # For now, let's copy the specific normalization block from ATOP logic to generic helper.

    normalized_counts <- normalize_median_of_ratios(raw_data_df, c(gRNA_col, gene_col), c(condition1_replicate_cols, condition2_replicate_cols))

    return(list(
        processed_sg_data = processed_sg_data,
        gene_summary_data = gene_summary_df,
        lfc_details = raw_data_df,
        normalized_counts = normalized_counts,
        sgrna_count_validation_report = NULL,
        # Raw MAGeCK outputs for display
        raw_mageck_gene_summary = gene_res,
        raw_mageck_sgrna_summary = sgrna_res,
        params = list(
            gRNA_col = gRNA_col,
            gene_col = gene_col,
            analysis_tool = "MAGeCK RRA"
        )
    ))
}

# Wrapper for MAGeCK MLE
perform_mageck_mle_analysis_wrapper <- function(raw_data_df,
                                                gRNA_col,
                                                gene_col,
                                                condition1_replicate_cols, # Treatment (Test)
                                                condition2_replicate_cols, # Control
                                                output_prefix_base = "mageck_mle_run",
                                                shiny_session = NULL) {
    # 1. Prepare Workspace
    temp_dir <- file.path(tempdir(), paste0("mageck_mle_", format(Sys.time(), "%H%M%S")))
    if (!dir.exists(temp_dir)) dir.create(temp_dir)

    input_file_path <- file.path(temp_dir, "counts.txt")
    design_file_path <- file.path(temp_dir, "design.txt")
    output_prefix <- file.path(temp_dir, "output")

    samples <- c(condition1_replicate_cols, condition2_replicate_cols)

    # 2. Write Input File
    prepare_mageck_input_file(raw_data_df, c(gRNA_col, gene_col), samples, input_file_path)

    # 3. Create Design Matrix
    # Format: Sample, baseline, beta_variable (e.g. condition)
    # Baseline: 1 (always)
    # Condition: 0 for Control, 1 for Treatment

    design_df <- data.frame(
        Samples = samples,
        baseline = 1,
        Condition = ifelse(samples %in% condition1_replicate_cols, 1, 0)
    )
    write.table(design_df, file = design_file_path, sep = "\t", quote = FALSE, row.names = FALSE)

    # 4. Construct Command
    # mageck mle -k counts.txt -d design.txt -n output --design-id samples baseline Condition

    mageck_cmd <- find_mageck_executable()

    # Detect optimal thread count (same logic as ATOP C++)
    n_cores <- parallel::detectCores()
    optimal_threads <- max(1, n_cores - 1) # Leave one core for system

    args <- c(
        "mle",
        "-k", input_file_path,
        "-d", design_file_path,
        "-n", output_prefix,
        "--threads", as.character(optimal_threads)
    )

    if (!is.null(shiny_session)) {
        msg <- paste0("Executing MAGeCK MLE with ", optimal_threads, " threads: ", mageck_cmd)
        shiny::setProgress(value = 0.2, message = "Running MAGeCK MLE...", detail = msg)
    }

    # 5. Run Command with Progress
    exit_code <- execute_mageck_with_progress(
        mageck_cmd,
        args,
        shiny_session,
        progress_message = "Running MAGeCK MLE...",
        temp_dir = temp_dir
    )

    if (exit_code != 0) {
        err_msg <- paste("MAGeCK MLE execution failed with code", exit_code)
        if (file.exists(file.path(temp_dir, "stderr.log"))) {
            err_log <- readLines(file.path(temp_dir, "stderr.log"))
            err_msg <- paste(err_msg, "\nLast stderr:", paste(tail(err_log, 5), collapse = "\n"))
            err_msg <- paste(err_msg, "\nAttempted binary:", mageck_cmd)
        }
        stop(err_msg)
    }

    if (!is.null(shiny_session)) shiny::setProgress(value = 0.7, message = "Parsing MAGeCK Results...", detail = "Importing output...")

    # 6. Parse Output
    gene_summary_file <- paste0(output_prefix, ".gene_summary.txt")
    if (!file.exists(gene_summary_file)) stop("MAGeCK MLE output file not found: gene_summary.txt")

    gene_res <- read.delim(gene_summary_file, stringsAsFactors = FALSE)

    # Columns: Gene, sgRNA, Condition|beta, Condition|z, Condition|p-value, ...
    # Column names in R might become Condition.beta, Condition.z etc.

    # Identify the beta/p info columns
    # We named variable "Condition" in design matrix
    beta_col <- "Condition.beta"
    pval_col <- "Condition.wald.p.value" # Use Wald p-value usually for MLE
    fdr_col <- "Condition.wald.fdr"

    if (!beta_col %in% names(gene_res)) {
        # Fallback search
        beta_col <- grep("\\.beta$", names(gene_res), value = TRUE)[1]
        pval_col <- grep("wald\\.p\\.value$", names(gene_res), value = TRUE)[1]
        fdr_col <- grep("wald\\.fdr$", names(gene_res), value = TRUE)[1]
    }

    gene_summary_df <- data.frame(
        Gene = gene_res$Gene,
        n_sgrnas = gene_res$sgRNA,
        GenePositiveScore = gene_res[[beta_col]], # Beta as score
        GeneNegativeScore = gene_res[[beta_col]], # Beta as score (Plotting logic handles sign)
        P_positive = gene_res[[pval_col]],
        P_negative = gene_res[[pval_col]],
        P_positive_adj_bh = gene_res[[fdr_col]],
        P_negative_adj_bh = gene_res[[fdr_col]],

        # Ranks inferred from Score (Beta)
        GenePositiveRank = rank(-gene_res[[beta_col]]),
        GeneNegativeRank = rank(gene_res[[beta_col]])
    )
    names(gene_summary_df)[1] <- gene_col

    # Processed sgRNA data? MLE doesn't output per-sgRNA LFC in summary easily
    # It assumes sgRNA counts.
    # We can just return the input parsed sgrna?
    # Or create dummy diff_score for compatibility

    processed_sg_data <- data.frame(
        gRNA_col = paste0("dummy_", 1:10), # Impossible to reconstruct easily without another file
        gene_col = "dummy",
        diff_score = 0
    )
    # Actually, visualizers need sgRNA data.
    # We can compute LFCs using ATOP normalizer to populate this field, even if testing is MLE.
    # This ensures 'sgRNA Paired Plot' in App works (which uses simple LFC).

    # Compute Normalized Counts & Simple LFC for plotting compatibility
    normalized_counts <- normalize_median_of_ratios(raw_data_df, c(gRNA_col, gene_col), samples)

    # Simple LFC for processed_sg_data (needed for volcano tooltips or sgrna plots)
    # Reuse simple LFC logic
    lfc_matrix <- calculate_simple_lfc(normalized_counts, condition1_replicate_cols, condition2_replicate_cols, samples)

    processed_sg_data <- raw_data_df[, c(gRNA_col, gene_col)]
    processed_sg_data$diff_score <- rowMeans(lfc_matrix, na.rm = TRUE)

    # Add rankings
    processed_sg_data <- processed_sg_data %>%
        dplyr::group_by(!!dplyr::sym(gene_col)) %>%
        dplyr::mutate(
            sgPositiveRank = rank(-diff_score, ties.method = "first"),
            sgNegativeRank = rank(diff_score, ties.method = "first")
        ) %>%
        dplyr::ungroup()

    return(list(
        processed_sg_data = processed_sg_data,
        gene_summary_data = gene_summary_df,
        lfc_details = raw_data_df,
        normalized_counts = normalized_counts,
        sgrna_count_validation_report = NULL,
        # Raw MAGeCK outputs for display
        raw_mageck_gene_summary = gene_res,
        raw_mageck_sgrna_summary = processed_sg_data, # MLE doesn't have separate sgrna summary, use computed
        params = list(
            gRNA_col = gRNA_col,
            gene_col = gene_col,
            analysis_tool = "MAGeCK MLE"
        )
    ))
}

# --- Shared Utils ---

normalize_median_of_ratios <- function(data_df, id_cols, count_cols) {
    # Simple R implementation of Median of Ratios
    # Assumes data_df has id_cols and count_cols

    counts_mat <- as.matrix(data_df[, count_cols])
    counts_mat[is.na(counts_mat)] <- 0
    counts_mat <- counts_mat + 1 # Pseudo count

    # GeoMean
    log_counts <- log(counts_mat)
    log_geom_means <- rowMeans(log_counts)
    geom_means <- exp(log_geom_means)

    # Ratios
    # Handle geom_mean 0? (Shouldn't happen with +1, except if empty?)
    ratios <- sweep(counts_mat, 1, geom_means, "/")

    # Size Factors
    size_factors <- apply(ratios, 2, median)

    # Normalize
    norm_counts <- sweep(counts_mat, 2, size_factors, "/")

    # Reassemble
    res_df <- data_df[, id_cols, drop = FALSE]
    res_df <- cbind(res_df, as.data.frame(norm_counts))
    return(res_df)
}

calculate_simple_lfc <- function(norm_counts_df, group1, group2, all_samples) {
    # Extract
    g1_data <- norm_counts_df[, group1, drop = FALSE]
    g2_data <- norm_counts_df[, group2, drop = FALSE]

    # +1 already done in normalization input if using regular counts,
    # but norm_counts are output. Safe to add small epsilon or assume no 0s if derived from +1 input.
    # Usually +1 is done before log.

    # Pairwise LFC or Mean LFC?
    # ATOP uses pairwise if paired? No, just all vs all or mean.
    # Let's do simple mean(log2(T)) - mean(log2(C)) for robustness

    # Return single column approx or replicate matrix?
    # Returning matrix for compatibility if wrapper expectation is matrix.
    # But usually just need rowMeans for diff_score.

    # Let's emulate ATOP's n_reps structure if possible, usually n_reps = len(group1)
    # Delta LFC i k = log2( (T_k + 1) / (C_k + 1) ) if paired.
    # If not paired, this is arbitrary.
    # For compatibility, return matrix of length n_reps (min length)

    n_reps <- min(length(group1), length(group2))
    lfc_list <- list()
    for (k in 1:n_reps) {
        lfc_list[[k]] <- log2((g1_data[[k]] + 1) / (g2_data[[k]] + 1))
    }
    return(do.call(cbind, lfc_list))
}
