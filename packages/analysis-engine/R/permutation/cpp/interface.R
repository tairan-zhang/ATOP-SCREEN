setup_cpp_environment <- function() {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("The bundled Rcpp package is required by the C++ permutation engine")
    }
    invisible(TRUE)
}

compile_cpp_engine <- function() {
    cpp_file <- "src/cpp_permutation_engine.cpp"
    if (!file.exists(cpp_file)) {
        stop("C++ source file not found: ", cpp_file)
    }

    cat("Compiling C++ Permutation Engine...\n")
    cache_dir <- file.path(
        tryCatch(
            tools::R_user_dir("ATOP-SCREEN", which = "cache"),
            error = function(e) file.path(tempdir(), "ATOP-SCREEN")
        ),
        "cpp-engine-2026-07-20-o3"
    )
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    r_command <- file.path(
        R.home("bin"),
        if (.Platform$OS.type == "windows") "R.exe" else "R"
    )
    base_flags <- tryCatch(
        system2(r_command, c("CMD", "config", "CXX17FLAGS"), stdout = TRUE, stderr = FALSE),
        error = function(e) character()
    )
    if (length(base_flags) > 0L && nzchar(paste(base_flags, collapse = " "))) {
        makevars_file <- tempfile("atop-cpp17-makevars-")
        writeLines(
            paste("CXX17FLAGS =", paste(base_flags, collapse = " "), "-O3 -DNDEBUG"),
            makevars_file
        )
        old_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
        Sys.setenv(R_MAKEVARS_USER = makevars_file)
        on.exit({
            if (is.na(old_makevars)) {
                Sys.unsetenv("R_MAKEVARS_USER")
            } else {
                Sys.setenv(R_MAKEVARS_USER = old_makevars)
            }
            unlink(makevars_file)
        }, add = TRUE)
    }

    tryCatch(
        {
            suppressMessages(
                Rcpp::sourceCpp(
                    cpp_file,
                    cacheDir = cache_dir,
                    rebuild = FALSE,
                    showOutput = FALSE,
                    verbose = FALSE
                )
            )
            if (!exists("perform_cpp_permutation_test", mode = "function")) {
                stop("C++ function is unavailable after compilation")
            }
            cat("C++ Engine compiled successfully\n")
            TRUE
        },
        error = function(e) {
            cat("C++ Engine compilation failed:", conditionMessage(e), "\n")
            FALSE
        }
    )
}

load_precompiled_cpp_engine <- function() {
    library_name <- if (.Platform$OS.type == "windows") {
        "atop_permutation_engine.dll"
    } else {
        "atop_permutation_engine.so"
    }
    candidates <- unique(c(
        Sys.getenv("ATOP_CPP_ENGINE_LIBRARY", unset = ""),
        file.path(getwd(), "lib", library_name)
    ))
    candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
    if (length(candidates) == 0L) {
        return(FALSE)
    }

    for (candidate in candidates) {
        loaded <- tryCatch(
            {
                dll_info <- dyn.load(normalizePath(candidate, mustWork = TRUE))
                calculate_function <- Rcpp:::sourceCppFunction(
                    function(diff_scores, gene_labels, unique_genes, min_sgrna_threshold = 3L) {},
                    FALSE,
                    dll_info,
                    "sourceCpp_1_calculate_observed_scores_cpp"
                )
                permutation_function <- Rcpp:::sourceCppFunction(
                    function(
                        diff_scores,
                        gene_labels,
                        unique_genes,
                        n_permutations = 1000L,
                        min_sgrna_threshold = 3L,
                        seed = 42L,
                        show_progress = TRUE,
                        n_threads = 0L,
                        permutation_offset = 0L
                    ) {},
                    FALSE,
                    dll_info,
                    "sourceCpp_1_perform_cpp_permutation_test"
                )
                target_environment <- environment(initialize_cpp_engine)
                ncs_function <- Rcpp:::sourceCppFunction(
                    function(lfc, delta0, delta1, delta2) {}, FALSE, dll_info,
                    "sourceCpp_1_calculate_ncs_cpp"
                )
                assign("calculate_ncs_cpp", ncs_function, envir = target_environment)
                assign("calculate_observed_scores_cpp", calculate_function, envir = target_environment)
                assign("perform_cpp_permutation_test", permutation_function, envir = target_environment)
                TRUE
            },
            error = function(e) FALSE
        )
        if (isTRUE(loaded)) {
            cat("Precompiled C++ Engine loaded successfully\n")
            return(TRUE)
        }
    }
    FALSE
}

detect_system_config <- function() {
    detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
    if (length(detected) != 1L || is.na(detected) || detected < 1L) {
        detected <- 1L
    }
    detected <- as.integer(detected)
    recommended_threads <- max(1L, detected - 1L)

    memory_gb <- tryCatch(
        {
            if (identical(Sys.info()[["sysname"]], "Darwin")) {
                bytes <- as.numeric(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE))
                round(bytes / 1024^3, 1)
            } else {
                NA_real_
            }
        },
        error = function(e) NA_real_
    )

    list(
        n_cores = detected,
        optimal_threads = recommended_threads,
        memory_gb = memory_gb
    )
}

resolve_cpp_threads <- function(n_threads = NULL, n_permutations = 1L) {
    config <- detect_system_config()
    if (
        is.null(n_threads) ||
        length(n_threads) != 1L ||
        is.na(n_threads) ||
        !is.finite(n_threads) ||
        n_threads <= 0
    ) {
        n_threads <- config$optimal_threads
    }
    as.integer(max(1L, min(as.integer(n_threads), as.integer(n_permutations), config$n_cores)))
}

normalize_cpp_seed <- function(seed = NULL) {
    if (is.null(seed)) {
        return(sample.int(.Machine$integer.max, 1L))
    }
    if (length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
        stop("seed must be one finite integer")
    }
    as.integer(seed)
}

prepare_cpp_permutation_data <- function(
    sgrna_data,
    gene_col = "Gene",
    score_col = "diff_score",
    min_sgrna_threshold = 3L
) {
    if (!gene_col %in% colnames(sgrna_data)) {
        stop("Gene column not found: ", gene_col)
    }
    if (!score_col %in% colnames(sgrna_data)) {
        stop("Score column not found: ", score_col)
    }
    if (length(min_sgrna_threshold) != 1L || is.na(min_sgrna_threshold) || min_sgrna_threshold < 1L) {
        stop("min_sgrna_threshold must be at least 1")
    }

    valid_rows <- !is.na(sgrna_data[[gene_col]]) &
        !is.na(sgrna_data[[score_col]]) &
        is.finite(sgrna_data[[score_col]])
    if (!any(valid_rows)) {
        stop("No valid data rows")
    }

    clean_data <- sgrna_data[valid_rows, , drop = FALSE]
    gene_counts <- table(clean_data[[gene_col]])
    valid_genes <- names(gene_counts)[gene_counts >= min_sgrna_threshold]
    if (length(valid_genes) == 0L) {
        stop("No genes have at least ", min_sgrna_threshold, " valid sgRNAs")
    }

    final_data <- clean_data[clean_data[[gene_col]] %in% valid_genes, , drop = FALSE]
    gene_labels <- as.character(final_data[[gene_col]])
    list(
        diff_scores = as.numeric(final_data[[score_col]]),
        gene_labels = gene_labels,
        unique_genes = sort(unique(gene_labels)),
        n_sgrnas = nrow(final_data)
    )
}

format_cpp_permutation_result <- function(
    cpp_result,
    gene_summary_data = NULL,
    gene_col = "Gene",
    elapsed_time = NA_real_
) {
    p_positive <- cpp_result$P_positive
    p_negative <- cpp_result$P_negative
    gene_names <- cpp_result$gene_names

    if (!is.null(gene_summary_data)) {
        if (!gene_col %in% colnames(gene_summary_data)) {
            stop("Gene column not found in gene summary data: ", gene_col)
        }
        matched <- match(as.character(gene_summary_data[[gene_col]]), gene_names)
        p_positive <- cpp_result$P_positive[matched]
        p_negative <- cpp_result$P_negative[matched]
    }

    elapsed_for_rate <- max(as.numeric(elapsed_time), .Machine$double.eps, na.rm = TRUE)
    list(
        P_positive = p_positive,
        P_negative = p_negative,
        gene_names = if (is.null(gene_summary_data)) gene_names else NULL,
        elapsed_time = as.numeric(elapsed_time),
        n_threads = as.integer(cpp_result$n_threads),
        seed = as.integer(cpp_result$seed),
        engine = "C++17",
        performance_stats = list(
            n_valid_genes = as.integer(cpp_result$n_valid_genes),
            n_permutations = as.integer(cpp_result$n_permutations),
            speed = paste(round(cpp_result$n_permutations / elapsed_for_rate), "perms/sec")
        )
    )
}

perform_cpp_permutation <- function(
    sgrna_data,
    gene_summary_data = NULL,
    gene_col = "Gene",
    score_col = "diff_score",
    n_permutations = 1000L,
    min_sgrna_threshold = 3L,
    n_threads = NULL,
    progress_callback = NULL,
    seed = NULL
) {
    n_permutations <- as.integer(n_permutations)
    if (length(n_permutations) != 1L || is.na(n_permutations) || n_permutations < 1L) {
        stop("n_permutations must be at least 1")
    }

    if (!is.null(progress_callback)) {
        progress_callback("Starting C++17 permutation engine")
        progress_callback("Preprocessing data")
    }

    prepared <- prepare_cpp_permutation_data(
        sgrna_data,
        gene_col,
        score_col,
        min_sgrna_threshold
    )
    resolved_threads <- resolve_cpp_threads(n_threads, n_permutations)
    resolved_seed <- normalize_cpp_seed(seed)

    if (!is.null(progress_callback)) {
        progress_callback(paste0(
            "Executing ", n_permutations, " permutations with ",
            resolved_threads, " C++ thread", if (resolved_threads == 1L) "" else "s"
        ))
    }

    start_time <- proc.time()[["elapsed"]]
    cpp_result <- tryCatch(
        perform_cpp_permutation_test(
            diff_scores = prepared$diff_scores,
            gene_labels = prepared$gene_labels,
            unique_genes = prepared$unique_genes,
            n_permutations = n_permutations,
            min_sgrna_threshold = as.integer(min_sgrna_threshold),
            seed = resolved_seed,
            show_progress = FALSE,
            n_threads = resolved_threads,
            permutation_offset = 0L
        ),
        error = function(e) stop("C++ permutation test failed: ", conditionMessage(e))
    )
    elapsed_time <- proc.time()[["elapsed"]] - start_time

    if (!is.null(progress_callback)) {
        progress_callback(paste0(
            "Permutation test completed in ", round(elapsed_time, 2),
            " seconds using ", cpp_result$n_threads, " C++ thread",
            if (cpp_result$n_threads == 1L) "" else "s"
        ))
    }

    format_cpp_permutation_result(
        cpp_result,
        gene_summary_data,
        gene_col,
        elapsed_time
    )
}

initialize_cpp_engine <- function() {
    cat("Initializing C++ Engine...\n")
    if (exists("perform_cpp_permutation_test", mode = "function")) {
        cat("C++ Engine available, skipping compilation\n")
        return(TRUE)
    }
    if (!isTRUE(tryCatch(setup_cpp_environment(), error = function(e) FALSE))) {
        return(FALSE)
    }
    if (isTRUE(load_precompiled_cpp_engine())) {
        return(TRUE)
    }
    if (!identical(Sys.getenv("ATOP_ALLOW_CPP_COMPILE", unset = "0"), "1")) {
        cat("Precompiled C++ Engine unavailable; runtime compilation is disabled\n")
        return(FALSE)
    }
    compile_cpp_engine()
}
