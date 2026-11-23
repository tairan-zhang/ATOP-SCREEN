health_result <- function() {
    source_modules(c("analysis", "gsea", "plotting"))
    packages <- list(
        jsonlite = TRUE,
        dplyr = requireNamespace("dplyr", quietly = TRUE),
        readxl = requireNamespace("readxl", quietly = TRUE),
        Rcpp = requireNamespace("Rcpp", quietly = TRUE),
        clusterProfiler = requireNamespace("clusterProfiler", quietly = TRUE),
        ggplot2 = requireNamespace("ggplot2", quietly = TRUE),
        ggpubr = requireNamespace("ggpubr", quietly = TRUE),
        ggrepel = requireNamespace("ggrepel", quietly = TRUE),
        patchwork = requireNamespace("patchwork", quietly = TRUE),
        svglite = requireNamespace("svglite", quietly = TRUE)
    )
    list(
        status = "ready",
        analysisRoot = analysis_root,
        runtimeRoot = runtime_root,
        rVersion = R.version.string,
        loadedModules = state$loaded,
        functions = list(
            process = exists("perform_crispr_screen_analysis", mode = "function"),
            gsea = exists("perform_gsea_analysis_clusterprofiler", mode = "function"),
            plots = exists("generate_volcano_plot", mode = "function")
        ),
        packages = packages
    )
}

inspect_file <- function(request) {
    preview <- read_screen_file(request$path, preview_rows = 30L)
    list(
        path = normalizePath(request$path, mustWork = TRUE),
        fileName = basename(request$path),
        columns = names(preview),
        preview = preview,
        previewRows = nrow(preview),
        fileSize = as.numeric(file.info(request$path)$size)
    )
}

analysis_defaults <- function() {
    list(
        gRNA_col = "gRNA",
        gene_col = "Gene",
        sequence_col = NULL,
        condition1_replicate_cols = character(),
        condition2_replicate_cols = character(),
        N_perm = 100000L,
        pseudo_count_lfc = 1,
        min_sgrna_threshold = 3L,
        user_engine_choice = "cpp",
        skip_normalization = FALSE,
        diff_score_col1 = NULL,
        diff_score_col2 = NULL
    )
}

merge_defaults <- function(values, defaults) {
    values <- values %||% list()
    for (key in names(defaults)) if (is.null(values[[key]])) values[[key]] <- defaults[[key]]
    values
}

process_data <- function(request) {
    source_modules("analysis")
    params <- merge_defaults(request$parameters, analysis_defaults())
    with_mageck <- isTRUE(params$run_mageck)
    params$run_mageck <- NULL
    raw <- read_screen_file(request$path)
    if (identical(params$user_engine_choice, "cpp") && exists("initialize_cpp_engine", mode = "function")) {
        try(with_captured_progress(initialize_cpp_engine()), silent = TRUE)
    }
    result <- with_captured_progress(do.call(
        perform_crispr_screen_analysis,
        c(list(raw_data_df = raw), params)
    ))
    mageck <- if (with_mageck) run_mageck_rra(raw, result$params, request$artifactDir)
    result$params$run_mageck <- with_mageck
    result$params$mageck <- mageck$settings
    result$params$r_version <- as.character(getRversion())
    result$params$software_version <- read.dcf(file.path(analysis_root, "DESCRIPTION"), fields = "Version")[[1L]]
    result$params$input_sha256 <- digest::digest(file = request$path, algo = "sha256", serialize = FALSE)
    lock_path <- if (file.exists(file.path(runtime_root, "renv.lock"))) file.path(runtime_root, "renv.lock") else file.path(analysis_root, "..", "..", "renv.lock")
    if (file.exists(lock_path)) result$params$renv_lock_sha256 <- digest::digest(file = lock_path, algo = "sha256", serialize = FALSE)
    state$data_processing <- result
    state$last_input_path <- normalizePath(request$path, mustWork = TRUE)
    artifacts <- list(
        normalizedCountsRds = save_rds_artifact(result$normalized_counts, request$artifactDir, "normalized_counts.rds"),
        analysisStateRds = save_rds_artifact(result, request$artifactDir, "analysis_state.rds")
    )
    artifacts <- c(artifacts, mageck$files)
    if (!is.null(request$artifactDir)) {
        target <- file.path(request$artifactDir, "analysis-settings.json")
        jsonlite::write_json(result$params, target, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
        artifacts$settings <- target
    }
    list(
        status = "completed",
        summary = list(
            processedSgrnaRows = nrow(result$processed_sg_data),
            geneRows = nrow(result$gene_summary_data),
            filteredRows = if (is.null(result$filtered_genes_sgrna_data)) 0L else nrow(result$filtered_genes_sgrna_data),
            parameters = result$params
        ),
        tables = list(
            sgrna = compact_table(result$processed_sg_data),
            ncsExcluded = compact_table(result$ncs_excluded),
            mageckGenes = compact_table(mageck$genes),
            mageckGuides = compact_table(mageck$guides),
            ncsQc = compact_table(result$sg_ncs_filter_qc),
            genes = compact_table(result$gene_summary_data),
            filtered = compact_table(result$filtered_genes_sgrna_data)
        ),
        artifacts = artifacts[!vapply(artifacts, is.null, logical(1))]
    )
}

gsea_defaults <- function() {
    list(
        clp_n_permutations = 1000L,
        clp_pvalue_cutoff = 1,
        clp_min_gene_set_size = 15L,
        max_gene_set_size = 500L,
        clp_seed_logical = TRUE,
        permutation_seed = 123L,
        clp_by = "multilevel"
    )
}

run_gsea <- function(request) {
    source_modules("gsea")
    if (is.null(state$data_processing)) stop("Complete data processing before running GSEA")
    if (is.null(request$gmtPath) || !file.exists(request$gmtPath)) stop("Select a valid GMT file before running GSEA")
    params <- merge_defaults(request$parameters, gsea_defaults())
    result <- with_captured_progress(do.call(
        perform_gsea_analysis_clusterprofiler,
        c(
            list(
                gene_summary_df = state$data_processing$gene_summary_data,
                gmt_file_path = request$gmtPath,
                gene_id_column = state$data_processing$params$gene_col
            ),
            params
        )
    ))
    params$method_used <- if (identical(params$clp_by, "permute")) "permute" else "multilevel"
    params$permutation_count_used <- if (params$method_used == "permute") params$clp_n_permutations else NULL
    params$multilevel_simple_permutations <- if (params$method_used == "multilevel") 1000L else NULL
    state$gsea <- result
    artifacts <- list(
        gseaPositiveRds = save_rds_artifact(result$gsea_object_positive, request$artifactDir, "gsea_positive.rds"),
        gseaNegativeRds = save_rds_artifact(result$gsea_object_negative, request$artifactDir, "gsea_negative.rds"),
        gseaStateRds = save_rds_artifact(result, request$artifactDir, "gsea_state.rds")
    )
    list(
        status = "completed",
        parameters = params,
        positive = compact_table(result$gsea_results_positive_df),
        negative = compact_table(result$gsea_results_negative_df),
        artifacts = artifacts[!vapply(artifacts, is.null, logical(1))]
    )
}

restore_project_state <- function(request) {
    params <- request$parameters %||% list()
    analysis_path <- params$analysisStateRds %||% NULL
    gsea_path <- params$gseaStateRds %||% NULL
    if (is.null(analysis_path) || !file.exists(analysis_path)) stop("The saved ATOP analysis state is unavailable")
    state$data_processing <- readRDS(normalizePath(analysis_path, mustWork = TRUE))
    state$gsea <- if (!is.null(gsea_path) && file.exists(gsea_path)) readRDS(normalizePath(gsea_path, mustWork = TRUE)) else NULL
    list(status = "restored", analysis = TRUE, gsea = !is.null(state$gsea))
}
