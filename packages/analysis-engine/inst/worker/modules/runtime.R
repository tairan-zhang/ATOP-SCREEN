suppressPackageStartupMessages(library(jsonlite))

analysis_root <- normalizePath(enc2utf8(Sys.getenv("ATOP_ANALYSIS_ROOT", unset = getwd())), mustWork = FALSE)
runtime_root <- enc2utf8(Sys.getenv("ATOP_RUNTIME_ROOT", unset = ""))
state <- new.env(parent = emptyenv())
state$loaded <- character()
state$data_processing <- NULL
state$gsea <- NULL
state$last_input_path <- NULL

`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0L) rhs else lhs

emit <- function(payload) {
    cat(jsonlite::toJSON(payload, digits = NA, auto_unbox = TRUE, null = "null", na = "null"), "\n", sep = "")
    flush.console()
}

emit_error <- function(id, error) {
    emit(list(id = id, ok = FALSE, error = conditionMessage(error), class = class(error)[[1]] %||% "error"))
}

safe_library <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(FALSE)
    suppressPackageStartupMessages(library(package, character.only = TRUE))
    TRUE
}

load_common_packages <- function() {
    packages <- c(
        "dplyr", "tidyr", "stringr", "data.table", "readxl", "Rcpp",
        "ggplot2", "ggpubr", "ggrepel", "patchwork", "RColorBrewer",
        "svglite", "ragg"
    )
    loaded <- setNames(logical(length(packages)), packages)
    for (index in seq_along(packages)) {
        loaded[[index]] <- safe_library(packages[[index]])
    }
    invisible(loaded)
}

module_paths <- list(
    progress = file.path("R", "runtime", "progress_tracker.R"),
    r_permutation = file.path("R", "permutation", "r_parallel.R"),
    cpp_interface = file.path("R", "permutation", "cpp", "interface.R"),
    engine_selector = file.path("R", "permutation", "engine_selector.R"),
    atop_validation = file.path("R", "atop", "validation.R"),
    atop_preprocessing = file.path("R", "atop", "preprocessing.R"),
    atop_scoring = file.path("R", "atop", "scoring.R"),
    atop_statistics = file.path("R", "atop", "statistics.R"),
    atop_orchestrator = file.path("R", "atop", "orchestrator.R"),
    gsea_clusterprofiler = file.path("R", "gsea", "clusterprofiler.R"),
    plotting_typography = file.path("R", "plotting", "common", "typography.R"),
    gene_common = file.path("R", "plotting", "gene", "common.R"),
    gene_volcano = file.path("R", "plotting", "gene", "volcano.R"),
    gene_ranking = file.path("R", "plotting", "gene", "ranking.R"),
    gsea_lollipop = file.path("R", "plotting", "gsea", "lollipop.R"),
    gsea_running_helpers = file.path("R", "plotting", "gsea", "running_helpers.R"),
    gsea_running_gseavis = file.path("R", "plotting", "gsea", "running_gseavis.R"),
    gsea_running_style = file.path("R", "plotting", "gsea", "running_style.R"),
    gsea_running = file.path("R", "plotting", "gsea", "running.R"),
    sgrna_paired = file.path("R", "plotting", "sgrna", "paired.R")
)

source_modules <- function(groups = "analysis") {
    load_common_packages()
    requested <- unique(unlist(lapply(groups, function(group) {
        switch(group,
            analysis = c("progress", "r_permutation", "cpp_interface", "engine_selector", "atop_validation", "atop_preprocessing", "atop_scoring", "atop_statistics", "atop_orchestrator"),
            gsea = "gsea_clusterprofiler",
            plotting = c(
                "plotting_typography", "gene_common", "gene_volcano", "gene_ranking", "gsea_lollipop",
                "gsea_running_helpers", "gsea_running_gseavis", "gsea_running_style", "gsea_running",
                "sgrna_paired"
            ),
            character()
        )
    })))
    for (name in requested) {
        if (name %in% state$loaded) next
        file <- file.path(analysis_root, module_paths[[name]])
        if (!file.exists(file)) stop("Analysis module not found: ", file)
        sys.source(file, envir = .GlobalEnv)
        state$loaded <- c(state$loaded, name)
    }
    invisible(TRUE)
}

read_screen_file <- function(path, preview_rows = NULL) {
    if (!is.character(path) || length(path) != 1L || !nzchar(path)) stop("A data file path is required")
    path <- enc2utf8(path)
    path <- normalizePath(path, mustWork = TRUE)
    extension <- tolower(tools::file_ext(path))
    n_max <- if (is.null(preview_rows)) NULL else as.integer(preview_rows)
    read_delimited <- function(reader, ...) {
        args <- list(...)
        if (!is.null(n_max)) args$nrows <- n_max
        do.call(reader, args)
    }
    if (extension == "csv") return(read_delimited(read.csv, file = path, stringsAsFactors = FALSE, check.names = FALSE))
    if (extension %in% c("txt", "tsv")) return(read_delimited(read.delim, file = path, stringsAsFactors = FALSE, check.names = FALSE))
    if (extension == "xlsx") {
        if (!requireNamespace("readxl", quietly = TRUE)) stop("The readxl R package is required for XLSX files")
        if (is.null(n_max)) return(readxl::read_excel(path, .name_repair = "minimal"))
        return(readxl::read_excel(path, .name_repair = "minimal", n_max = n_max))
    }
    stop("Unsupported file type: .", extension)
}

compact_table <- function(value, rows = NULL) {
    if (is.null(value)) return(NULL)
    frame <- as.data.frame(value, stringsAsFactors = FALSE)
    if (!is.null(rows) && nrow(frame) > rows) frame <- frame[seq_len(rows), , drop = FALSE]
    frame
}

save_rds_artifact <- function(value, directory, file_name) {
    if (is.null(value) || is.null(directory) || !nzchar(as.character(directory))) return(NULL)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    target <- file.path(directory, file_name)
    saveRDS(value, target)
    if (!file.exists(target)) stop("Failed to persist analysis artifact: ", file_name)
    normalizePath(target, mustWork = TRUE)
}

with_captured_progress <- function(expression) {
    result <- NULL
    output <- capture.output(result <- force(expression), type = "output")
    for (line in output[nzchar(trimws(output))]) {
        cat(jsonlite::toJSON(list(event = "progress", message = line), auto_unbox = TRUE), "\n", sep = "", file = stderr())
    }
    result
}
