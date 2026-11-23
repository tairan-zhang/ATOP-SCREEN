#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)

argument_value <- function(name, default = NULL) {
    prefix <- paste0("--", name, "=")
    match <- arguments[startsWith(arguments, prefix)]
    if (length(match) == 0L) return(default)
    sub(prefix, "", match[[1L]], fixed = TRUE)
}

script_arguments <- commandArgs(trailingOnly = FALSE)
script_path <- sub(
    "^--file=",
    "",
    script_arguments[grepl("^--file=", script_arguments)][1L]
)
if (is.na(script_path) || !nzchar(script_path)) stop("Run this file with Rscript")

project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source_path <- file.path(project_root, "packages", "analysis-engine", "src", "cpp_permutation_engine.cpp")
output_path <- argument_value(
    "output",
    file.path(project_root, "packaging", "runtime", "macos-arm64", "analysis", "lib", "atop_permutation_engine.so")
)
if (!grepl("^/", output_path)) output_path <- file.path(project_root, output_path)

if (!file.exists(source_path)) stop("C++ engine source not found: ", source_path)
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required to build the C++ engine")

build_root <- tempfile("atop-cpp-engine-")
dir.create(build_root, recursive = TRUE)
on.exit(unlink(build_root, recursive = TRUE, force = TRUE), add = TRUE)

makevars_path <- file.path(build_root, "Makevars")
base_flags <- system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX17FLAGS"), stdout = TRUE)
base_flags <- base_flags[!startsWith(base_flags, "WARNING:")]
writeLines(
    paste("CXX17FLAGS =", paste(base_flags, collapse = " "), "-O3 -DNDEBUG"),
    makevars_path
)
old_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
Sys.setenv(R_MAKEVARS_USER = makevars_path)
on.exit({
    if (is.na(old_makevars)) Sys.unsetenv("R_MAKEVARS_USER") else Sys.setenv(R_MAKEVARS_USER = old_makevars)
}, add = TRUE)

loaded_before <- vapply(getLoadedDLLs(), `[[`, character(1), "path")
Rcpp::sourceCpp(
    source_path,
    cacheDir = file.path(build_root, "cache"),
    rebuild = TRUE,
    showOutput = TRUE,
    verbose = FALSE
)
loaded_after <- vapply(getLoadedDLLs(), `[[`, character(1), "path")
engine_candidates <- setdiff(loaded_after, loaded_before)
engine_candidates <- engine_candidates[file.exists(engine_candidates)]
engine_candidates <- engine_candidates[grepl("sourceCpp.*[.]so$", basename(engine_candidates))]
if (length(engine_candidates) != 1L) {
    stop("Unable to identify the compiled C++ engine library")
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
if (!file.copy(engine_candidates[[1L]], output_path, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Failed to copy the C++ engine to: ", output_path)
}
cat("Built ATOP C++ engine with R ", as.character(getRversion()), ": ", normalizePath(output_path), "\n", sep = "")
