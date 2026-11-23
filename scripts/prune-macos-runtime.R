#!/usr/bin/env Rscript

command_arguments <- commandArgs(trailingOnly = TRUE)
apply_changes <- "--apply" %in% command_arguments

argument_value <- function(name, default = NULL) {
    prefix <- paste0("--", name, "=")
    match <- command_arguments[startsWith(command_arguments, prefix)]
    if (length(match) == 0L) return(default)
    sub(prefix, "", match[[1L]], fixed = TRUE)
}

script_argument <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_argument[grepl("^--file=", script_argument)][1])
if (is.na(script_path) || !nzchar(script_path)) stop("Run this file with Rscript")

project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
manifest_path <- file.path(project_root, "packaging", "runtime-manifest.json")
if (!file.exists(manifest_path)) stop("Runtime manifest not found: ", manifest_path)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required to read the runtime manifest")

manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
runtime_override <- argument_value("runtime-root")
runtime_root <- if (is.null(runtime_override)) {
    normalizePath(file.path(project_root, manifest$runtimeRoot), mustWork = TRUE)
} else {
    normalizePath(runtime_override, mustWork = TRUE)
}
library_root <- normalizePath(file.path(runtime_root, manifest$libraryPath), mustWork = TRUE)
resources_root <- dirname(library_root)

if (basename(library_root) != "library" || !file.exists(file.path(library_root, "base", "DESCRIPTION"))) {
    stop("Refusing to prune an unrecognised R library: ", library_root)
}

installed <- as.data.frame(utils::installed.packages(lib.loc = library_root), stringsAsFactors = FALSE)
installed_packages <- unique(installed$Package)
direct_packages <- unique(as.character(manifest$directPackages))
missing_direct <- setdiff(direct_packages, installed_packages)
if (length(missing_direct) > 0L) {
    stop("Required runtime packages are missing: ", paste(missing_direct, collapse = ", "))
}

dependency_map <- tools::package_dependencies(
    packages = direct_packages,
    db = installed,
    which = c("Depends", "Imports", "LinkingTo"),
    recursive = TRUE
)
recommended_packages <- rownames(installed)[installed$Priority %in% c("base", "recommended")]
retained_packages <- sort(unique(intersect(
    c(direct_packages, unlist(dependency_map, use.names = FALSE), recommended_packages),
    installed_packages
)))
removed_packages <- sort(setdiff(installed_packages, retained_packages))

directory_size <- function(path) {
    files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, include.dirs = FALSE)
    sum(file.info(files)$size, na.rm = TRUE)
}

removed_bytes <- sum(vapply(file.path(library_root, removed_packages), directory_size, numeric(1)))
cat(sprintf("Runtime package audit: retain %d of %d packages; remove %d packages (%.1f MiB).\n",
    length(retained_packages), length(installed_packages), length(removed_packages), removed_bytes / 1024^2))
cat("Packages outside the declared closure:\n", paste(removed_packages, collapse = ", "), "\n", sep = "")

if (!apply_changes) {
    cat("Dry run only. Re-run with --apply to prune this bundled runtime.\n")
    quit(status = 0L)
}

for (package_name in removed_packages) {
    unlink(file.path(library_root, package_name), recursive = TRUE, force = TRUE)
}

if (isTRUE(manifest$prune$removeRCoreTests)) {
    unlink(file.path(resources_root, "tests"), recursive = TRUE, force = TRUE)
}
if (isTRUE(manifest$prune$removeDebugSymbols)) {
    debug_symbol_dirs <- list.dirs(resources_root, recursive = TRUE, full.names = TRUE)
    debug_symbol_dirs <- debug_symbol_dirs[grepl("\\.dSYM$", debug_symbol_dirs)]
    for (directory in debug_symbol_dirs) unlink(directory, recursive = TRUE, force = TRUE)
}

remove_paths <- unique(as.character(manifest$prune$removePaths))
remove_paths <- remove_paths[nzchar(remove_paths)]
for (relative_path in remove_paths) {
    if (grepl("^/", relative_path) || grepl("(^|/)\\.\\.(/|$)", relative_path)) {
        stop("Refusing unsafe runtime prune path: ", relative_path)
    }
    target_path <- file.path(runtime_root, relative_path)
    target_parent <- normalizePath(dirname(target_path), mustWork = TRUE)
    runtime_prefix <- paste0(runtime_root, .Platform$file.sep)
    if (!startsWith(paste0(target_parent, .Platform$file.sep), runtime_prefix)) {
        stop("Runtime prune path escapes the runtime root: ", relative_path)
    }
    if (file.exists(target_path) || dir.exists(target_path)) {
        unlink(target_path, recursive = TRUE, force = TRUE)
        cat("Removed unused runtime path: ", relative_path, "\n", sep = "")
    }
}

remaining <- sort(utils::installed.packages(lib.loc = library_root)[, "Package"])
missing_retained <- setdiff(retained_packages, remaining)
unexpected_remaining <- setdiff(remaining, retained_packages)
if (length(missing_retained) > 0L || length(unexpected_remaining) > 0L) {
    stop(
        "Runtime package validation failed after pruning. Missing: ",
        paste(missing_retained, collapse = ", "),
        "; unexpected: ",
        paste(unexpected_remaining, collapse = ", ")
    )
}
writeLines(remaining, file.path(resources_root, "ATOP_SCREEN_R_PACKAGE_MANIFEST.txt"))
remaining_installed <- as.data.frame(
    utils::installed.packages(lib.loc = library_root),
    stringsAsFactors = FALSE
)
manifest_columns <- c("Package", "Version", "Built", "Repository")
for (column_name in setdiff(manifest_columns, colnames(remaining_installed))) {
    remaining_installed[[column_name]] <- ""
}
remaining_metadata <- remaining_installed[remaining, manifest_columns, drop = FALSE]
utils::write.csv(
    remaining_metadata,
    file.path(resources_root, "ATOP_SCREEN_R_PACKAGE_VERSIONS.csv"),
    row.names = FALSE,
    na = ""
)
cat(sprintf("Runtime pruned successfully. Remaining package count: %d.\n", length(remaining)))
