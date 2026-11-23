#!/usr/bin/env Rscript
# Restore the exact source archives selected in renv.lock into a dedicated library.
args <- commandArgs(trailingOnly = TRUE)
library_arg <- grep("^--library=", args, value = TRUE)
if (length(library_arg) != 1L) stop("Supply --library=/absolute/path/to/a/dedicated/R/library")
dest <- sub("^--library=", "", library_arg)
if (!grepl("^(/|[A-Za-z]:[/\\\\])", dest)) stop("Use an absolute library path.")
if (as.character(getRversion()) != "4.6.1") stop("The lockfile requires R 4.6.1.")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
dest <- normalizePath(dest, mustWork = TRUE)
.libPaths(c(dest, .Library), include.site = FALSE)
Sys.setenv(R_LIBS = dest, R_LIBS_USER = dest, R_LIBS_SITE = "")
options(repos = NULL)
sources <- normalizePath("packaging/sources", mustWork = TRUE)
sums <- utils::read.table(file.path(sources, "SHA256SUMS"), col.names = c("sha", "file"),
                          colClasses = "character", comment.char = "", stringsAsFactors = FALSE)
paths <- file.path(sources, sums$file)
bad <- which(!file.exists(paths) | tolower(tools::sha256sum(paths)) != sums$sha)
if (length(bad)) stop("Missing or modified source archive: ", paste(sums$file[bad], collapse = ", "))
archive <- function(record) {
    name <- paste0(record$Package, "_", record$Version, ".tar.gz")
    hit <- which(basename(paths) == name)
    if (length(hit) != 1L) stop("Missing or ambiguous pinned archive: ", name)
    paths[[hit]]
}
install_source <- function(path) {
    name <- basename(path)
    log <- tempfile(pattern = "atop-r-install-", fileext = ".log")
    status <- system2(file.path(R.home("bin"), "R"),
                      c("CMD", "INSTALL", "-l", shQuote(dest), shQuote(path)),
                      stdout = log, stderr = log)
    if (status != 0L) {
        cat(readLines(log, warn = FALSE), sep = "\n")
        stop("Installation failed for ", name, call. = FALSE)
    }
    unlink(log)
    message("Installed ", name)
    TRUE
}
installed <- function(pkg, library) {
    path <- file.path(library, pkg, "DESCRIPTION")
    if (file.exists(path)) read.dcf(path, fields = "Version")[[1L]] else ""
}
if (installed("renv", dest) != "1.2.1")
    install_source(archive(list(Package = "renv", Version = "1.2.1")))
invisible(loadNamespace("renv", lib.loc = dest))
records <- renv::lockfile_read("renv.lock")$Packages
# Fail before compiling if the release omits a locked source or dependency.
invisible(lapply(records, archive))
unknown <- setdiff(unique(unlist(lapply(records, "[[", "Requirements"))), names(records))
if (length(unknown)) stop("Incomplete lockfile: ", paste(unknown, collapse = ", "))
matches <- function(record) {
    priority <- if (file.exists(file.path(.Library, record$Package, "DESCRIPTION")))
        read.dcf(file.path(.Library, record$Package, "DESCRIPTION"), fields = "Priority")[[1L]]
    actual <- installed(record$Package, dest)
    if (!nzchar(actual) && identical(priority, "recommended")) actual <- installed(record$Package, .Library)
    identical(actual, record$Version)
}
pending <- names(Filter(Negate(matches), records))
jobs <- if (.Platform$OS.type == "windows") 1L else as.integer(Sys.getenv("ATOP_BUILD_JOBS", "2"))
if (is.na(jobs) || jobs < 1L) stop("ATOP_BUILD_JOBS must be a positive integer.")
while (length(pending)) {
    ready <- pending[vapply(pending, function(pkg) !any(records[[pkg]]$Requirements %in% pending), logical(1))]
    if (!length(ready)) stop("Circular dependencies in renv.lock: ", paste(pending, collapse = ", "))
    results <- parallel::mclapply(ready, function(pkg) install_source(archive(records[[pkg]])),
                                  mc.cores = min(jobs, length(ready)), mc.preschedule = FALSE)
    if (!all(vapply(results, identical, logical(1), TRUE)))
        stop("Could not restore ", paste(ready[!vapply(results, identical, logical(1), TRUE)], collapse = ", "))
    pending <- setdiff(pending, ready)
}
if (!all(vapply(records, matches, logical(1)))) stop("Restored versions differ from renv.lock.")
cat("Verified", length(records), "locked R packages.\n")
