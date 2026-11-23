# Run upstream MAGeCK on raw counts, independently of ATOP normalization and NCS.
run_mageck_rra <- function(raw, params, directory) {
    if (isTRUE(params$skip_normalization)) stop("MAGeCK RRA requires raw counts.")
    endpoint <- params$condition1_replicate_cols
    baseline <- params$condition2_replicate_cols
    columns <- c(endpoint, baseline)
    counts <- raw[columns]
    if (!all(vapply(counts, function(x) is.numeric(x) && all(is.finite(x) & x >= 0 & x == floor(x)), logical(1)))) {
        stop("MAGeCK requires nonnegative integer counts for every selected sample.")
    }
    ids <- raw[[params$gRNA_col]]
    genes <- raw[[params$gene_col]]
    if (anyNA(ids) || anyNA(genes) || any(!nzchar(ids)) || any(!nzchar(genes)) || anyDuplicated(ids) ||
        any(grepl("[[:space:]]", c(ids, genes)))) stop("MAGeCK requires unique guide IDs and whitespace-free guide and gene IDs.")
    if (is.null(directory)) directory <- tempfile("atop-mageck-")
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    directory <- normalizePath(directory, mustWork = TRUE)
    bin <- Sys.getenv("ATOP_MAGECK_BIN", file.path(runtime_root, "mageck", "bin"))
    python <- file.path(bin, "python")
    cli <- file.path(bin, "mageck")
    if (!all(file.exists(c(python, cli, file.path(bin, "RRA"))))) stop("Bundled MAGeCK runtime is missing.")
    # MAGeCK starts its native RRA helper by name. Keep generated filenames and
    # sample arguments independent of input paths or column names containing spaces.
    previous_dir <- setwd(directory)
    previous_path <- Sys.getenv("PATH")
    on.exit({setwd(previous_dir); Sys.setenv(PATH = previous_path)}, add = TRUE)
    Sys.setenv(PATH = paste(bin, previous_path, sep = .Platform$path.sep))
    samples <- paste0("sample", seq_along(columns))
    names(counts) <- samples
    utils::write.table(data.frame(sgRNA = ids, Gene = genes, counts, check.names = FALSE),
                       "mageck-input.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    run <- function(args) system2(python, c("-IB", shQuote(cli), args), stdout = TRUE, stderr = TRUE)
    version <- run("--version")
    if (!identical(trimws(tail(version, 1)), "0.5.9.5")) stop("Expected MAGeCK 0.5.9.5.")
    response <- run(c("test", "-k", "mageck-input.tsv", "-t", paste(head(samples, length(endpoint)), collapse = ","),
                      "-c", paste(tail(samples, length(baseline)), collapse = ","), "-n", "mageck"))
    if (!is.null(attr(response, "status"))) stop("MAGeCK failed: ", paste(tail(response, 8), collapse = "\n"))
    files <- list(mageckGenes = "mageck.gene_summary.txt", mageckGuides = "mageck.sgrna_summary.txt", mageckLog = "mageck.log")
    # Preserve gene and sgRNA identifiers even when they look numeric.
    read <- function(file) {
        labels <- strsplit(readLines(file, n = 1L), "\t", fixed = TRUE)[[1]]
        classes <- rep(NA_character_, length(labels))
        classes[tolower(labels) %in% c("id", "sgrna", "gene")] <- "character"
        utils::read.delim(file, check.names = FALSE, colClasses = classes, stringsAsFactors = FALSE)
    }
    list(genes = read(files$mageckGenes), guides = read(files$mageckGuides),
         files = lapply(files, normalizePath, mustWork = TRUE),
         settings = list(version = "0.5.9.5", method = "test (RRA)", normalization = "median",
                         rra_seed = 123456L, rra_permutation_rounds = 100L,
                         samples = stats::setNames(as.list(columns), samples)))
}
