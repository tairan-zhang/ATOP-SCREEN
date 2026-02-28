# SPDX-License-Identifier: GPL-3.0-or-later
# install_dependencies.R
# Part of ATOP-SCREEN
#
# Copyright (C) 2026 Tairan Zhang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# --- CRISPR Screen Analysis System - Dependency Auto-Installer ---

cat("CRISPR Screen Analysis System - Dependency Installer\n")
cat("===========================================\n\n")

# --- Zero-Config Bootstrap ---
# 1. Set CRAN Mirror if not set (Prevents interactive hang)
if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
    cat("  [BOOT] Setting default CRAN mirror (cloud.r-project.org)...\n")
    r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org/"
    options(repos = r)
}

# 2. Helper for silent package checking
is_installed <- function(pkg) {
    return(requireNamespace(pkg, quietly = TRUE))
}

# Check R Version
r_version <- R.Version()$version.string
cat("Current R Version:", r_version, "\n")

if (as.numeric(R.Version()$major) < 4) {
    stop("Requires R version 4.0.0 or higher, current version is too old")
}

cat("R version check passed\n\n")

# --- Helper Functions ---
log_info <- function(msg) cat(sprintf("  [INFO] %s\n", msg))
log_ok <- function(pkg) cat(sprintf("  [OK] %s\n", pkg))
log_miss <- function(pkg) cat(sprintf("  [..] %s missing. Installing...\n", pkg))
log_fail <- function(pkg, err = "") cat(sprintf("  [!!] %s FAILED. %s\n", pkg, err))

install_if_missing <- function(packages, install_func = install.packages, desc = "Package") {
    success_count <- 0
    total_count <- length(packages)

    cat(sprintf("> Checking %s dependencies (%d total)...\n", desc, total_count))

    for (pkg in packages) {
        if (!is_installed(pkg)) {
            log_miss(pkg)
            tryCatch(
                {
                    install_func(pkg, dependencies = TRUE)
                    if (is_installed(pkg)) {
                        log_ok(pkg)
                        success_count <- success_count + 1
                    } else {
                        log_fail(pkg, "Verification failed")
                    }
                },
                error = function(e) {
                    log_fail(pkg, e$message)
                }
            )
        } else {
            log_ok(pkg)
            success_count <- success_count + 1
        }
    }

    if (success_count == total_count) {
        cat(sprintf("> %s check/install complete. All OK.\n\n", desc))
    } else {
        cat(sprintf("> %s check/install incomplete. (%d/%d OK)\n\n", desc, success_count, total_count))
    }

    return(success_count == total_count)
}

# --- Step 1: Basic CRAN Packages ---
cat("Step 1: Basic CRAN Packages\n")

cran_packages <- c(
    # === Shiny App Framework ===
    "shiny", # Web App Framework
    "shinyjs", # JavaScript Interaction

    # === Data Processing Core ===
    "readxl", # Excel File Reading
    "dplyr", # Data Manipulation Grammar
    "tibble", # Modern Data Frames
    "rlang", # Metaprogramming Support
    "tidyr", # Data Tidy
    "stringr", # String Processing

    # === User Interface ===
    "DT", # Interactive Data Tables

    # === Scientific Plotting ===
    "ggplot2", # Plotting Grammar
    "ggpubr", # Publication Ready Plots
    "ggpmisc", # Plot Statistics
    "ggnewscale", # Multiple Color Scales
    "patchwork", # Plot Layouts
    "RColorBrewer", # Color Palettes
    "ggrepel", # Text Repulsion
    "ggpp", # Grammar Extensions

    # === High Performance Computing ===
    "parallel", # Parallel Computing
    "data.table", # Fast Data Manipulation

    # === C++ Interface ===
    "Rcpp", # R-C++ Interface
    "RcppParallel", # Parallel C++ Computing

    # === File System ===
    "zip", # Zip File Handling

    # === Async Processing ===
    "later" # Deferred Execution
)

cran_success <- install_if_missing(cran_packages, desc = "CRAN Package")

# --- Step 2: Bioconductor Packages ---
cat("Step 2: Bioconductor Packages\n")

# Robust BiocManager Bootstrap
if (!is_installed("BiocManager")) {
    log_miss("BiocManager")
    install.packages("BiocManager")
    if (!is_installed("BiocManager")) {
        log_fail("BiocManager", "Critical failure: Could not bootstrap Bioconductor.")
    } else {
        log_ok("BiocManager")
    }
}

bioc_packages <- c(
    "clusterProfiler", # Gene Enrichment Analysis
    "enrichplot" # Enrichment Visualization
)

bioc_install_func <- function(pkg, ...) {
    BiocManager::install(pkg, ...)
}

bioc_success <- install_if_missing(bioc_packages, bioc_install_func, "Bioconductor Package")

# --- Step 3: Optional Enhanced Packages ---
cat("Step 3: Optional Enhanced Packages\n")

optional_packages <- c(
    "GseaVis" # Advanced GSEA Visualization
)

optional_success <- TRUE
cat(sprintf("> Checking Optional Packages (%d total)...\n", length(optional_packages)))

for (pkg in optional_packages) {
    if (!is_installed(pkg)) {
        log_miss(pkg)
        tryCatch(
            {
                install.packages(pkg, dependencies = TRUE)
                if (is_installed(pkg)) {
                    log_ok(pkg)
                } else {
                    log_info(paste(pkg, "optional install skipped/failed (non-critical)"))
                }
            },
            error = function(e) {
                log_info(paste(pkg, "optional install failed (non-critical):", e$message))
            }
        )
    } else {
        log_ok(pkg)
    }
}
cat("\n")

# --- Step 4: Validate Critical Components ---
cat("Step 4: Validate Critical Components\n")

critical_tests <- list(
    "Shiny Framework" = "shiny",
    "C++ Interface" = "Rcpp",
    "Bioinformatics" = "clusterProfiler",
    "Scientific Plotting" = "ggplot2",
    "Fast Computing" = "data.table",
    "Parallel Computing" = "parallel"
)

all_critical_ok <- TRUE
for (test_name in names(critical_tests)) {
    pkg <- critical_tests[[test_name]]
    if (is_installed(pkg)) {
        cat(sprintf("  [OK] %-20s (%s)\n", test_name, pkg))
    } else {
        cat(test_name, "validation failed\n")
        all_critical_ok <- FALSE
    }
}

# --- Step 5: System Environment Check ---
cat("\nStep 5: System Environment Check\n")

# Detect CPU Cores
n_cores <- parallel::detectCores()
cat("CPU Cores:", n_cores, "\n")

# Detect Memory (Approx)
if (Sys.info()["sysname"] == "Windows") {
    memory_info <- "Requires extra tool to detect"
} else if (Sys.info()["sysname"] == "Darwin") {
    memory_gb <- tryCatch(
        {
            round(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 1024^3, 1)
        },
        error = function(e) "Unknown"
    )
    memory_info <- paste(memory_gb, "GB")
} else {
    memory_info <- "Linux System"
}
cat("System Memory:", memory_info, "\n")

# Detect C++ Compiler
cpp_available <- tryCatch(
    {
        system("R CMD config CXX", intern = TRUE)
        TRUE
    },
    error = function(e) {
        FALSE
    }
)

if (cpp_available) {
    cat("C++ Compiler: Available\n")
} else {
    cat("C++ Compiler: Not Available (Will use R engine)\n")
}

# --- Installation Summary ---
cat("\nInstallation Summary\n")
cat("====================\n")

if (cran_success && bioc_success && all_critical_ok) {
    cat("  [OK] System Ready. All dependencies installed.\n")
    cat("  [..] Usage: source('app.R')\n")
} else {
    cat("  [!!] System verification FAILED.\n")
    cat("  [..] Please check network connection and try again.\n")
}

cat("\nScript Completed.\n")
