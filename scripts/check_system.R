# SPDX-License-Identifier: GPL-3.0-or-later
# check_system.R
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

cat("CRISPR Screen Analysis System - System Check\n")
cat("===============================\n\n")

cat("Basic Environment Check:\n")
cat("R Version:", R.Version()$version.string, "\n")
cat("Platform:", R.Version()$platform, "\n")
cat("OS:", Sys.info()["sysname"], "\n\n")

cat("Critical Package Check:\n")
critical_packages <- c(
    "shiny", "shinyjs", "readxl", "dplyr", "tidyr", "stringr", "DT",
    "ggplot2", "ggpubr", "ggnewscale", "patchwork", "RColorBrewer",
    "ggrepel", "ggpp", "enrichplot", "clusterProfiler",
    "data.table", "parallel", "zip", "later", "processx", "Rcpp"
)

all_ok <- TRUE
for (pkg in critical_packages) {
    status <- if (requireNamespace(pkg, quietly = TRUE)) "[OK]" else "[MISSING]"
    cat(sprintf("  %-15s %s\n", pkg, status))
    if (status == "[MISSING]") all_ok <- FALSE
}

cat("\nSystem Resources:\n")
cat("CPU Cores:", parallel::detectCores(), "\n")

cat("\nProject File Check:\n")
required_files <- c(
    "app.R", "src/cpp_permutation_engine.cpp",
    "R/analysis/crispr_analysis_functions.R", "R/plotting/gsea_functions.R"
)

files_ok <- TRUE
for (file in required_files) {
    status <- if (file.exists(file)) "[OK]" else "[MISSING]"
    cat(sprintf("  %-30s %s\n", file, status))
    if (status == "[MISSING]") files_ok <- FALSE
}

cat("\nCheck Summary:\n")
if (all_ok && files_ok) {
    cat("System ready, you can run source('app.R')\n")
} else {
    cat("System NOT ready:\n")
    if (!all_ok) cat("  - Please run source('scripts/install_dependencies.R')\n")
    if (!files_ok) cat("  - Please check project file integrity\n")
}

cat("\nFor detailed installation guide, please check README.md\n")
