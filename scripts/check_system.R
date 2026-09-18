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

source("scripts/dependencies.R")

cat("ATOP-SCREEN system check\n")
cat("R:", R.version.string, "\nPlatform:", R.version$platform, "\n\n")
missing_packages <- check_packages()
required_files <- c("app.R", "VERSION", "www/app.css", "www/app.js", "www/atop-logo.png",
                    "src/cpp_permutation_engine.cpp", "R/analysis/crispr_analysis_functions.R",
                    "R/plotting/gsea_functions.R")
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_packages) || length(missing_files)) {
    stop(paste("System check failed.",
        if (length(missing_packages)) paste("Install packages:", paste(missing_packages, collapse = ", ")),
        if (length(missing_files)) paste("Missing files:", paste(missing_files, collapse = ", "))), call. = FALSE)
}
cat("\nRequired R packages and application files are available.\n")
cat("Check optional C++ compilation: Rscript scripts/verify_cpp_compilation.R\n")
cat("Launch from the repository root with the Rscript command in README.md.\n")
