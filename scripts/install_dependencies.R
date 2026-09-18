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

source("scripts/dependencies.R")

cat("ATOP-SCREEN dependency installer\n")
options(timeout = max(600, getOption("timeout", 60)))
repos <- getOption("repos")
if (is.null(repos) || !"CRAN" %in% names(repos) || is.na(repos[["CRAN"]]) || repos[["CRAN"]] == "@CRAN@") {
    repos <- c(CRAN = "https://cloud.r-project.org")
}
options(repos = repos)

install_missing <- function(packages, installer) {
    for (package in packages) {
        if (package_ready(package)) next
        cat("Installing", package, "\n")
        tryCatch(installer(package), error = function(e) {
            message("Installation failed for ", package, ": ", conditionMessage(e))
        })
    }
}

install_missing(dependency_groups$CRAN, function(package) {
    install.packages(package, dependencies = NA)
})

bioc_missing <- dependency_groups$Bioconductor[
    !vapply(dependency_groups$Bioconductor, package_ready, logical(1))]
github_missing <- names(dependency_groups$GitHub)[
    !vapply(names(dependency_groups$GitHub), package_ready, logical(1))]

if (length(bioc_missing) || length(github_missing)) {
    install_missing("BiocManager", function(package) install.packages(package))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        stop("BiocManager could not be installed. Check the installation log and rerun this script.", call. = FALSE)
    }
    options(repos = BiocManager::repositories())
    install_missing(bioc_missing, function(package) {
        BiocManager::install(package, ask = FALSE, update = FALSE, dependencies = NA)
    })
}

if (length(github_missing)) {
    install_missing("remotes", function(package) install.packages(package))
    if (!requireNamespace("remotes", quietly = TRUE)) {
        stop("remotes could not be installed. Check the installation log and rerun this script.", call. = FALSE)
    }
    install_missing(github_missing, function(package) {
        remotes::install_github(dependency_groups$GitHub[[package]],
            dependencies = NA, upgrade = "never", build_vignettes = FALSE)
    })
}

missing <- check_packages()
if (length(missing)) {
    stop(paste("Installation incomplete:", paste(missing, collapse = ", "),
               "\nReview the errors above and rerun in a fresh R session."), call. = FALSE)
}
cat("All required R packages are available.\n")
cat("Next: Rscript scripts/check_system.R\n")
cat("MAGeCK, compiler toolchains and Arial are installed separately; see README.md.\n")
