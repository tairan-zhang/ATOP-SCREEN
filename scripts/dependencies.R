dependency_groups <- list(
    CRAN = c("shiny", "shinyjs", "readxl", "dplyr", "tibble", "rlang", "tidyr",
             "stringr", "DT", "ggplot2", "svglite", "ggpubr", "ggnewscale",
             "patchwork", "RColorBrewer", "ggrepel", "ggpp", "data.table",
             "Rcpp", "zip", "processx", "later", "httpuv"),
    Bioconductor = c("clusterProfiler", "enrichplot"),
    GitHub = c(GseaVis = "junjunlab/GseaVis")
)

required_packages <- unique(c(dependency_groups$CRAN, dependency_groups$Bioconductor,
                             names(dependency_groups$GitHub)))

package_ready <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(FALSE)
    if (package == "GseaVis") {
        return(utils::packageVersion(package) >= "0.1.1" &&
               all(c("rankCol", "geneCol", "pvalSize") %in% names(formals(GseaVis::gseaNb))))
    }
    TRUE
}

check_packages <- function() {
    ready <- vapply(required_packages, package_ready, logical(1))
    for (package in required_packages) {
        cat(sprintf("[%s] %s\n", if (ready[[package]]) "OK" else "MISSING/INCOMPATIBLE", package))
    }
    names(ready)[!ready]
}
