# ATOP-SCREEN

**Version 0.0.1** · An R/Shiny application for pooled screen analysis.

ATOP-SCREEN is an integrated graphical user interface (GUI) platform for pooled screen analysis based on the Adaptive Top-N aggregation algorithm. It integrates data processing, permutation-based significance testing, and interactive visualization within a unified workflow, linking gene-level scores to guide-level diagnostic evidence and pathway-level interpretation through Gene Set Enrichment Analysis (GSEA).

## Quick start

Requirements: R 4.0 or later and compatible R/Bioconductor packages. A C++ compiler enables the Rcpp permutation engine; an R fallback is available.

```sh
git clone https://github.com/tairan-zhang/ATOP-SCREEN.git
cd ATOP-SCREEN
Rscript scripts/install_dependencies.R
Rscript -e "shiny::runApp('.', launch.browser = TRUE)"
```

For MAGeCK, install it separately and configure its executable path in the application or `config/mageck_path.txt`. GSEA enrichment curves require `GseaVis`. Install Arial on your system for consistent plot typography.

## Workflow

1. **Analysis:** Upload a count table (CSV, TSV, TXT or Excel), select gene/sgRNA identifiers and treatment/control columns, then run Screen analysis. Upload a GMT gene-set file to run GSEA.
2. **Results:** Review and download sgRNA, gene and pathway tables.
3. **Visualization:** Create paired sgRNA, volcano, gene-ranking, pathway lollipop and enrichment plots. Set dimensions in centimetres, select **Apply** to update the preview, and download the figure.

## Example data

- `test_data/A375_rawcount.txt`: A375 vemurafenib screen, PLX7 versus D7, from [Shalem et al., Science (2014)](https://doi.org/10.1126/science.1247005).
- `test_data/c2.all.v2024.1.Hs.symbols.gmt`: MSigDB C2 gene sets for pathway enrichment.

## Development

`app.R` defines the application; `R/` contains analysis, plotting and server modules; `src/` contains the C++ permutation engine; `www/` contains interface assets.

Run checks from the repository root:

```sh
Rscript scripts/test_permutation.R
Rscript scripts/test_modules.R
Rscript scripts/test_app.R
Rscript scripts/test_figure_style.R
```

## License

[GPL-3.0](LICENSE).
