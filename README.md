# ATOP-SCREEN

**Version 0.0.1** · An R/Shiny application for pooled screen analysis.

ATOP-SCREEN is an integrated graphical user interface (GUI) platform for pooled screen analysis based on the Adaptive Top-N aggregation algorithm. It integrates data processing, permutation-based significance testing, and interactive visualization within a unified workflow, linking gene-level scores to guide-level diagnostic evidence and pathway-level interpretation through Gene Set Enrichment Analysis (GSEA).

## Quick start

Install a current release of [R](https://cran.r-project.org/) and Git. The application has been tested locally with R 4.6. Run the following commands in a terminal from the repository root.

```sh
git clone https://github.com/tairan-zhang/ATOP-SCREEN.git
cd ATOP-SCREEN
Rscript scripts/install_dependencies.R
Rscript scripts/check_system.R
Rscript -e "shiny::runApp('.', launch.browser = TRUE)"
```

The installer uses CRAN, the [Bioconductor release matched to your R version](https://bioconductor.org/install/), and the [official GseaVis GitHub repository](https://github.com/junjunlab/GseaVis) (version 0.1.1 or later). It installs required dependencies, including GseaVis's declared GitHub dependencies, and fails if any required package is missing or cannot load. Rerun it in a fresh R session after resolving an installation error. `requirements.txt` is a reference list, not a pip requirements file or a version lockfile.

- **Build tools:** Source packages may need Rtools on Windows, Xcode Command Line Tools on macOS, or development compilers/libraries on Linux. The application's C++ engine can be checked with `Rscript scripts/verify_cpp_compilation.R`; analysis has an R fallback.
- **MAGeCK (optional):** Install [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) separately and set its executable path in the application or `config/mageck_path.txt`.
- **Fonts and PDF:** Install Arial for consistent typography. PDF export uses Quartz on macOS and requires Cairo support in R on other systems.

These system components are not installed by the R dependency script.

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
