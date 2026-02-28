# ATOP-SCREEN

An Integrated Graphical User Interface Platform for Pooled CRISPR Screen Analysis

ATOP-SCREEN is an integrated R/Shiny-based graphical user interface (GUI) platform for pooled CRISPR screen analysis. The platform implements an Adaptive Top-N gene-level aggregation strategy to reduce dilution from low-effect guides under heterogeneous sgRNA signals, and supports permutation-based empirical significance testing and interactive visualization within a unified analytical workflow.

## Key Features

### Core Analysis Algorithm (ATOP)

*   **Adaptive Top-N Aggregation Algorithm**: Adaptively aggregates a top-ranked subset of sgRNAs per gene to reduce dilution from low-effect guides while retaining a majority of guides.
    *   **Dynamic Selection**: Uses Top-$N_g$ aggregation where $N_g=\lceil 2n_g/3\rceil$.
    *   **Configurable Filtering**: User-selectable sgRNA count threshold (default: 3). Genes with fewer sgRNAs than the threshold are excluded from analysis to ensure statistical reliability. In the manuscript benchmark analyses, the default $n_{min}=3$ setting was used.
    *   **Transparent Exclusion**: All filtered genes and their sgRNAs are displayed in a separate table for review and download.
    *   **Dual Scoring**: Reports separate GenePositiveScore ($S_g^{+}$) and GeneNegativeScore ($S_g^{-}$) as the mean of the top-$N_g$ highest or lowest sgRNA LFCs, respectively.
*   **Hybrid C++/R Permutation Backend**: Implements a high-performance C++ permutation engine (via Rcpp) when available, with a parallel R-based fallback. Supports multi-threading and can utilize multiple CPU cores (default: 10,000 permutations for empirical p-value estimation).

### Alternative Analysis Methods

*   **MAGeCK Integration**: Seamlessly integrated MAGeCK RRA and MLE algorithms with:
    *   Automatic MAGeCK executable detection or user-specified path configuration
    *   Progress monitoring with status updates
    *   Native MAGeCK output display with original column names
    *   One-click download of raw MAGeCK results

### Complete Analytical Workflow

*   **Data Processing**:
    *   Input validation and preprocessing (supports common tabular formats e.g., .csv, .txt, .tsv; .xlsx supported when appropriate readers are available)
    *   Median-of-ratios normalization. A fixed pseudocount of 1 is applied prior to median-of-ratios normalization to define geometric means, while a separate user-configurable stabilization pseudocount $\alpha$ (default: 1) is applied only in the log-ratio step.
    *   sgRNA-level log2 fold-change (LFC) calculation with replicate tracking
    *   Gene-level scoring and permutation-based statistical testing
*   **Streamlined Results Display**:
    *   **sgRNA Level Data**: All guide-level metrics including per-replicate LFCs
    *   **Gene Level Summary Data**: Scores, rankings, P-values, and sgRNA counts
    *   **Filtered Genes Data**: Excluded genes with insufficient sgRNA coverage
*   **Integrated GSEA**: Built-in Gene Set Enrichment Analysis (GSEA) based on ranked gene-level scores using the `clusterProfiler` framework.
*   **Integrated Visualization Modules**:
    *   Volcano plots and gene-ranking plots
    *   Paired sgRNA plots for guide-level consistency inspection
    *   GSEA enrichment and lollipop plots

## Project Structure


The codebase is organized to separate UI, server logic, and computational backends:

```text
ATOP/
├── app.R                       # Main entry point (UI definition & Server initialization)
├── R/                          # R Modular Logic
│   ├── analysis/               # Core statistical analysis & MAGeCK wrappers
│   ├── interface/              # R-C++ interfaces (Rcpp)
│   ├── plotting/               # Visualization functions (Volcano, GSEA, etc.)
│   └── server/                 # Shiny server modules and helper functions
├── src/                        # High-Performance Backend
│   └── cpp_permutation_engine.cpp # C++ implementation of permutation testing
├── scripts/                    # Utility Scripts
│   ├── install_dependencies.R  # Dependency installer
│   ├── check_system.R          # System integrity verifier
│   └── verify_cpp_compilation.R # C++ compilation verification
├── config/                     # Configuration Files
│   └── mageck_path.txt         # MAGeCK executable path (auto-generated)
├── test_data/                  # Sample Datasets
│   ├── A375_rawcount.txt       # Example CRISPR screen data
│   └── c2.all.v2024.1.Hs.symbols.gmt # MSigDB gene sets for GSEA
├── requirements.txt            # R package dependencies
└── LICENSE                     # GPL-3.0 License
```


## Installation

### Prerequisites

*   R (>= 4.0.0)
*   A C++ compiler (e.g., RTools on Windows, Xcode Command Line Tools on macOS, or GCC on Linux) for compiling the optional C++ backend
*   **(Optional) MAGeCK**: If you plan to use the MAGeCK RRA or MLE modules within the application, you must have [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) installed on your system.

### Setup

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/tairan-zhang/ATOP-SCREEN.git
    cd ATOP-SCREEN
    ```

2.  **Install Dependencies:**
    Run the provided helper script to install all required R packages:
    ```r
    Rscript scripts/install_dependencies.R
    ```
    *Alternatively, package dependencies can be reviewed in `requirements.txt` for reproducibility.*

3.  **Configure MAGeCK Path (Optional):**
    If you intend to use the MAGeCK analysis methods, you must configure the path to the executable. Open the file `config/mageck_path.txt` and replace its contents with the absolute path to your `mageck` executable. Ensure there are no extra spaces or empty lines.
    
    Example (`config/mageck_path.txt`):
    ```
    /opt/miniconda3/envs/mageck-vispr/bin/mageck
    ```

## Usage

1.  **Launch the Application:**
    Open the project in RStudio or run from the terminal:
    ```r
    R -e "shiny::runApp()"
    ```

2.  **Workflow Steps:**
    *   **Upload**: Load your raw CRISPR screen count matrix (.csv, .txt, .tsv, .xlsx).
    *   **Define Columns**: Select columns for sgRNA ID, Gene Symbol, and replicate conditions (Treatment vs. Control).
    *   **Configure Parameters**: 
        *   Choose analysis method (ATOP, MAGeCK RRA, or MAGeCK MLE)
        *   Set sgRNA filtering threshold (genes with fewer sgRNAs will be excluded)
        *   Configure permutation parameters (default: 10,000 permutations)
    *   **Process**: Run the analysis pipeline and monitor progress.
    *   **Review Results**: Examine three consolidated tables:
        *   sgRNA-level data with per-replicate fold changes
        *   Gene-level summary with statistical significance
        *   Filtered genes (if any were excluded)
    *   **Analyze**: Perform GSEA or visualize ranked candidate genes.
    *   **Export**: Download result tables and plots.

## Test Data

The `test_data/` directory contains sample datasets to demonstrate the analysis pipeline:

1.  **`A375_rawcount.txt`**
    *   **Description**: Public genome-scale CRISPR-Cas9 screen in A375 melanoma cells treated with vemurafenib. The included example corresponds to PLX7 (vemurafenib, 7 days) versus D7 (vehicle control), as described in the manuscript.
    *   **Source**: Shalem, O., Sanjana, N.E., Hartenian, E., Shi, X., Scott, D.A., Mikkelson, T., Heckl, D., Ebert, B.L., Root, D.E., Doench, J.G., Zhang, F. (2014). Genome-scale CRISPR-Cas9 knockout screening in human cells. *Science*, 343, 84-87. [https://doi.org/10.1126/science.1247005](https://doi.org/10.1126/science.1247005).

2.  **`c2.all.v2024.1.Hs.symbols.gmt`**
    *   **Description**: Curated gene sets (C2 collection) from MSigDB used for Gene Set Enrichment Analysis (GSEA).

## License

This project is licensed under the **GPL-3.0 License**. See the `LICENSE` file for details.
