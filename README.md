# ATOP-SCREEN

ATOP-SCREEN analyzes CRISPR screen count matrices on your Mac. It ranks genes for enrichment and depletion using Adaptive Top-N scoring and directional permutation tests. From a gene-level result, you can inspect individual sgRNA counts, compare MAGeCK RRA results, or explore pathways with GSEA.

**[Download ATOP-SCREEN for Mac (Apple silicon, DMG)](https://github.com/tairan-zhang/ATOP-SCREEN/releases/download/v0.1.3/ATOP-SCREEN-0.1.3-macos-arm64.dmg)** · [Installation and first analysis](docs/INSTALL.md) · [中文安装指南](docs/INSTALL.zh-CN.md)

Open the disk image, drag **ATOP SCREEN** to **Applications**, then open it from Applications. R, the analysis packages, MAGeCK and its Python runtime are included. The app runs locally and accepts your own count matrices and gene sets. The disk image is for Apple silicon Macs. For the first-launch macOS security prompt, see the [installation guide](docs/INSTALL.md#first-launch).

## Analyze a screen

1. Open a project and choose **Analyze counts**. Import a table with one row per sgRNA, a guide ID, a target gene and one raw-count column per sample.
2. Select the endpoint and baseline columns. Use **Matched pairs** for known pairs, **Shared baseline** for several endpoints measured against a single starting sample, or **Ratio of group means** for independent groups, including unequal group sizes.
3. Choose **Start analysis**. ATOP ranks the guide effects within each gene and averages the upper and lower two-thirds separately, rounding the retained count up. A guide-label permutation test and Benjamini–Hochberg adjustment are applied separately to each direction. An optional concordance filter (NCS) acts before ranking when replicate contrasts are available.
4. Inspect the guide and gene tables or render a count, ranking or volcano plot. For pathway analysis, choose **Enrich pathways**, select a GMT file with identifiers matching the count table, and run GSEA. Figures export as SVG, PNG and PDF.

Try the [small example count matrix](packaging/container/example-counts.tsv) with `E1`, `E2` and `E3` as endpoints and `T0` as the shared baseline. The [A375 vemurafenib screen](https://raw.githubusercontent.com/tairan-zhang/ATOP-SCREEN/v0.1.3/test_data/A375_rawcount.txt) provides a paired example using `PLX7_R1`/`D7_R1` and `PLX7_R2`/`D7_R2`. [Shalem et al., Science (2014)](https://doi.org/10.1126/science.1247005)

## Reproducing an analysis

Desktop **0.1.3** includes the **0.1.2** analysis engine, R **4.6.1**, Bioconductor **3.23** and MAGeCK **0.5.9.5**. Project settings record the sample mapping, guide filters, Top-N rule, permutation count, seed, engine and input checksum. Keep these settings with the count matrix to reproduce a run. See the [R lockfile](renv.lock), [runtime manifest](packaging/runtime-manifest.json) and [developer build instructions](docs/DEVELOPMENT.md).

## Source and licensing

The [R/C++ analysis engine](packages/analysis-engine), [Electron desktop process](electron) and [React interface](web/src) are included here. The engine declares GPL-3.0-or-later, and the desktop code carries its MIT notice. Third-party components retain their own notices. [GPL license](LICENSE) · [Desktop license](LICENSE-desktop)
