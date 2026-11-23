# Desktop distribution

The [macOS installer](https://github.com/tairan-zhang/ATOP-SCREEN/releases/download/v0.1.3/ATOP-SCREEN-0.1.3-macos-arm64.dmg) contains the Electron app, ATOP R/C++ engine, R 4.6.1 and the locked analysis packages, and MAGeCK RRA with Python. The [installation guide](INSTALL.md) covers drag-and-drop installation and first launch. The [build guide](DEVELOPMENT.md) describes how to assemble the same runtime from the release and create the disk image.

For a distribution check, copy the app out of the DMG to a separate location and eject the image. Verify the copied app with `codesign --verify --deep --strict`. Launch it with an empty project and run the [example count table](../packaging/container/example-counts.tsv) against its shared baseline. Check the gene and guide tables, MAGeCK results, and SVG, PNG and PDF exports. Use a GMT file with matching gene identifiers to verify GSEA and a running enrichment plot.
