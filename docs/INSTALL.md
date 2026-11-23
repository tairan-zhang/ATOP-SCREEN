# Install ATOP-SCREEN on macOS

## Install the app

1. In **Apple menu → About This Mac**, check the **Chip** field. The installer is for Apple silicon (M-series).
2. [Download the macOS disk image](https://github.com/tairan-zhang/ATOP-SCREEN/releases/download/v0.1.3/ATOP-SCREEN-0.1.3-macos-arm64.dmg).
3. Double-click the DMG and drag **ATOP SCREEN** onto **Applications** in the same window.
4. Eject the disk image. In Finder, open **Applications → ATOP SCREEN**.

The app includes the R runtime, C++ permutation engine, pathway-analysis packages and MAGeCK RRA with Python. You can run the analysis offline after downloading the app and your input files.

### First launch

This build uses an ad hoc signature and has not been notarized by Apple. If macOS says it cannot verify the developer, dismiss the message after the first open attempt. For the app downloaded from this repository, go to **System Settings → Privacy & Security → Open Anyway**, authenticate if prompted and select **Open**. The exception applies to this app. [Apple’s first-launch instructions](https://support.apple.com/en-us/102445)

## Run the example

1. [Save the example count table](https://raw.githubusercontent.com/tairan-zhang/ATOP-SCREEN/v0.1.3/packaging/container/example-counts.tsv) as `example-counts.tsv`.
2. Create a project. Choose **Analyze counts** and select the saved table. Map the sgRNA ID to `guide` and the target gene to `gene`.
3. Choose **Shared baseline** under **LFC comparison**. Select `E1`, `E2`, `E3` under **Endpoint / treatment counts** and `T0` under **Baseline / control counts**.
4. Select **Run MAGeCK RRA** if you want its results alongside ATOP, then choose **Start analysis**.
5. Open **Screen analysis → Results** to see the guide and gene tables. Under **Visualization**, select a gene plot and render the figure. You can save it as SVG, PNG or PDF.

For your own screen, import a CSV, TSV, TXT or XLSX count table. Use unique sgRNA identifiers and nonnegative integer read counts. To analyze pathways, supply a GMT file whose gene identifiers match the target-gene column.
