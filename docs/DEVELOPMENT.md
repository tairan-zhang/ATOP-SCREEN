# Build from source on an Apple silicon Mac

For normal use, [install the Mac app](INSTALL.md). To work on the source, use Git, Node.js 22.19.0 with npm 10, and the Xcode Command Line Tools (including Python 3, `otool`, `install_name_tool` and `codesign`).

## Desktop development and packaging

The matching DMG supplies the prebuilt R framework, analysis library and MAGeCK runtime. Install the app in Applications, then run the following in Terminal:

```sh
git clone --branch v0.1.3 https://github.com/tairan-zhang/ATOP-SCREEN.git
cd ATOP-SCREEN
runtime=packaging/runtime/macos-arm64
resources="/Applications/ATOP SCREEN.app/Contents/Resources/atop"
mkdir -p "$runtime/analysis"
ditto "$resources/runtime" "$runtime/runtime"
ditto "$resources/analysis/lib" "$runtime/analysis/lib"
ditto "$resources/mageck" "$runtime/mageck"
npm ci
npm run verify
```

Launch the development interface using the copied runtime:

```sh
export ATOP_RSCRIPT="$PWD/packaging/runtime/macos-arm64/runtime/R.framework/Resources/bin/Rscript"
npm run dev
```

Build a local DMG and check the app bundle:

```sh
CSC_IDENTITY_AUTO_DISCOVERY=false npm run package:mac -- --publish never
codesign --verify --deep --strict "release/0.1.3/mac-arm64/ATOP SCREEN.app"
```

The installer is written to `release/0.1.3/ATOP-SCREEN-0.1.3-macos-arm64.dmg`. The packaging hook makes R subprocesses use the bundled runtime and applies an ad hoc signature. Apple Developer ID signing and notarization require separate credentials.

## Locked R sources

This step is for compiling the analysis dependencies or building the analysis-worker container. The pinned GseaVis, renv and MAGeCK sources are under `packaging/sources/`. The remaining R package archives are provided as a release asset. From the repository root:

```sh
curl --fail --location --output /tmp/ATOP-SCREEN-0.1.3-r-sources.tar \
  https://github.com/tairan-zhang/ATOP-SCREEN/releases/download/v0.1.3/ATOP-SCREEN-0.1.3-r-sources.tar
tar -xf /tmp/ATOP-SCREEN-0.1.3-r-sources.tar
(cd packaging/sources && shasum -a 256 -c SHA256SUMS)
```

With R 4.6.1 and the required build toolchain, `Rscript --vanilla scripts/restore-r-runtime.R --library=/absolute/path/to/library` restores the pinned dependency set into a dedicated library. The [runtime manifest](../packaging/runtime-manifest.json), [R lockfile](../renv.lock) and [source checksums](../packaging/sources/SHA256SUMS) identify the versions and source archives.

## Analysis checks

`python3 scripts/check-worker-designs.py` exercises the paired, shared-baseline and group-mean designs using the bundled runtimes. It checks repeatability for the R and C++ permutation engines and runs MAGeCK RRA. Source-level R tests use `Rscript --vanilla scripts/run-source-r-tests.R` with the development package `testthat` installed.
