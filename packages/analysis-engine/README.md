# ATOP-SCREEN analysis engine

This is the canonical internal R package for ATOP-SCREEN. The Electron worker
under `inst/worker/` translates typed JSON-lines requests into calls to the
domain functions here.

## R modules

- `R/atop/orchestrator.R`: the stable ATOP workflow entry point; it coordinates normalization, scoring, and permutation testing without owning desktop transport.
- `R/permutation/`: the R-parallel fallback, engine selection, and isolated C++ bridge (`cpp/`).
- `R/gsea/clusterprofiler.R`: the preserved `clusterProfiler` GSEA workflow.
- `R/plotting/gene/`: shared gene-figure primitives plus independent volcano and ranking builders.
- `R/plotting/gsea/`: separate lollipop and `GseaVis` running-plot builders.
- `R/plotting/sgrna/paired.R`: sgRNA paired-count figure builder.
- `R/runtime/progress_tracker.R`: UI-neutral progress bookkeeping retained for long-running calculations.
- `src/cpp_permutation_engine.cpp`: optional C++17 permutation implementation.
- `tests/testthat/`: source-parse and regression tests for the engine.
- `inst/extdata/`: package-owned data resources when a fixture must ship with the engine.

## React adapter

`inst/worker/analysis_worker.R` is a JSON-lines process. Each request has an
`id` and a `command`; each response is one JSON object with `ok`, `result`, or
`error`. It supports:

- `health`: R version, loaded modules, functions, and package availability.
- `inspect_file`: safe column and first-30-row preview for CSV/TXT/TSV/XLSX.
- `process_data`: calls `perform_crispr_screen_analysis()` and keeps full results in the worker session while returning compact tables to React.
- `gsea`: calls `perform_gsea_analysis_clusterprofiler()` using the stored analysis result.

The public GSEA analysis and lollipop entry points are
`perform_gsea_analysis_clusterprofiler()`, `generate_gsea_lollipop_plot()`,
and `generate_gsea_enrichment_plot_atop()`.
Long-running domain functions accept an optional `progress_callback`; the
worker owns conversion of callback messages to renderer progress events. Native
permutation dispatch uses only `R/permutation/cpp/interface.R` as its R/C++
adapter.

The worker uses `ATOP_ANALYSIS_ROOT`, `ATOP_RUNTIME_ROOT`, and
`ATOP_CPP_ENGINE_LIBRARY` supplied by Electron. Its stdout is reserved for JSON;
captured R progress is forwarded through stderr and shown in the renderer
status line. `packaging/runtime-manifest.json` is the authoritative packaged
R-runtime dependency declaration; it explicitly retains clusterProfiler and
GseaVis for the supported GSEA workflow.

## Module rules

Each file owns one domain concern and exposes a small function surface. The
worker's `module_paths` registry is the only place that defines source order;
add a new module there rather than relying on `list.files()` ordering. Public
function names remain stable so existing projects, artifacts, and Electron IPC
calls are unaffected by internal reorganization.
