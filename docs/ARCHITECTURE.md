# ATOP SCREEN architecture

## 1. Scope and goals

ATOP SCREEN is a local-first Electron application. Its stable product flow is:

1. inspect a screen matrix and choose columns from the detected header;
2. run the original ATOP analysis through the R worker and optional C++
   permutation library;
3. persist the project workspace and analysis snapshots locally;
4. run optional GSEA and ATOP plotting functions after analysis succeeds;
5. preview and download generated SVG, PNG, and PDF assets.

The architecture optimizes for clear boundaries, one-way dependencies, testable
pure logic, and a small desktop attack surface.

## 2. System context

```text
User
  │
  ▼
React/Vite renderer ── typed preload bridge ── Electron main process
       │                                      │
       │                                      ├─ project/file repositories
       │                                      └─ R JSON-lines worker ── ATOP R/C++ core
       │
       └─ preview/download requests ─────────────── local project files
```

No network service is required for analysis. Input data, project manifests,
intermediate RDS files, result JSON, and generated figures remain on the local
machine unless the user explicitly downloads an asset elsewhere.

## 3. Component responsibilities

| Boundary | Responsibility | Must not own |
| --- | --- | --- |
| `web/src/app` | Compose workbench state and route views | R algorithms or filesystem calls |
| `web/src/features` | Business workflows and feature-level validation | Electron IPC implementation |
| `web/src/components` | Reusable ATOP presentation and interaction primitives | Project orchestration |
| `web/src/domain` | Pure formatting, defaults, and browser storage helpers | UI rendering or Electron |
| `web/src/services` | Renderer-side API facade and error boundary | Direct `ipcRenderer` access |
| `electron/preload` | Expose the minimum typed desktop surface | Business rules |
| `electron/ipc` | Validate requests and dispatch use cases | UI state |
| `electron/main` | Worker lifecycle, file/project repositories, window | React concerns |
| `packages/analysis-engine` | ATOP R/C++ algorithms and plot functions | Electron or JSON transport |
| `packages/ipc-contracts` | Shared request/result shapes | Implementation logic |

## 4. Repository structure

```text
ATOP-SCREEN/
├── web/src/
│   ├── app/                 # application composition
│   ├── features/            # business modules
│   ├── components/          # shared ATOP UI
│   ├── domain/              # pure domain helpers
│   ├── services/            # renderer API facade
│   ├── types/               # renderer contracts
│   └── styles/              # single ATOP style entry
├── electron/
│   ├── main/                # desktop runtime and repositories
│   ├── ipc/                 # protocol handlers and validation
│   └── preload/             # isolated bridge
├── packages/
│   ├── analysis-engine/     # R/C++ ATOP core
│   └── ipc-contracts/       # shared IPC types
├── fixtures/                # small deterministic test inputs
├── packaging/               # app icon and bundled runtime assets
├── docs/                    # maintained product documentation
└── package-lock.json        # reproducible JavaScript dependency graph
```

### Analysis-engine module layout

```text
packages/analysis-engine/R/
├── atop/                    # ATOP workflow only
│   ├── validation.R          # request and column contracts
│   ├── preprocessing.R       # count preparation, normalization and LFCs
│   ├── scoring.R             # adaptive Top-N and BH correction
│   ├── statistics.R          # permutation dispatch
│   └── orchestrator.R        # stable public analysis entry point
├── permutation/              # implementation detail of significance testing
│   ├── r_parallel.R          # compatibility fallback
│   ├── engine_selector.R     # C++/R selection and fallback policy
│   └── cpp/                  # sole compiled-engine adapter
├── gsea/                     # preserved clusterProfiler workflow
├── plotting/
│   ├── gene/                # volcano and ranking builders
│   ├── gsea/                # lollipop and GseaVis running plots
│   └── sgrna/               # paired-count builder
└── runtime/                  # worker-neutral progress support
```

The worker uses an explicit module registry rather than recursively sourcing
files. A new analysis module must declare its dependencies through that
registry, expose a narrow function surface, and must not import Electron or
write JSON protocol messages directly.

## 5. Dependency rules

Dependencies point inward toward stable contracts and pure logic:

- `app → features → components/domain/services/types`.
- `components` may use `domain` formatting, but never import Electron.
- `services` imports only `packages/ipc-contracts` and `window.atop` types.
- `electron/ipc` validates at the boundary before calling `electron/main` use
  cases or repositories.
- `electron/main` may launch the R worker, but R code never imports Electron.
- The R package depends on R/C++ packages only and prints JSON protocol output
  through its worker adapter.
- New cross-feature types belong in `packages/ipc-contracts` or
  `web/src/types/atop.ts`; do not create feature-to-feature imports.
- New UI controls must use the shared component and style tokens instead of
  introducing a second design system.

## 6. Data flow

`AnalysisView` asks `services/desktop-api.ts` to inspect a selected path. The
main process validates the command, queues it for the worker, and returns a
typed result. The view derives column options from `FileInspection.columns`,
then submits only selected column names and ATOP parameters.

The worker sources the analysis modules on demand. Its stdout is reserved for
JSON responses; progress and diagnostics are sent through stderr and forwarded
to the renderer status line. Successful analysis results are passed to the
project store, which copies source/assets into the project workspace and writes
JSON snapshots atomically.

R domain functions expose an optional `progress_callback` and do not own worker
or renderer state. The worker supplies a callback that converts status messages
to protocol-safe stderr events. The same callback contract is used by the ATOP
orchestrator, GSEA analysis, R-parallel fallback and the sole C++ adapter at
`R/permutation/cpp/interface.R`.

GSEA and visualization follow the same path and are available only after a
valid analysis result exists. The artifact panel reads only paths granted by
the main-process asset policy.

## 7. IPC contract

The canonical command union is defined in
[`packages/ipc-contracts/src/index.ts`](../packages/ipc-contracts/src/index.ts):

- `health`
- `inspect_file`
- `process_data`
- `gsea`
- `visualize`
- `restore_state`

Project channels are `projects:list`, `projects:create`, `projects:open`,
`projects:save-analysis`, `projects:save-gsea`,
`projects:save-visualization`, and `projects:delete`. File channels are limited
to choosing inputs, reading approved assets, and downloading an approved asset.

The preload bridge exposes these operations through `window.atop`; renderer
features must not call `ipcRenderer` directly.

## 8. Data model

Each project has a UUID and uses the folder selected in the desktop interface. A completed screen analysis writes its input and results below the project manifest:

```text
project-folder/
├── manifest.json
├── analyses/<analysis>/input/
├── analyses/<analysis>/screen/      # settings, result tables, RDS files and runs/
└── plots/
```

`manifest.json` is the index. `registry.json` under Electron `userData`
contains only project ids and manifest paths. Writes use a temporary file and a
rename so an interrupted write does not replace a valid manifest.

## 9. Error model

- Invalid IPC payloads fail at validation with a descriptive `TypeError`.
- Missing or unsupported files fail before entering R.
- Worker failures reject the pending request and reset the worker state.
- R errors are returned as a user-facing message and mirrored in the status
  line; raw stdout is never mixed into the JSON protocol.
- Asset access rejects paths outside registered project/artifact roots,
  unsupported extensions, non-files, and oversized files.
- A completed calculation is not silently discarded if project persistence
  fails; the UI reports that distinction explicitly.

## 10. Reliability and performance

- The analysis process owns one worker and a bounded in-memory request queue.
- Cancellation rejects the request and terminates the active worker so the next
  request starts from a clean state.
- Project snapshots are written atomically and source files are copied into the
  project before the manifest is updated.
- Resize handlers update the DOM on `requestAnimationFrame` and commit React
  state once at the end of the gesture.
- Large result tables use page-size controls and render only the selected page.

## 11. Observability

The renderer exposes runtime health, loaded R modules, function availability,
and the latest progress message. The main process forwards worker diagnostics
as progress events. Project manifests record timestamps, status, input name,
analysis count, and generated plot metadata. There is no network telemetry.

## 12. Security

- Electron uses `contextIsolation`, disabled Node integration, and a sandboxed
  renderer.
- The preload bridge is the only renderer-to-main boundary.
- Browser window creation is denied.
- File reads and downloads use an allowlist of project and artifact roots.
- User-controlled names are normalized before they become workspace names.
- The CSP blocks external scripts and object embedding.
- Do not commit personal paths, credentials, generated project workspaces, or
  runtime secrets.

## 13. Build, deployment, and test

`npm run verify` runs type-checking, tests, and the renderer build. `npm run package:dir`
creates an unpackaged Electron application under `release/`; macOS installer
creation is available through `npm run package:mac` when signing/runtime
requirements are present.

The macOS R runtime is governed by `packaging/runtime-manifest.json`; its
package closure deliberately retains `clusterProfiler` and `GseaVis` because
they are part of the supported GSEA result and plotting pipeline. The runtime
prune script removes only packages outside that declared closure plus the exact
unused paths listed in the manifest.

Quality gates are `npm run typecheck`, `npm test`, and `npm run build`. The test
suite covers IPC validation and repositories in Node, plus pure/SSR-safe
renderer components and column mapping in Vitest.

## 14. ADR index

- [ADR-0001: Electron and local ATOP boundary](adr/0001-electron-atop-boundary.md)

Changes that alter IPC, project persistence, worker lifecycle, or dependency
direction require a new ADR.
