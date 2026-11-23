# Contributing to ATOP SCREEN

Keep contributions inside the ATOP SCREEN product boundary: local CRISPR
screen analysis, project persistence, pathway enrichment, and visualization.

Before opening a change:

```bash
npm install
npm run typecheck
npm test
npm run build
```

Use `npm run verify` before release work. Runtime dependency changes require an
update to `packaging/runtime-manifest.json` and a clean run of the bundled R
worker; do not remove `clusterProfiler` or `GseaVis` without a separately
reviewed GSEA migration.

Use the existing `src/app`, `src/features`, `src/components`, `src/domain`,
`src/services`, and `src/types` boundaries. Desktop operations belong behind
the preload bridge; R algorithms belong in `packages/analysis-engine`.

Describe user-visible changes, tests run, and any documentation or contract
changes in the pull request. Do not commit generated project workspaces,
private datasets, credentials, or machine-specific paths.
