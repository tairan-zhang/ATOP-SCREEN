# ATOP SCREEN security

ATOP SCREEN is local-first. The main security boundaries are Electron IPC,
filesystem access, the R worker, and packaged runtime assets.

- Keep `contextIsolation`, renderer sandboxing, and disabled Node integration
  enabled.
- Add new desktop operations to the typed preload contract and validate their
  payloads in `electron/ipc`.
- Keep asset reads inside registered project/artifact roots and enforce file
  type and size limits.
- Normalize project names and reject unsafe workspace paths.
- Never commit credentials, tokens, personal research data, or generated
  project directories.
- Review Electron, npm, R, and packaged runtime dependencies before release.

For a suspected vulnerability, do not publish sensitive details in a public
issue. Contact the project maintainer with a reproducible description,
affected version, and a minimal safe proof of impact.
