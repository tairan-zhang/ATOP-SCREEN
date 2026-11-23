# ADR-0001: Electron and local ATOP boundary

- Status: accepted
- Date: 2026-08-28

## Context

ATOP SCREEN needs a native file picker, a long-running R process, local project
persistence, and preview/download access. Mixing those operations into React
made the product boundary difficult to test and made filesystem access too
easy to widen accidentally.

## Decision

Use Electron as the desktop shell. Keep the renderer behind a narrow preload
bridge, validate every IPC payload in the main process, and run the original
ATOP R/C++ implementation in a JSON-lines worker. Keep project persistence and
asset allowlisting in the main process.

## Consequences

- React remains focused on UI and feature state.
- R remains reusable outside Electron and can be tested independently.
- IPC contracts become explicit and reviewable.
- The packaged runtime is larger because it includes R, but analysis is local
  and reproducible.
- Any future protocol or persistence change must update the contract and an ADR.
