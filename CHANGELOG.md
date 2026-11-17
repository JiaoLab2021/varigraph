# Changelog

## v1.1-r2 (overlapped index merging, parameter errors and logs) (2026-11-16)

### Added

- Genotyping and GPU read collection report progress.

### Changed

- Merge graph-index results in graph order while workers keep claiming nodes (on one extra thread when `-t` > 1), instead of after all nodes; CPU and GPU share this code.
- Show the haplotype selection of each window only with `--debug`, fold the loaded-graph summary into its step, and drop redundant success lines.
- Print log numbers with thousands separators except in the final resource line, and use the same command descriptions in all help pages.
- Build the host code of `varigraph-gpu` with `-O3`, as the CPU build.

### Fixed

- Report an invalid numeric option value as a parameter error followed by the command help.

## v1.1-r1 (CLI, logging and scheduling) (2026-11-11)

### Added

- Bash completion for commands, files and directories.

### Changed

- Simplify CPU/GPU command handling, help and program metadata.
- Unify logs, debug control and progress reporting.
- Improve CPU task scheduling and resource cleanup.

### Fixed

- Fix build compatibility, argument validation and input/output error handling.
