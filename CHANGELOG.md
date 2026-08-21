# Changelog

All notable changes to this project are documented in this file.

## [1.1.0] - 2026-08-21

### Added

- Exploratory cell-level metrics for reversal potential, local slope
  conductance, apparent input resistance, rectification, and current retention.
- Exploratory apparent chord-conductance curves with configurable exclusion near
  each cell's interpolated reversal potential.
- Supplementary CSV, Excel, figure, JSON, and Markdown report outputs.
- Tests for nonstandard voltage grids, variable sampling and sweep dimensions,
  and raw-versus-Holm supplementary p-value reporting.

### Changed

- Adapted voltage targets to the observed voltage set instead of assuming a
  fixed sweep count or a -100 to +100 mV protocol.
- Expanded statistical reports with group means, SEM, 95% confidence intervals,
  Hedges' g, exact p values, and an explicit significance conclusion.
- Applied readable formatting to generated Excel workbooks, including frozen
  headers, filters, styled headings, and bounded column widths.
- Updated the skill instructions, example configuration, output schema,
  methodology notes, README, and Codex interface metadata.
- Added repository line-ending rules for stable cross-platform diffs.

### Validation

- Skill structure validation passed.
- All 10 unit tests passed.
- A 20-cell real-data regression produced 42 non-empty outputs while preserving
  all source ABF hashes.
