# Output schema

```text
analysis_workflow/
|-- 01_extracted/
|   |-- per_cell/
|   |   `-- <sample>_IV.xlsx
|   |-- <prefix>_IV_Cell_Data.csv
|   |-- <prefix>_IV_Group_Summary.csv
|   |-- <prefix>_IV_Supplementary_Metrics.csv
|   |-- <prefix>_IV_Apparent_Conductance.csv
|   `-- sample_registry.json
|-- 02_figures/
|   |-- <prefix>_IV_GroupMean_Traces.{png,pdf,svg}
|   |-- <prefix>_IV_Curve.{png,pdf,svg}
|   |-- <prefix>_IV_Supplementary_Metrics.{png,pdf,svg}
|   |-- <prefix>_Apparent_Conductance_Curve.{png,pdf,svg}
|   |-- <prefix>_Normalized_IV_Curve.{png,pdf,svg}  # only when enabled
|   `-- <prefix>_IV_Figure_Data.xlsx
|-- 03_statistics/
|   |-- <prefix>_IV_Statistics.xlsx
|   `-- <prefix>_IV_Results.json
|-- 04_report/
|   `-- <prefix>_IV_Report.md
|-- resolved_config.yaml
`-- analysis_manifest.json
```

## Primary tables

`Cell_Data` has one row per cell and voltage after repeated same-voltage sweeps have
been averaged. `endpoint_value` and `endpoint_unit` identify the configured primary
measurement. The raw, baseline-subtracted, and peak pA columns remain available
for audit.

`Group_Summary` contains cell count, mean, SD, and SEM for each group and voltage.
`Global_Test` contains at most one whole-curve two-group result. `Pointwise_Tests`
contains supplementary Welch comparisons and always retains `raw_p_value`; a Holm
column exists only when explicitly requested.

`Supplementary_Metrics` contains one row per cell for reversal potential, local
slope conductance, apparent input resistance, rectification, and current retention.
`Supplementary_Tests` reports group means, SEM, mean-difference confidence intervals,
raw p values, and Hedges' g. `Apparent_Conductance` contains the exploratory chord-
conductance curve after excluding voltages too close to each cell's interpolated
reversal potential.

`sample_registry.json` links each cell to the source ABF, SHA-256 hash, acquisition
metadata, selected group, and per-cell workbook. `analysis_manifest.json` records
the config, input hashes, generated output hashes, and software versions.

## Per-cell workbook

- `Raw_Traces`: full-resolution time and current traces for every sweep;
- `IV_Params`: command voltage, detected windows, baseline, steady-state, and peak
  endpoints per sweep;
- `Metadata`: source filename/hash, protocol, group, sample, sampling rate, units,
  and consensus step window.

Python-generated statistics are authoritative. The wide figure-data workbook is
provided for independent plotting or import into Origin; it is not a second
statistical analysis.
