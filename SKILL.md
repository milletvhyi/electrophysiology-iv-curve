---
name: electrophysiology-iv-curve
description: Analyze voltage-step patch-clamp recordings from Axon Binary Format (.abf) files and produce reproducible cell-level I-V curves, statistics, figures, spreadsheets, and a bounded Markdown report. Use for I-V or current-voltage analysis; do not use for MET/fluid-jet, action-potential, or other protocols.
---

# Electrophysiology I-V curve

Keep raw ABF files immutable and write every derived artifact to a separate output
directory. Treat each sample folder as one cell unless the user supplies a different,
traceable cell identifier.

## Workflow

1. Copy `assets/analysis_config.example.yaml` beside the dataset and edit paths,
   group prefixes/order, file-selection regex, extraction windows, and current mode.
2. Inventory the recordings before analysis:

```powershell
python "<skill-dir>\scripts\run_iv_analysis.py" --config "<analysis_config.yaml>" --inventory-only
```

3. Confirm that every included cell has exactly one intended voltage-step ABF, the
   command waveform contains a consistent step window, and the configured group and
   current measurement match the experiment.
4. Run the analysis into a new or empty output directory:

```powershell
python "<skill-dir>\scripts\run_iv_analysis.py" --config "<analysis_config.yaml>"
```

Use `--resume` only to regenerate the same configured run. It overwrites known
derived filenames but never deletes files.

Read `references/methodology.md` when choosing raw versus baseline-subtracted
current, current density, normalization, endpoint windows, or statistical options.
Read `references/output-schema.md` when interpreting or importing the exported
tables.

## Analysis rules

- Never infer an I-V protocol from ABF file size. Select by acquisition protocol or
  filename regex, then validate the command waveform.
- Use cells, not sweeps, as independent statistical units. Average repeated sweeps
  at the same voltage within each cell before any group comparison.
- Build trace panels and curves from every included cell. Do not select a
  representative cell unless the user explicitly asks for one and it is labeled as
  selected.
- Use raw signed steady-state current by default. Call a result current density only
  when a positive capacitance value is supplied for every included cell.
- Treat the raw/current-density I-V curve as primary. Enable per-cell normalization
  only as an explicitly labeled exploratory analysis; it does not replace the
  primary curve.
- When enabled, export exploratory reversal potential, local slope conductance,
  apparent input resistance, rectification, current retention, and an apparent
  conductance-voltage curve. Resolve targets from the observed voltage set; never
  assume 21 sweeps or a fixed -100 to +100 mV protocol. Label apparent conductance
  as non-channel-specific.
- With two adequately replicated groups, compare whole cell profiles by relabeling
  intact cell curves. Pointwise Welch tests are supplementary and report raw p
  values by default. Apply Holm correction only when explicitly configured.
- With one group or inadequate replication, export descriptive results without
  inventing an inferential comparison.
- Do not infer animal counts, biological replication, genotype, protocol details,
  or mechanisms that are absent from ABFs, configuration, or supplied metadata.

## Report

State the source and output directories, config path, included/excluded ABFs, cell
counts, voltage range, endpoint windows, current mode and units, group order, global
test method, pointwise correction policy, normalization status, and any validation
failure or manual exclusion. For every performed comparison, state its exact p
value and whether it meets the configured alpha; include group means, confidence
intervals, and effect size when available. Distinguish successful local execution
from biological or experimental validation.
