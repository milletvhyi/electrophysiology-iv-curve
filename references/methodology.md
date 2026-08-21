# I-V analysis methodology

## Contents

- [Input and protocol gate](#input-and-protocol-gate)
- [Voltage-step detection](#voltage-step-detection)
- [Current endpoints](#current-endpoints)
- [Cell-level aggregation](#cell-level-aggregation)
- [Statistics](#statistics)
- [Interpretation boundaries](#interpretation-boundaries)

## Input and protocol gate

The workflow expects one immediate child directory per cell/sample. Group membership
is assigned from configured folder prefixes. A directory that contains ABFs but does
not match exactly one prefix is an error. Selection uses the configured regular
expression against both the acquisition protocol name and relative file path; file
size is never used as a protocol classifier.

Every included file must contain at least two sweeps, one recorded current channel,
and a command waveform in mV. The default workflow requires exactly one selected I-V
ABF per cell so that an accidental duplicate acquisition is not silently chosen.

## Voltage-step detection

For each sweep, the holding command is the median of the first baseline-window
samples. Command points differing from that value by more than
`command_threshold_mV` form candidate step segments. The longest segment meeting
`min_step_duration_ms` is retained. A consensus start and end are calculated across
sweeps and all detected segments must agree within
`step_alignment_tolerance_ms`.

The consensus window is used for every sweep. This is important when a test voltage
equals the holding voltage and therefore produces no command transition in that
sweep. The test voltage is the median command within the consensus window.

## Current endpoints

For every sweep the workflow records:

- baseline current: mean during the configured interval immediately before the
  voltage step;
- raw signed steady-state current: mean during the final configured interval of the
  voltage step;
- baseline-subtracted steady-state current: raw steady-state minus baseline;
- peak current: the largest absolute deviation from baseline after excluding the
  configured artifact guard at both step edges.

`raw_steady_state` is the compatibility default and retains the sign of the measured
current. It is not an absolute-value transform. `baseline_subtracted` is available
when it matches the prespecified experimental endpoint. The choice must be made
before interpreting group differences.

Current density modes divide the selected pA endpoint by capacitance in pF. They
require one unique, finite, positive capacitance value for every included cell.
Without that denominator, report pA and do not call it current density.

## Cell-level aggregation

Repeated sweeps at the same voltage are averaged within a cell. The resulting cell
profiles—not individual sweeps—enter group means, SEMs, plots, and statistics. All
cells must share the same voltage set for the compact group comparison. Trace panels
show group means calculated from every included cell at every voltage.

## Statistics

For exactly two groups with at least two cells each, the primary global statistic is
the mean squared difference between group curves after standardizing each voltage by
the across-cell standard deviation. Labels are permuted for whole cell profiles, so
the repeated voltage observations remain together.

All labelings are enumerated when their count does not exceed
`max_exact_labelings`. Otherwise a seeded Monte Carlo permutation test is used and
the plus-one correction is applied to its p value. The JSON and spreadsheet outputs
identify which method ran and how many labelings/permutations were evaluated.

Pointwise Welch tests are supplementary. Raw p values are the default. Optional Holm
family-wise correction is produced only when `multiple_comparison_correction: holm`
is explicitly set. With one group, more than two groups, fewer than two cells in a
group, or an incomplete cell-by-voltage matrix, inferential statistics are omitted
with a machine-readable reason.

Per-cell normalization divides the selected endpoint by the same cell's value at a
configured reference voltage. It is off by default, requires a non-negligible
denominator for every cell, and is labeled exploratory in every output.

## Exploratory cell-level metrics

When enabled, reversal potential is linearly interpolated at the best-supported
zero-current crossing in each cell. Local slope conductance is fitted from the
configured number of voltage points nearest that crossing; 1 pA/mV is reported as
1 nS. Apparent input resistance is 1000 divided by this slope and is omitted for
current-density endpoints. Rectification uses an observed symmetric voltage pair
nearest the configured target. Current retention uses the raw signed steady-state
to peak ratio at the selected positive voltage, defaulting to the largest positive
voltage recorded.

The apparent chord-conductance curve divides selected current by `V - Erev` and
omits points within the configured minimum driving force. It is an exploratory
algebraic transform, not evidence for channel identity or a valid single-channel
conductance model. These rules adapt to voltage count, spacing, sampling rate,
trace length, and detected step timing. A common voltage set remains required
within one group comparison so whole-cell profiles are comparable.

## Interpretation boundaries

An ABF-derived curve establishes the recorded current-voltage relationship under the
configured extraction rule. It does not by itself establish current identity,
channel mechanism, animal-level replication, successful blinding, or causal biology.
The generated report intentionally omits claims that are not traceable to the input
files, configuration, supplied capacitance metadata, or computed results.
