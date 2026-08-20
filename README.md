# Electrophysiology I-V Curve Skill

A standalone Codex skill for reproducible, cell-level I-V analysis of Axon Binary
Format (`.abf`) voltage-step recordings.

It preserves raw ABFs, detects the shared voltage-step window from command waveforms,
extracts raw signed and baseline-subtracted endpoints, handles repeated sweeps within
cells, and exports figures, CSV/Excel tables, JSON provenance, and a bounded Markdown
report. Two-group inference permutes intact cell curves; pointwise tests retain raw p
values by default.

## Install as a Codex skill

```powershell
git clone https://github.com/milletvhyi/electrophysiology-iv-curve.git "$env:USERPROFILE\.codex\skills\electrophysiology-iv-curve"
python -m pip install -r "$env:USERPROFILE\.codex\skills\electrophysiology-iv-curve\requirements.txt"
```

For project-scoped discovery, clone or copy the folder to
`<project>/.agents/skills/electrophysiology-iv-curve` instead.

## Quick start

1. Copy `assets/analysis_config.example.yaml` beside your data.
2. Edit the data/output paths, sample-folder prefixes, and protocol regex.
3. Run a read-only inventory first:

```powershell
python scripts/run_iv_analysis.py --config path\to\analysis_config.yaml --inventory-only
```

4. After confirming the selected files and protocol, run:

```powershell
python scripts/run_iv_analysis.py --config path\to\analysis_config.yaml
```

The input layout is one immediate child folder per cell. Each selected cell must have
exactly one intended I-V ABF. One-group datasets receive descriptive exports; exactly
two groups with adequate replication also receive a whole-curve permutation test and
supplementary pointwise Welch tests.

## Scientific defaults

- cells, never sweeps, are independent statistical units;
- all cells contribute to group curves and trace panels;
- raw signed steady-state current in pA is the default endpoint;
- current density requires supplied capacitance for every cell;
- normalized I-V is optional and explicitly exploratory;
- pointwise p values are unadjusted unless Holm correction is explicitly configured;
- unsupported biological, animal-count, and mechanism claims are omitted.

See [methodology](references/methodology.md) and the
[output schema](references/output-schema.md) for details.

## Validate the repository

```powershell
python -m unittest discover -s tests -v
python path\to\skill-creator\scripts\quick_validate.py .
```

This software supports research analysis and does not replace protocol review,
experimental quality control, or biological interpretation.

## 中文简介

该 skill 将 ABF 电压阶跃记录整理为可审计的细胞级 I-V 曲线。默认不把 sweep
当作独立样本，不自动做 Holm 校正，不把标准化 I-V 当成主分析，也不会根据文件
大小猜测协议。正式运行前请先执行 `--inventory-only`，核对入选文件、分组和电压
阶跃窗口。

## License

MIT
