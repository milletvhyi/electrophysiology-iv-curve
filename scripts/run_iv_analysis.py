#!/usr/bin/env python
"""Reproducible cell-level I-V analysis for Axon ABF voltage-step recordings."""

from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
import math
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyabf
import scipy
import yaml
from scipy import stats


SCRIPT_VERSION = "1.0.0"
CURRENT_MEASURES = {
    "raw_steady_state": ("steady_state_raw_pA", "pA", False),
    "baseline_subtracted": ("steady_state_baseline_subtracted_pA", "pA", False),
    "current_density_raw_steady_state": (
        "steady_state_raw_pA",
        "pA/pF",
        True,
    ),
    "current_density_baseline_subtracted": (
        "steady_state_baseline_subtracted_pA",
        "pA/pF",
        True,
    ),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_filename(value: str) -> str:
    cleaned = re.sub(r'[<>:"/\\|?*\x00-\x1f]+', "_", str(value)).strip(" .")
    if not cleaned:
        raise ValueError("A project prefix or sample name produced an empty filename")
    return cleaned


def resolve_path(value: Any, config_path: Path) -> Path:
    if value is None or str(value).strip() == "":
        raise ValueError("A required path is empty")
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = config_path.parent / path
    return path.resolve()


def require_mapping(container: dict[str, Any], key: str) -> dict[str, Any]:
    value = container.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Config field '{key}' must be a mapping")
    return value


def require_positive(mapping: dict[str, Any], key: str) -> float:
    try:
        value = float(mapping[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Config field '{key}' must be a positive number") from exc
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"Config field '{key}' must be a positive number")
    return value


def load_and_validate_config(config_path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ValueError("The YAML root must be a mapping")

    project = require_mapping(raw, "project")
    paths = require_mapping(raw, "paths")
    files = require_mapping(raw, "files")
    extraction = require_mapping(raw, "extraction")
    analysis = require_mapping(raw, "analysis")

    project_name = str(project.get("name", "I-V experiment")).strip()
    prefix = safe_filename(str(project.get("prefix", project_name)))
    raw_root = resolve_path(paths.get("raw_data"), config_path)
    output_root = resolve_path(paths.get("output"), config_path)
    if not raw_root.is_dir():
        raise ValueError(f"Raw-data directory does not exist: {raw_root}")
    if (
        raw_root == output_root
        or raw_root in output_root.parents
        or output_root in raw_root.parents
    ):
        raise ValueError("Output and raw-data directories must be separate and non-nested")

    raw_groups = raw.get("groups")
    if not isinstance(raw_groups, list) or not raw_groups:
        raise ValueError("Config field 'groups' must be a non-empty list")
    groups: list[dict[str, str]] = []
    seen_names: set[str] = set()
    seen_prefixes: set[str] = set()
    for index, item in enumerate(raw_groups):
        if not isinstance(item, dict):
            raise ValueError(f"groups[{index}] must be a mapping")
        folder_prefix = str(item.get("folder_prefix", ""))
        name = str(item.get("name", "")).strip()
        role = str(item.get("role", "unspecified")).strip() or "unspecified"
        color = str(item.get("color", "")).strip()
        if not folder_prefix or not name:
            raise ValueError(f"groups[{index}] needs folder_prefix and name")
        if folder_prefix in seen_prefixes or name in seen_names:
            raise ValueError("Group prefixes and names must be unique")
        if not re.fullmatch(r"#[0-9A-Fa-f]{6}", color):
            raise ValueError(f"Invalid HEX color for group {name}: {color}")
        seen_prefixes.add(folder_prefix)
        seen_names.add(name)
        groups.append(
            {
                "folder_prefix": folder_prefix,
                "name": name,
                "role": role,
                "color": color,
            }
        )

    include_regex = files.get("include_regex")
    accept_all = bool(files.get("accept_all_abf", False))
    if not accept_all and not isinstance(include_regex, str):
        raise ValueError("files.include_regex is required unless accept_all_abf is true")
    try:
        include_pattern = re.compile(include_regex) if include_regex else None
        exclude_pattern = (
            re.compile(str(files["exclude_regex"]))
            if files.get("exclude_regex")
            else None
        )
    except re.error as exc:
        raise ValueError(f"Invalid file-selection regular expression: {exc}") from exc

    try:
        channel = int(extraction.get("channel", 0))
        voltage_round_decimals = int(extraction.get("voltage_round_decimals", 3))
    except (TypeError, ValueError) as exc:
        raise ValueError("extraction.channel and voltage_round_decimals must be integers") from exc
    if channel < 0:
        raise ValueError("extraction.channel must be zero or greater")
    if not 0 <= voltage_round_decimals <= 6:
        raise ValueError("voltage_round_decimals must be between 0 and 6")

    extraction_resolved = {
        "channel": channel,
        "command_threshold_mV": require_positive(extraction, "command_threshold_mV"),
        "min_step_duration_ms": require_positive(extraction, "min_step_duration_ms"),
        "step_alignment_tolerance_ms": require_positive(
            extraction, "step_alignment_tolerance_ms"
        ),
        "baseline_window_ms": require_positive(extraction, "baseline_window_ms"),
        "steady_state_window_ms": require_positive(
            extraction, "steady_state_window_ms"
        ),
        "artifact_guard_ms": require_positive(extraction, "artifact_guard_ms"),
        "voltage_round_decimals": voltage_round_decimals,
    }

    if str(analysis.get("statistical_unit", "cell")).lower() != "cell":
        raise ValueError("This skill requires analysis.statistical_unit: cell")
    current_measure = str(analysis.get("current_measure", "raw_steady_state")).lower()
    if current_measure not in CURRENT_MEASURES:
        raise ValueError(
            "analysis.current_measure must be one of "
            + ", ".join(CURRENT_MEASURES)
        )
    correction = str(
        analysis.get("multiple_comparison_correction", "none")
    ).lower()
    if correction not in {"none", "holm"}:
        raise ValueError("multiple_comparison_correction must be none or holm")
    try:
        max_exact = int(analysis.get("max_exact_labelings", 100000))
        monte_carlo = int(analysis.get("monte_carlo_permutations", 10000))
        random_seed = int(analysis.get("random_seed", 20260820))
    except (TypeError, ValueError) as exc:
        raise ValueError("Permutation settings must be integers") from exc
    if max_exact < 1 or monte_carlo < 1:
        raise ValueError("Permutation counts must be positive")

    normalization = analysis.get("normalization", {})
    if not isinstance(normalization, dict):
        raise ValueError("analysis.normalization must be a mapping")
    normalization_resolved = {
        "enabled": bool(normalization.get("enabled", False)),
        "reference_voltage_mV": float(
            normalization.get("reference_voltage_mV", 100.0)
        ),
        "minimum_denominator_abs": float(
            normalization.get("minimum_denominator_abs", 1.0)
        ),
    }
    if normalization_resolved["minimum_denominator_abs"] <= 0:
        raise ValueError("normalization.minimum_denominator_abs must be positive")

    metadata_path = None
    if analysis.get("cell_metadata_csv"):
        metadata_path = resolve_path(analysis["cell_metadata_csv"], config_path)
        if not metadata_path.is_file():
            raise ValueError(f"Cell metadata CSV does not exist: {metadata_path}")
    if CURRENT_MEASURES[current_measure][2] and metadata_path is None:
        raise ValueError("Current-density analysis requires analysis.cell_metadata_csv")

    analysis_resolved = {
        "statistical_unit": "cell",
        "current_measure": current_measure,
        "cell_metadata_csv": str(metadata_path) if metadata_path else None,
        "sample_column": str(analysis.get("sample_column", "sample")),
        "capacitance_column": str(
            analysis.get("capacitance_column", "capacitance_pF")
        ),
        "pointwise_tests": bool(analysis.get("pointwise_tests", True)),
        "multiple_comparison_correction": correction,
        "max_exact_labelings": max_exact,
        "monte_carlo_permutations": monte_carlo,
        "random_seed": random_seed,
        "normalization": normalization_resolved,
    }

    resolved = copy.deepcopy(raw)
    resolved["project"] = {"name": project_name, "prefix": prefix}
    resolved["paths"] = {"raw_data": str(raw_root), "output": str(output_root)}
    resolved["groups"] = groups
    resolved["files"] = {
        "include_regex": include_regex,
        "exclude_regex": files.get("exclude_regex"),
        "accept_all_abf": accept_all,
    }
    resolved["extraction"] = extraction_resolved
    resolved["analysis"] = analysis_resolved

    runtime = {
        "project_name": project_name,
        "prefix": prefix,
        "raw_root": raw_root,
        "output_root": output_root,
        "groups": groups,
        "include_pattern": include_pattern,
        "exclude_pattern": exclude_pattern,
        "accept_all": accept_all,
        "extraction": extraction_resolved,
        "analysis": analysis_resolved,
        "metadata_path": metadata_path,
    }
    return resolved, runtime


def abf_header(path: Path) -> dict[str, Any]:
    abf = pyabf.ABF(str(path), loadData=False)
    return {
        "protocol": str(abf.protocol or ""),
        "sweep_count": int(abf.sweepCount),
        "channel_count": int(abf.channelCount),
        "sampling_rate_Hz": int(abf.dataRate),
        "points_per_sweep": int(abf.sweepPointCount),
        "adc_units": [str(item) for item in getattr(abf, "adcUnits", [])],
        "dac_units": [str(item) for item in getattr(abf, "dacUnits", [])],
        "size_bytes": int(path.stat().st_size),
    }


def group_for_sample(sample: str, groups: list[dict[str, str]]) -> dict[str, str]:
    matches = [item for item in groups if sample.startswith(item["folder_prefix"])]
    if len(matches) != 1:
        prefixes = [item["folder_prefix"] for item in matches]
        raise ValueError(
            f"Sample folder {sample!r} matched {len(matches)} group prefixes: {prefixes}"
        )
    return matches[0]


def discover_files(runtime: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    raw_root: Path = runtime["raw_root"]
    direct_abfs = [
        item for item in raw_root.iterdir() if item.is_file() and item.suffix.lower() == ".abf"
    ]
    if direct_abfs:
        raise ValueError(
            "ABFs were found directly under raw_data; create one child folder per cell: "
            + ", ".join(item.name for item in direct_abfs[:5])
        )

    selected: list[dict[str, Any]] = []
    inventory: list[dict[str, Any]] = []
    for sample_dir in sorted(
        (item for item in raw_root.iterdir() if item.is_dir()),
        key=lambda item: item.name.lower(),
    ):
        abfs = sorted(
            (
                item
                for item in sample_dir.rglob("*")
                if item.is_file() and item.suffix.lower() == ".abf"
            ),
            key=lambda item: item.as_posix().lower(),
        )
        if not abfs:
            continue
        group = group_for_sample(sample_dir.name, runtime["groups"])
        sample_selected: list[dict[str, Any]] = []
        for path in abfs:
            header = abf_header(path)
            relative = path.relative_to(raw_root).as_posix()
            search_text = relative + "\n" + header["protocol"]
            included = runtime["accept_all"] or bool(
                runtime["include_pattern"]
                and runtime["include_pattern"].search(search_text)
            )
            excluded = bool(
                runtime["exclude_pattern"]
                and runtime["exclude_pattern"].search(search_text)
            )
            is_selected = included and not excluded
            record = {
                "sample": sample_dir.name,
                "group": group["name"],
                "role": group["role"],
                "relative_path": relative,
                **header,
                "selected": is_selected,
                "selection_reason": (
                    "excluded_by_regex"
                    if excluded
                    else "accept_all_abf"
                    if runtime["accept_all"]
                    else "include_regex_match"
                    if included
                    else "no_include_regex_match"
                ),
            }
            inventory.append(record)
            if is_selected:
                chosen = {**record, "path": path, "group_spec": group}
                selected.append(chosen)
                sample_selected.append(chosen)
        if len(sample_selected) != 1:
            names = [item["relative_path"] for item in sample_selected]
            raise ValueError(
                f"Sample {sample_dir.name!r} must have exactly one selected I-V ABF; "
                f"found {len(sample_selected)}: {names}"
            )

    if not selected:
        raise ValueError("No I-V ABFs were selected")
    observed = {item["group"] for item in selected}
    empty_groups = [item["name"] for item in runtime["groups"] if item["name"] not in observed]
    if empty_groups:
        raise ValueError(f"Configured groups with no selected cells: {empty_groups}")
    return selected, inventory


def normalize_unit(unit: str) -> str:
    return str(unit).strip().replace("μ", "u").replace("µ", "u").lower()


def current_to_pa_factor(unit: str) -> float:
    factors = {"pa": 1.0, "na": 1e3, "ua": 1e6, "ma": 1e9, "a": 1e12}
    key = normalize_unit(unit)
    if key not in factors:
        raise ValueError(f"Recorded channel is not in a recognized current unit: {unit!r}")
    return factors[key]


def command_to_mv_factor(unit: str) -> float:
    factors = {"uv": 1e-3, "mv": 1.0, "v": 1e3}
    key = normalize_unit(unit)
    if key not in factors:
        raise ValueError(f"Command channel is not in a recognized voltage unit: {unit!r}")
    return factors[key]


def contiguous_segments(mask: np.ndarray) -> list[tuple[int, int]]:
    edges = np.diff(np.r_[False, mask, False].astype(np.int8))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]
    return [(int(start), int(end)) for start, end in zip(starts, ends)]


def detect_step_window(
    commands_mV: np.ndarray,
    sampling_rate_Hz: float,
    extraction: dict[str, Any],
) -> dict[str, Any]:
    commands = np.asarray(commands_mV, dtype=float)
    if commands.ndim != 2 or commands.shape[0] < 2:
        raise ValueError("An I-V recording must contain at least two command sweeps")
    baseline_n = max(
        1, int(round(extraction["baseline_window_ms"] * sampling_rate_Hz / 1000))
    )
    min_step_n = max(
        1, int(round(extraction["min_step_duration_ms"] * sampling_rate_Hz / 1000))
    )
    tolerance_n = max(
        1,
        int(round(extraction["step_alignment_tolerance_ms"] * sampling_rate_Hz / 1000)),
    )
    if baseline_n >= commands.shape[1]:
        raise ValueError("Baseline window is longer than the sweep")

    holds: list[float] = []
    detected: list[tuple[int, int] | None] = []
    usable: list[tuple[int, int]] = []
    for command in commands:
        hold = float(np.median(command[:baseline_n]))
        holds.append(hold)
        mask = np.abs(command - hold) > extraction["command_threshold_mV"]
        segments = [
            segment
            for segment in contiguous_segments(mask)
            if segment[1] - segment[0] >= min_step_n
        ]
        if not segments:
            detected.append(None)
            continue
        longest = max(segments, key=lambda item: item[1] - item[0])
        detected.append(longest)
        usable.append(longest)

    if not usable:
        raise ValueError("No voltage-step command segment was detected in any sweep")
    consensus_start = int(round(float(np.median([item[0] for item in usable]))))
    consensus_end = int(round(float(np.median([item[1] for item in usable]))))
    for start, end in usable:
        if (
            abs(start - consensus_start) > tolerance_n
            or abs(end - consensus_end) > tolerance_n
        ):
            raise ValueError(
                "Detected voltage-step windows are not aligned across sweeps: "
                f"consensus={consensus_start}:{consensus_end}, observed={start}:{end}"
            )
    if consensus_start < baseline_n:
        raise ValueError("The configured baseline window does not fit before the step")
    if consensus_end <= consensus_start or consensus_end > commands.shape[1]:
        raise ValueError("Detected voltage-step window is invalid")

    steady_n = max(
        1, int(round(extraction["steady_state_window_ms"] * sampling_rate_Hz / 1000))
    )
    guard_n = max(
        0, int(round(extraction["artifact_guard_ms"] * sampling_rate_Hz / 1000))
    )
    if steady_n > consensus_end - consensus_start:
        raise ValueError("Steady-state window is longer than the voltage step")
    if 2 * guard_n >= consensus_end - consensus_start:
        raise ValueError("Artifact guards consume the entire voltage step")
    return {
        "start_index": consensus_start,
        "end_index": consensus_end,
        "baseline_samples": baseline_n,
        "steady_samples": steady_n,
        "artifact_guard_samples": guard_n,
        "holding_mV": holds,
        "detected_segments": detected,
    }


def extract_recording(
    selected: dict[str, Any], extraction: dict[str, Any], source_hash: str
) -> dict[str, Any]:
    path: Path = selected["path"]
    abf = pyabf.ABF(str(path))
    channel = extraction["channel"]
    if channel >= abf.channelCount:
        raise ValueError(
            f"{path.name}: channel {channel} is unavailable; channel_count={abf.channelCount}"
        )
    adc_units = [str(item) for item in getattr(abf, "adcUnits", [])]
    dac_units = [str(item) for item in getattr(abf, "dacUnits", [])]
    if channel >= len(adc_units) or channel >= len(dac_units):
        raise ValueError(f"{path.name}: ABF channel units are unavailable")
    current_factor = current_to_pa_factor(adc_units[channel])
    command_factor = command_to_mv_factor(dac_units[channel])

    times: list[np.ndarray] = []
    currents: list[np.ndarray] = []
    commands: list[np.ndarray] = []
    for sweep in range(abf.sweepCount):
        abf.setSweep(sweep, channel=channel)
        times.append(np.array(abf.sweepX, dtype=float, copy=True) * 1000)
        currents.append(
            np.array(abf.sweepY, dtype=float, copy=True) * current_factor
        )
        commands.append(
            np.array(abf.sweepC, dtype=float, copy=True) * command_factor
        )
    lengths = {len(item) for item in times + currents + commands}
    if len(lengths) != 1:
        raise ValueError(f"{path.name}: sweep lengths are inconsistent")
    reference_time = times[0]
    for time in times[1:]:
        if not np.allclose(time, reference_time, rtol=0, atol=1e-9):
            raise ValueError(f"{path.name}: time bases differ across sweeps")

    current_array = np.vstack(currents)
    command_array = np.vstack(commands)
    window = detect_step_window(command_array, abf.dataRate, extraction)
    start = window["start_index"]
    end = window["end_index"]
    baseline_n = window["baseline_samples"]
    steady_n = window["steady_samples"]
    guard_n = window["artifact_guard_samples"]
    decimals = extraction["voltage_round_decimals"]

    rows: list[dict[str, Any]] = []
    voltages: list[float] = []
    for sweep in range(abf.sweepCount):
        current = current_array[sweep]
        command = command_array[sweep]
        voltage = round(float(np.median(command[start:end])), decimals)
        voltages.append(voltage)
        baseline = float(np.mean(current[start - baseline_n : start]))
        steady_raw = float(np.mean(current[end - steady_n : end]))
        steady_delta = steady_raw - baseline
        peak_region = current[start + guard_n : end - guard_n if guard_n else end]
        peak_index_local = int(np.argmax(np.abs(peak_region - baseline)))
        peak_raw = float(peak_region[peak_index_local])
        peak_delta = peak_raw - baseline
        detected = window["detected_segments"][sweep]
        rows.append(
            {
                "sweep": sweep,
                "voltage_mV": voltage,
                "holding_mV": window["holding_mV"][sweep],
                "detected_transition": detected is not None,
                "detected_step_start_index": detected[0] if detected else np.nan,
                "detected_step_end_index": detected[1] if detected else np.nan,
                "consensus_step_start_index": start,
                "consensus_step_end_index": end,
                "baseline_start_index": start - baseline_n,
                "baseline_end_index": start,
                "steady_state_start_index": end - steady_n,
                "steady_state_end_index": end,
                "baseline_current_pA": baseline,
                "steady_state_raw_pA": steady_raw,
                "steady_state_baseline_subtracted_pA": steady_delta,
                "peak_raw_pA": peak_raw,
                "peak_baseline_subtracted_pA": peak_delta,
            }
        )

    voltage_trace_means: dict[float, np.ndarray] = {}
    voltage_array = np.asarray(voltages, dtype=float)
    for voltage in sorted(set(voltages)):
        voltage_trace_means[voltage] = np.mean(
            current_array[np.isclose(voltage_array, voltage)], axis=0
        )

    return {
        "sample": selected["sample"],
        "group": selected["group"],
        "role": selected["role"],
        "color": selected["group_spec"]["color"],
        "path": path,
        "relative_path": selected["relative_path"],
        "source_sha256": source_hash,
        "protocol": str(abf.protocol or ""),
        "sampling_rate_Hz": int(abf.dataRate),
        "sweep_count": int(abf.sweepCount),
        "current_input_unit": adc_units[channel],
        "command_input_unit": dac_units[channel],
        "time_ms": reference_time,
        "currents_pA": current_array,
        "commands_mV": command_array,
        "voltages_mV": voltages,
        "voltage_trace_means": voltage_trace_means,
        "params": pd.DataFrame(rows),
        "step_start_index": start,
        "step_end_index": end,
        "step_start_ms": float(reference_time[start]),
        "step_end_ms": float(reference_time[end - 1]),
    }


def load_capacitance(
    runtime: dict[str, Any], samples: Iterable[str]
) -> dict[str, float]:
    path: Path | None = runtime["metadata_path"]
    if path is None:
        return {}
    analysis = runtime["analysis"]
    table = pd.read_csv(path)
    sample_column = analysis["sample_column"]
    capacitance_column = analysis["capacitance_column"]
    missing_columns = [
        item for item in (sample_column, capacitance_column) if item not in table.columns
    ]
    if missing_columns:
        raise ValueError(f"Cell metadata CSV is missing columns: {missing_columns}")
    table = table[[sample_column, capacitance_column]].copy()
    table[sample_column] = table[sample_column].astype(str)
    if table[sample_column].duplicated().any():
        duplicate = table.loc[table[sample_column].duplicated(), sample_column].tolist()
        raise ValueError(f"Duplicate sample rows in cell metadata CSV: {duplicate}")
    values = dict(zip(table[sample_column], table[capacitance_column]))
    result: dict[str, float] = {}
    for sample in samples:
        if sample not in values:
            raise ValueError(f"Cell metadata CSV has no row for sample {sample!r}")
        try:
            capacitance = float(values[sample])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid capacitance for sample {sample!r}") from exc
        if not np.isfinite(capacitance) or capacitance <= 0:
            raise ValueError(f"Capacitance must be positive for sample {sample!r}")
        result[sample] = capacitance
    return result


def validate_voltage_sets(extractions: list[dict[str, Any]]) -> list[float]:
    reference = sorted(set(extractions[0]["voltages_mV"]))
    mismatches: list[str] = []
    for item in extractions[1:]:
        observed = sorted(set(item["voltages_mV"]))
        if observed != reference:
            mismatches.append(f"{item['sample']}: {observed}")
    if mismatches:
        raise ValueError(
            "All cells must share the same voltage set. Reference "
            f"{extractions[0]['sample']}: {reference}; mismatches: {mismatches[:5]}"
        )
    return reference


def build_cell_data(
    extractions: list[dict[str, Any]],
    capacitance: dict[str, float],
    current_measure: str,
) -> pd.DataFrame:
    source_column, endpoint_unit, density = CURRENT_MEASURES[current_measure]
    rows: list[dict[str, Any]] = []
    numeric_columns = [
        "baseline_current_pA",
        "steady_state_raw_pA",
        "steady_state_baseline_subtracted_pA",
        "peak_raw_pA",
        "peak_baseline_subtracted_pA",
    ]
    for item in extractions:
        grouped = item["params"].groupby("voltage_mV", sort=True)
        for voltage, frame in grouped:
            row: dict[str, Any] = {
                "sample": item["sample"],
                "group": item["group"],
                "role": item["role"],
                "color": item["color"],
                "voltage_mV": float(voltage),
                "sweep_count_at_voltage": int(len(frame)),
                "capacitance_pF": capacitance.get(item["sample"], np.nan),
                "current_measure": current_measure,
                "endpoint_unit": endpoint_unit,
            }
            for column in numeric_columns:
                row[column] = float(frame[column].mean())
            endpoint = row[source_column]
            if density:
                endpoint = endpoint / row["capacitance_pF"]
            row["endpoint_value"] = float(endpoint)
            rows.append(row)
    return pd.DataFrame(rows)


def group_summary(
    cell_data: pd.DataFrame,
    group_order: list[str],
    value_column: str,
    value_unit: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for group in group_order:
        group_frame = cell_data.loc[cell_data["group"] == group]
        for voltage in sorted(group_frame["voltage_mV"].unique()):
            values = group_frame.loc[
                np.isclose(group_frame["voltage_mV"], voltage), value_column
            ].to_numpy(float)
            rows.append(
                {
                    "group": group,
                    "voltage_mV": float(voltage),
                    "n_cells": int(len(values)),
                    "mean": float(np.mean(values)),
                    "sd": float(np.std(values, ddof=1)) if len(values) > 1 else np.nan,
                    "sem": float(stats.sem(values)) if len(values) > 1 else np.nan,
                    "unit": value_unit,
                }
            )
    return pd.DataFrame(rows)


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    df = len(a) + len(b) - 2
    if df <= 0:
        return np.nan
    pooled_variance = (
        (len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)
    ) / df
    if not np.isfinite(pooled_variance) or pooled_variance <= 0:
        return np.nan
    correction = 1 - 3 / (4 * df - 1)
    return float(correction * (np.mean(a) - np.mean(b)) / np.sqrt(pooled_variance))


def welch_summary(a: np.ndarray, b: np.ndarray) -> dict[str, Any]:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    result = stats.ttest_ind(a, b, equal_var=False)
    variance_a = np.var(a, ddof=1) / len(a)
    variance_b = np.var(b, ddof=1) / len(b)
    difference = float(np.mean(a) - np.mean(b))
    denominator = variance_a**2 / (len(a) - 1) + variance_b**2 / (len(b) - 1)
    df = (variance_a + variance_b) ** 2 / denominator if denominator > 0 else np.nan
    se_difference = np.sqrt(variance_a + variance_b)
    critical = stats.t.ppf(0.975, df) if np.isfinite(df) else np.nan
    return {
        "n_group_1": int(len(a)),
        "n_group_2": int(len(b)),
        "mean_group_1": float(np.mean(a)),
        "sem_group_1": float(stats.sem(a)),
        "mean_group_2": float(np.mean(b)),
        "sem_group_2": float(stats.sem(b)),
        "mean_difference_group_1_minus_2": difference,
        "ci95_low": float(difference - critical * se_difference),
        "ci95_high": float(difference + critical * se_difference),
        "welch_t": float(result.statistic),
        "welch_df": float(df),
        "raw_p_value": float(result.pvalue),
        "hedges_g": hedges_g(a, b),
    }


def holm_adjust(p_values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(p_values), dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted_ranked = np.maximum.accumulate(
        (len(ranked) - np.arange(len(ranked))) * ranked
    )
    adjusted_ranked = np.clip(adjusted_ranked, 0, 1)
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = adjusted_ranked
    return adjusted


def curve_statistic(matrix: np.ndarray, group_1_indices: np.ndarray) -> float:
    matrix = np.asarray(matrix, dtype=float)
    group_1_indices = np.asarray(group_1_indices, dtype=int)
    mask = np.zeros(matrix.shape[0], dtype=bool)
    mask[group_1_indices] = True
    scale = np.std(matrix, axis=0, ddof=1)
    scale[~np.isfinite(scale) | (scale == 0)] = 1.0
    difference = matrix[mask].mean(axis=0) - matrix[~mask].mean(axis=0)
    return float(np.mean((difference / scale) ** 2))


def curve_permutation_test(
    matrix: np.ndarray,
    labels: np.ndarray,
    group_1: str,
    group_2: str,
    max_exact_labelings: int,
    monte_carlo_permutations: int,
    random_seed: int,
) -> dict[str, Any]:
    matrix = np.asarray(matrix, dtype=float)
    labels = np.asarray(labels)
    n_group_1 = int(np.sum(labels == group_1))
    observed_indices = np.where(labels == group_1)[0]
    observed = curve_statistic(matrix, observed_indices)
    total_labelings = math.comb(matrix.shape[0], n_group_1)

    if total_labelings <= max_exact_labelings:
        method = "exact_label_permutation"
        evaluated = total_labelings
        extreme = 0
        for indices in itertools.combinations(range(matrix.shape[0]), n_group_1):
            statistic = curve_statistic(matrix, np.asarray(indices, dtype=int))
            extreme += int(statistic >= observed - 1e-12)
        p_value = extreme / evaluated
        plus_one = False
    else:
        method = "monte_carlo_label_permutation"
        evaluated = monte_carlo_permutations
        extreme = 0
        rng = np.random.default_rng(random_seed)
        for _ in range(evaluated):
            indices = rng.permutation(matrix.shape[0])[:n_group_1]
            statistic = curve_statistic(matrix, indices)
            extreme += int(statistic >= observed - 1e-12)
        p_value = (extreme + 1) / (evaluated + 1)
        plus_one = True

    return {
        "status": "performed",
        "method": method,
        "statistic_standardized_mean_square": observed,
        "p_value": float(p_value),
        "extreme_or_equal": int(extreme),
        "evaluated_labelings_or_permutations": int(evaluated),
        "total_possible_labelings": int(total_labelings),
        "plus_one_correction": plus_one,
        "random_seed": random_seed if method.startswith("monte_carlo") else None,
        "group_1": group_1,
        "group_2": group_2,
        "n_group_1_cells": n_group_1,
        "n_group_2_cells": int(np.sum(labels == group_2)),
        "statistical_unit": "cell",
    }


def run_inference(
    cell_data: pd.DataFrame,
    value_column: str,
    group_order: list[str],
    analysis: dict[str, Any],
) -> tuple[dict[str, Any], pd.DataFrame]:
    observed_groups = [
        group for group in group_order if group in set(cell_data["group"].astype(str))
    ]
    if len(observed_groups) != 2:
        return (
            {
                "status": "skipped",
                "reason": "Exactly two groups are required for the compact inferential analysis",
                "observed_groups": observed_groups,
            },
            pd.DataFrame(),
        )
    group_1, group_2 = observed_groups
    counts = cell_data[["sample", "group"]].drop_duplicates()["group"].value_counts()
    if counts.get(group_1, 0) < 2 or counts.get(group_2, 0) < 2:
        return (
            {
                "status": "skipped",
                "reason": "At least two cells per group are required",
                "group_cell_counts": {item: int(counts.get(item, 0)) for item in observed_groups},
            },
            pd.DataFrame(),
        )

    pivot = cell_data.pivot(
        index=["sample", "group"], columns="voltage_mV", values=value_column
    ).sort_index(axis=1)
    if pivot.isna().any().any():
        return (
            {
                "status": "skipped",
                "reason": "The cell-by-voltage matrix is incomplete",
            },
            pd.DataFrame(),
        )
    labels = np.asarray([index[1] for index in pivot.index], dtype=object)
    global_result = curve_permutation_test(
        pivot.to_numpy(float),
        labels,
        group_1,
        group_2,
        analysis["max_exact_labelings"],
        analysis["monte_carlo_permutations"],
        analysis["random_seed"],
    )

    if not analysis["pointwise_tests"]:
        return global_result, pd.DataFrame()
    pointwise_rows: list[dict[str, Any]] = []
    for voltage in pivot.columns:
        a = pivot.xs(group_1, level="group")[voltage].to_numpy(float)
        b = pivot.xs(group_2, level="group")[voltage].to_numpy(float)
        pointwise_rows.append(
            {
                "voltage_mV": float(voltage),
                "group_1": group_1,
                "group_2": group_2,
                **welch_summary(a, b),
            }
        )
    pointwise = pd.DataFrame(pointwise_rows)
    if analysis["multiple_comparison_correction"] == "holm":
        pointwise["holm_adjusted_p_value"] = holm_adjust(
            pointwise["raw_p_value"].to_numpy(float)
        )
    return global_result, pointwise


def normalize_cell_data(
    cell_data: pd.DataFrame, reference_voltage_mV: float, minimum_abs: float
) -> pd.DataFrame:
    voltages = sorted(cell_data["voltage_mV"].unique())
    matches = [value for value in voltages if np.isclose(value, reference_voltage_mV)]
    if len(matches) != 1:
        raise ValueError(
            f"Normalization reference {reference_voltage_mV:g} mV is not uniquely present; "
            f"available voltages: {voltages}"
        )
    reference = matches[0]
    denominators = (
        cell_data.loc[np.isclose(cell_data["voltage_mV"], reference)]
        .set_index("sample")["endpoint_value"]
        .to_dict()
    )
    bad = {
        sample: value
        for sample, value in denominators.items()
        if not np.isfinite(value) or abs(value) < minimum_abs
    }
    if bad:
        raise ValueError(f"Normalization denominators are invalid or too small: {bad}")
    result = cell_data[
        ["sample", "group", "role", "color", "voltage_mV"]
    ].copy()
    result["reference_voltage_mV"] = float(reference)
    result["normalization_denominator"] = result["sample"].map(denominators)
    result["normalized_value"] = (
        cell_data["endpoint_value"].to_numpy(float)
        / result["normalization_denominator"].to_numpy(float)
    )
    result["analysis_role"] = "exploratory"
    return result


def prepare_output(output: Path, resume: bool) -> dict[str, Path]:
    if output.exists() and any(output.iterdir()) and not resume:
        raise ValueError(
            f"Output directory is not empty: {output}. Use a new directory or --resume."
        )
    directories = {
        "root": output,
        "extracted": output / "01_extracted",
        "per_cell": output / "01_extracted" / "per_cell",
        "figures": output / "02_figures",
        "statistics": output / "03_statistics",
        "report": output / "04_report",
    }
    for path in directories.values():
        path.mkdir(parents=True, exist_ok=True)
    return directories


def write_per_cell_workbook(
    item: dict[str, Any], output_path: Path, capacitance_pF: float | None
) -> None:
    raw = pd.DataFrame({"Time_ms": item["time_ms"]})
    for sweep, voltage in enumerate(item["voltages_mV"]):
        raw[f"sweep_{sweep:03d}_{voltage:+g}mV_pA"] = item["currents_pA"][sweep]
    metadata = pd.DataFrame(
        [
            {
                "sample": item["sample"],
                "group": item["group"],
                "role": item["role"],
                "source_abf": item["relative_path"],
                "source_sha256": item["source_sha256"],
                "protocol": item["protocol"],
                "sampling_rate_Hz": item["sampling_rate_Hz"],
                "sweep_count": item["sweep_count"],
                "input_current_unit": item["current_input_unit"],
                "input_command_unit": item["command_input_unit"],
                "exported_current_unit": "pA",
                "exported_command_unit": "mV",
                "consensus_step_start_index": item["step_start_index"],
                "consensus_step_end_index": item["step_end_index"],
                "consensus_step_start_ms": item["step_start_ms"],
                "consensus_step_end_ms": item["step_end_ms"],
                "capacitance_pF": capacitance_pF,
            }
        ]
    )
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        raw.to_excel(writer, sheet_name="Raw_Traces", index=False)
        item["params"].to_excel(writer, sheet_name="IV_Params", index=False)
        metadata.to_excel(writer, sheet_name="Metadata", index=False)


def build_trace_table(
    extractions: list[dict[str, Any]], group_order: list[str], voltages: list[float]
) -> pd.DataFrame:
    reference_time = extractions[0]["time_ms"]
    for item in extractions[1:]:
        if len(item["time_ms"]) != len(reference_time) or not np.allclose(
            item["time_ms"], reference_time, rtol=0, atol=1e-6
        ):
            raise ValueError("Group-mean trace plotting requires a common time base")
    table = pd.DataFrame({"Time_ms": reference_time})
    for group in group_order:
        members = [item for item in extractions if item["group"] == group]
        for voltage in voltages:
            traces = [item["voltage_trace_means"][voltage] for item in members]
            table[f"{group}|{voltage:+g}mV|mean_pA"] = np.mean(traces, axis=0)
    return table


def style_axis(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(direction="out", length=4, width=1)


def save_figure(figure: plt.Figure, directory: Path, stem: str) -> list[str]:
    paths: list[str] = []
    for suffix, kwargs in (
        ("png", {"dpi": 300}),
        ("pdf", {}),
        ("svg", {}),
    ):
        path = directory / f"{stem}.{suffix}"
        figure.savefig(path, bbox_inches="tight", **kwargs)
        paths.append(str(path))
    plt.close(figure)
    return paths


def plot_group_mean_traces(
    trace_table: pd.DataFrame,
    groups: list[dict[str, str]],
    voltages: list[float],
    directory: Path,
    prefix: str,
) -> list[str]:
    figure, axes = plt.subplots(
        len(groups), 1, figsize=(7.2, max(3.1, 2.65 * len(groups))), sharex=True, sharey=True
    )
    axes_array = np.atleast_1d(axes)
    for axis, group in zip(axes_array, groups):
        color = group["color"]
        for index, voltage in enumerate(voltages):
            alpha = 0.25 + 0.75 * index / max(1, len(voltages) - 1)
            column = f"{group['name']}|{voltage:+g}mV|mean_pA"
            axis.plot(trace_table["Time_ms"], trace_table[column], color=color, alpha=alpha, lw=0.9)
        axis.axhline(0, color="#777777", lw=0.6, alpha=0.6)
        axis.set_ylabel("Current (pA)")
        axis.set_title(
            f"{group['name']} ({group['role']}): all-cell group means; "
            f"{voltages[0]:g} to {voltages[-1]:g} mV",
            fontsize=10,
        )
        style_axis(axis)
    axes_array[-1].set_xlabel("Time (ms)")
    figure.suptitle("Group-mean voltage-step current families", y=1.01)
    figure.tight_layout()
    return save_figure(figure, directory, f"{prefix}_IV_GroupMean_Traces")


def plot_iv_curve(
    cell_data: pd.DataFrame,
    summary: pd.DataFrame,
    groups: list[dict[str, str]],
    value_column: str,
    ylabel: str,
    title: str,
    global_result: dict[str, Any],
    directory: Path,
    stem: str,
) -> list[str]:
    figure, axis = plt.subplots(figsize=(7.0, 5.2))
    for group in groups:
        frame = cell_data.loc[cell_data["group"] == group["name"]]
        for _, cell in frame.groupby("sample", sort=True):
            cell = cell.sort_values("voltage_mV")
            axis.plot(
                cell["voltage_mV"],
                cell[value_column],
                color=group["color"],
                alpha=0.18,
                lw=0.8,
            )
        group_summary_frame = summary.loc[summary["group"] == group["name"]].sort_values(
            "voltage_mV"
        )
        sem = group_summary_frame["sem"].to_numpy(float)
        sem = np.where(np.isfinite(sem), sem, 0.0)
        axis.errorbar(
            group_summary_frame["voltage_mV"],
            group_summary_frame["mean"],
            yerr=sem,
            color=group["color"],
            marker="o",
            ms=4,
            lw=1.8,
            capsize=2.5,
            label=f"{group['name']} (n={int(group_summary_frame['n_cells'].max())})",
        )
    axis.axhline(0, color="#777777", lw=0.7)
    axis.axvline(0, color="#777777", lw=0.7)
    axis.set_xlabel("Voltage (mV)")
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    if global_result.get("status") == "performed":
        annotation = (
            f"Whole-curve {global_result['method'].replace('_', ' ')}\n"
            f"p={global_result['p_value']:.4g}"
        )
    else:
        annotation = "Descriptive curve; inferential test not performed"
    axis.text(
        0.02,
        0.98,
        annotation,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
    )
    axis.legend(frameon=False)
    style_axis(axis)
    figure.tight_layout()
    return save_figure(figure, directory, stem)


def json_clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(json_clean(value), ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )


def inventory_payload(
    selected: list[dict[str, Any]],
    inventory: list[dict[str, Any]],
    probes: list[dict[str, Any]],
    runtime: dict[str, Any],
) -> dict[str, Any]:
    counts = {
        group["name"]: sum(item["group"] == group["name"] for item in selected)
        for group in runtime["groups"]
    }
    return {
        "mode": "read_only_inventory",
        "raw_data": str(runtime["raw_root"]),
        "output_not_written": str(runtime["output_root"]),
        "group_cell_counts": counts,
        "selected_abf_count": len(selected),
        "excluded_abf_count": sum(not item["selected"] for item in inventory),
        "selected": probes,
        "all_candidates": inventory,
    }


def build_report(
    runtime: dict[str, Any],
    config_path: Path,
    inventory: list[dict[str, Any]],
    extractions: list[dict[str, Any]],
    voltages: list[float],
    global_result: dict[str, Any],
    pointwise: pd.DataFrame,
    normalized_result: dict[str, Any] | None,
    output_files: list[str],
) -> str:
    analysis = runtime["analysis"]
    extraction = runtime["extraction"]
    lines = [
        f"# {runtime['project_name']} I-V Analysis Report",
        "",
        f"Generated: {utc_now()}",
        "",
        "## Scope and provenance",
        "",
        f"- Raw ABF directory: `{runtime['raw_root']}`",
        f"- Output directory: `{runtime['output_root']}`",
        f"- Configuration: `{config_path}`",
        f"- Selected I-V ABFs/cells: {len(extractions)}",
        f"- Excluded ABFs by selection rules: {sum(not item['selected'] for item in inventory)}",
        "- Raw ABFs were opened read-only and verified by SHA-256 before and after analysis.",
        "",
        "## Dataset",
        "",
    ]
    for group in runtime["groups"]:
        count = sum(item["group"] == group["name"] for item in extractions)
        lines.append(
            f"- {group['name']} ({group['role']}, {group['color']}): n={count} cells"
        )
    lines.extend(
        [
            f"- Voltage set: {voltages[0]:g} to {voltages[-1]:g} mV "
            f"({len(voltages)} levels)",
            "- Statistical unit: cell; repeated same-voltage sweeps were averaged within cell.",
            "",
            "## Extraction",
            "",
            f"- Current measure: `{analysis['current_measure']}`",
            f"- Baseline window: {extraction['baseline_window_ms']:g} ms immediately before the step",
            f"- Steady-state window: final {extraction['steady_state_window_ms']:g} ms of the step",
            f"- Peak artifact guard: {extraction['artifact_guard_ms']:g} ms at each step edge",
            "- Trace panels are group means calculated from every included cell.",
            "",
            "## Results",
            "",
        ]
    )
    if global_result.get("status") == "performed":
        lines.append(
            f"- Whole-curve {global_result['method']}: "
            f"p={global_result['p_value']:.6g}; "
            f"evaluated={global_result['evaluated_labelings_or_permutations']}."
        )
    else:
        lines.append(f"- Whole-curve inference skipped: {global_result.get('reason')}.")
    if analysis["pointwise_tests"] and not pointwise.empty:
        correction = analysis["multiple_comparison_correction"]
        lines.append(
            f"- Supplementary pointwise Welch tests: {len(pointwise)} voltages; "
            f"raw p values retained; multiple-comparison correction=`{correction}`."
        )
    else:
        lines.append("- Supplementary pointwise tests were not produced.")
    if normalized_result is None:
        lines.append("- Per-cell normalized I-V analysis: disabled.")
    elif normalized_result.get("status") == "performed":
        lines.append(
            "- Exploratory per-cell normalized I-V: enabled; whole-curve "
            f"p={normalized_result['p_value']:.6g}. This does not replace the primary curve."
        )
    else:
        lines.append(
            "- Exploratory per-cell normalized I-V: enabled; inference skipped: "
            f"{normalized_result.get('reason')}."
        )
    lines.extend(
        [
            "",
            "## Deliverables",
            "",
        ]
    )
    lines.extend(f"- `{item}`" for item in output_files)
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "These outputs describe currents recorded under the configured extraction rule. "
            "They do not establish channel identity, mechanism, animal-level replication, "
            "or causal biology. Python-generated statistics are authoritative; spreadsheet "
            "and figure-data exports are derived plotting resources.",
            "",
        ]
    )
    return "\n".join(lines)


def run(config_path: Path, inventory_only: bool, resume: bool) -> int:
    resolved, runtime = load_and_validate_config(config_path)
    selected, inventory = discover_files(runtime)
    source_hashes_before = {item["relative_path"]: sha256(item["path"]) for item in selected}
    extractions = [
        extract_recording(item, runtime["extraction"], source_hashes_before[item["relative_path"]])
        for item in selected
    ]
    voltages = validate_voltage_sets(extractions)

    if inventory_only:
        probes = [
            {
                "sample": item["sample"],
                "group": item["group"],
                "relative_path": item["relative_path"],
                "protocol": item["protocol"],
                "sweep_count": item["sweep_count"],
                "sampling_rate_Hz": item["sampling_rate_Hz"],
                "voltage_levels_mV": sorted(set(item["voltages_mV"])),
                "step_start_ms": item["step_start_ms"],
                "step_end_ms": item["step_end_ms"],
                "source_sha256": item["source_sha256"],
            }
            for item in extractions
        ]
        print(
            json.dumps(
                inventory_payload(selected, inventory, probes, runtime),
                ensure_ascii=False,
                indent=2,
            )
        )
        return 0

    directories = prepare_output(runtime["output_root"], resume)
    resolved_config_path = directories["root"] / "resolved_config.yaml"
    resolved_config_path.write_text(
        yaml.safe_dump(resolved, sort_keys=False, allow_unicode=True), encoding="utf-8"
    )

    capacitance = load_capacitance(runtime, [item["sample"] for item in extractions])
    expected_workbooks: dict[str, Path] = {}
    used_names: set[str] = set()
    for item in extractions:
        filename = f"{safe_filename(item['sample'])}_IV.xlsx"
        if filename.lower() in used_names:
            raise ValueError(f"Sample filenames collide after sanitization: {filename}")
        used_names.add(filename.lower())
        expected_workbooks[item["sample"]] = directories["per_cell"] / filename
    stale = {
        path.name for path in directories["per_cell"].glob("*_IV.xlsx")
    } - {path.name for path in expected_workbooks.values()}
    if stale:
        raise ValueError(
            "--resume found stale per-cell workbooks; use a new output directory: "
            + ", ".join(sorted(stale))
        )

    for item in extractions:
        write_per_cell_workbook(
            item,
            expected_workbooks[item["sample"]],
            capacitance.get(item["sample"]),
        )

    group_order = [item["name"] for item in runtime["groups"]]
    current_measure = runtime["analysis"]["current_measure"]
    endpoint_unit = CURRENT_MEASURES[current_measure][1]
    cell_data = build_cell_data(extractions, capacitance, current_measure)
    summary = group_summary(cell_data, group_order, "endpoint_value", endpoint_unit)
    global_result, pointwise = run_inference(
        cell_data, "endpoint_value", group_order, runtime["analysis"]
    )
    global_result["current_measure"] = current_measure
    global_result["endpoint_unit"] = endpoint_unit

    prefix = runtime["prefix"]
    cell_csv = directories["extracted"] / f"{prefix}_IV_Cell_Data.csv"
    summary_csv = directories["extracted"] / f"{prefix}_IV_Group_Summary.csv"
    cell_data.to_csv(cell_csv, index=False, encoding="utf-8-sig")
    summary.to_csv(summary_csv, index=False, encoding="utf-8-sig")

    registry_rows: list[dict[str, Any]] = []
    for item in extractions:
        registry_rows.append(
            {
                "sample": item["sample"],
                "group": item["group"],
                "role": item["role"],
                "source_abf": item["relative_path"],
                "source_sha256": item["source_sha256"],
                "protocol": item["protocol"],
                "sampling_rate_Hz": item["sampling_rate_Hz"],
                "sweep_count": item["sweep_count"],
                "voltage_levels_mV": sorted(set(item["voltages_mV"])),
                "per_cell_workbook": expected_workbooks[item["sample"]].relative_to(
                    directories["root"]
                ).as_posix(),
                "capacitance_pF": capacitance.get(item["sample"]),
            }
        )
    registry_path = directories["extracted"] / "sample_registry.json"
    write_json(
        registry_path,
        {
            "statistical_unit": "cell",
            "current_measure": current_measure,
            "samples": registry_rows,
        },
    )

    trace_table = build_trace_table(extractions, group_order, voltages)
    trace_paths = plot_group_mean_traces(
        trace_table, runtime["groups"], voltages, directories["figures"], prefix
    )
    curve_paths = plot_iv_curve(
        cell_data,
        summary,
        runtime["groups"],
        "endpoint_value",
        f"{current_measure.replace('_', ' ')} ({endpoint_unit})",
        "Steady-state I-V relationship",
        global_result,
        directories["figures"],
        f"{prefix}_IV_Curve",
    )

    normalized_data = pd.DataFrame()
    normalized_summary = pd.DataFrame()
    normalized_pointwise = pd.DataFrame()
    normalized_result: dict[str, Any] | None = None
    normalized_paths: list[str] = []
    normalization = runtime["analysis"]["normalization"]
    if normalization["enabled"]:
        normalized_data = normalize_cell_data(
            cell_data,
            normalization["reference_voltage_mV"],
            normalization["minimum_denominator_abs"],
        )
        normalized_summary = group_summary(
            normalized_data, group_order, "normalized_value", "ratio"
        )
        normalized_result, normalized_pointwise = run_inference(
            normalized_data,
            "normalized_value",
            group_order,
            runtime["analysis"],
        )
        normalized_result["analysis_role"] = "exploratory"
        normalized_result["reference_voltage_mV"] = normalization[
            "reference_voltage_mV"
        ]
        normalized_paths = plot_iv_curve(
            normalized_data,
            normalized_summary,
            runtime["groups"],
            "normalized_value",
            "Normalized current (ratio)",
            "Exploratory normalized I-V relationship",
            normalized_result,
            directories["figures"],
            f"{prefix}_Normalized_IV_Curve",
        )

    figure_data_path = directories["figures"] / f"{prefix}_IV_Figure_Data.xlsx"
    with pd.ExcelWriter(figure_data_path, engine="openpyxl") as writer:
        cell_data.to_excel(writer, sheet_name="Main_Cell_Data", index=False)
        summary.to_excel(writer, sheet_name="Main_Group_Summary", index=False)
        trace_table.to_excel(writer, sheet_name="Trace_GroupMean", index=False)
        if not normalized_data.empty:
            normalized_data.to_excel(writer, sheet_name="Normalized_Cell_Data", index=False)
            normalized_summary.to_excel(
                writer, sheet_name="Normalized_Summary", index=False
            )

    statistics_path = directories["statistics"] / f"{prefix}_IV_Statistics.xlsx"
    with pd.ExcelWriter(statistics_path, engine="openpyxl") as writer:
        cell_data.to_excel(writer, sheet_name="Cell_Data", index=False)
        summary.to_excel(writer, sheet_name="Group_Summary", index=False)
        pd.DataFrame([json_clean(global_result)]).to_excel(
            writer, sheet_name="Global_Test", index=False
        )
        pointwise.to_excel(writer, sheet_name="Pointwise_Tests", index=False)
        if not normalized_data.empty:
            normalized_data.to_excel(writer, sheet_name="Norm_Cell_Data", index=False)
            normalized_summary.to_excel(writer, sheet_name="Norm_Summary", index=False)
            pd.DataFrame([json_clean(normalized_result)]).to_excel(
                writer, sheet_name="Norm_Global_Test", index=False
            )
            normalized_pointwise.to_excel(
                writer, sheet_name="Norm_Pointwise", index=False
            )

    results = {
        "script_version": SCRIPT_VERSION,
        "statistical_unit": "cell",
        "group_order": group_order,
        "group_cell_counts": {
            group: int(cell_data.loc[cell_data["group"] == group, "sample"].nunique())
            for group in group_order
        },
        "voltage_levels_mV": voltages,
        "current_measure": current_measure,
        "endpoint_unit": endpoint_unit,
        "primary_global_test": global_result,
        "pointwise_test_count": int(len(pointwise)),
        "pointwise_multiple_comparison_correction": runtime["analysis"][
            "multiple_comparison_correction"
        ],
        "normalization": {
            "enabled": normalization["enabled"],
            "analysis_role": "exploratory" if normalization["enabled"] else None,
            "global_test": normalized_result,
        },
        "trace_panels": "group means calculated from every included cell",
    }
    results_path = directories["statistics"] / f"{prefix}_IV_Results.json"
    write_json(results_path, results)

    report_path = directories["report"] / f"{prefix}_IV_Report.md"
    deliverable_paths = [
        cell_csv,
        summary_csv,
        registry_path,
        figure_data_path,
        statistics_path,
        results_path,
        *[Path(item) for item in trace_paths + curve_paths + normalized_paths],
    ]
    report_text = build_report(
        runtime,
        config_path,
        inventory,
        extractions,
        voltages,
        global_result,
        pointwise,
        normalized_result,
        [path.relative_to(directories["root"]).as_posix() for path in deliverable_paths],
    )
    report_path.write_text(report_text, encoding="utf-8")

    source_hashes_after = {item["relative_path"]: sha256(item["path"]) for item in selected}
    if source_hashes_after != source_hashes_before:
        raise RuntimeError("A source ABF hash changed during analysis")

    output_records: list[dict[str, Any]] = []
    for path in sorted(
        (item for item in directories["root"].rglob("*") if item.is_file()),
        key=lambda item: item.as_posix().lower(),
    ):
        if path.name == "analysis_manifest.json":
            continue
        output_records.append(
            {
                "path": path.relative_to(directories["root"]).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    manifest = {
        "generated_at_utc": utc_now(),
        "script_version": SCRIPT_VERSION,
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "raw_data": str(runtime["raw_root"]),
        "output": str(runtime["output_root"]),
        "inputs": [
            {"path": key, "sha256_before": value, "sha256_after": source_hashes_after[key]}
            for key, value in sorted(source_hashes_before.items())
        ],
        "outputs": output_records,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "pyabf": getattr(pyabf, "__version__", "unknown"),
            "pyyaml": getattr(yaml, "__version__", "unknown"),
        },
    }
    manifest_path = directories["root"] / "analysis_manifest.json"
    write_json(manifest_path, manifest)

    print(f"Completed I-V analysis: {runtime['output_root']}")
    print(
        "Cells by group: "
        + ", ".join(
            f"{group}={results['group_cell_counts'][group]}" for group in group_order
        )
    )
    if global_result.get("status") == "performed":
        print(
            f"Whole-curve {global_result['method']}: p={global_result['p_value']:.6g}"
        )
    else:
        print(f"Whole-curve inference skipped: {global_result.get('reason')}")
    print(f"Manifest: {manifest_path}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Analyze cell-level I-V curves from ABF voltage-step recordings"
    )
    parser.add_argument("--config", required=True, help="Path to analysis YAML")
    parser.add_argument(
        "--inventory-only",
        action="store_true",
        help="Read and validate ABFs, print the selection/protocol inventory, and write nothing",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Regenerate known outputs in the same configured output directory without deleting files",
    )
    args = parser.parse_args()
    if args.inventory_only and args.resume:
        parser.error("--inventory-only and --resume cannot be combined")
    config_path = Path(args.config).expanduser().resolve()
    if not config_path.is_file():
        parser.error(f"Config file does not exist: {config_path}")
    try:
        return run(config_path, args.inventory_only, args.resume)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    raise SystemExit(main())
