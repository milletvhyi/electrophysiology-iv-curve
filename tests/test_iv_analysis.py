import importlib.util
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_iv_analysis.py"
SPEC = importlib.util.spec_from_file_location("run_iv_analysis", MODULE_PATH)
iv = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(iv)


class StepDetectionTests(unittest.TestCase):
    def settings(self):
        return {
            "command_threshold_mV": 0.5,
            "min_step_duration_ms": 20.0,
            "step_alignment_tolerance_ms": 1.0,
            "baseline_window_ms": 10.0,
            "steady_state_window_ms": 20.0,
            "artifact_guard_ms": 2.0,
        }

    def test_consensus_window_handles_step_equal_to_holding(self):
        command = np.full((3, 300), -70.0)
        command[0, 100:250] = -100.0
        command[1, 100:250] = -70.0
        command[2, 100:250] = 100.0

        result = iv.detect_step_window(command, 1000, self.settings())

        self.assertEqual(result["start_index"], 100)
        self.assertEqual(result["end_index"], 250)
        self.assertIsNone(result["detected_segments"][1])
        self.assertEqual(result["steady_samples"], 20)

    def test_misaligned_steps_fail_closed(self):
        command = np.full((2, 300), -70.0)
        command[0, 100:250] = -100.0
        command[1, 110:260] = 100.0

        with self.assertRaisesRegex(ValueError, "not aligned"):
            iv.detect_step_window(command, 1000, self.settings())


class StatisticsTests(unittest.TestCase):
    def test_exact_curve_permutation_preserves_cell_profiles(self):
        matrix = np.array(
            [
                [0.0, 0.0],
                [0.0, 0.0],
                [10.0, 10.0],
                [10.0, 10.0],
            ]
        )
        labels = np.array(["A", "A", "B", "B"])

        result = iv.curve_permutation_test(
            matrix,
            labels,
            "A",
            "B",
            max_exact_labelings=100,
            monte_carlo_permutations=50,
            random_seed=7,
        )

        self.assertEqual(result["method"], "exact_label_permutation")
        self.assertEqual(result["evaluated_labelings_or_permutations"], 6)
        self.assertAlmostEqual(result["p_value"], 1 / 3)

    def test_holm_adjustment_is_monotone_in_rank_order(self):
        adjusted = iv.holm_adjust([0.01, 0.04, 0.03])
        np.testing.assert_allclose(adjusted, [0.03, 0.06, 0.06])


class CellAggregationTests(unittest.TestCase):
    def test_repeated_sweeps_are_averaged_within_cell(self):
        params = pd.DataFrame(
            {
                "voltage_mV": [0.0, 0.0, 100.0],
                "baseline_current_pA": [1.0, 3.0, 5.0],
                "steady_state_raw_pA": [10.0, 14.0, 30.0],
                "steady_state_baseline_subtracted_pA": [9.0, 11.0, 25.0],
                "peak_raw_pA": [20.0, 24.0, 40.0],
                "peak_baseline_subtracted_pA": [19.0, 21.0, 35.0],
            }
        )
        extractions = [
            {
                "sample": "cell1",
                "group": "Control",
                "role": "control",
                "color": "#111111",
                "params": params,
            }
        ]

        result = iv.build_cell_data(extractions, {}, "raw_steady_state")
        zero = result.loc[np.isclose(result["voltage_mV"], 0)].iloc[0]

        self.assertEqual(zero["sweep_count_at_voltage"], 2)
        self.assertAlmostEqual(zero["endpoint_value"], 12.0)

    def test_normalization_is_per_cell_and_guarded(self):
        data = pd.DataFrame(
            {
                "sample": ["c1", "c1", "c2", "c2"],
                "group": ["A", "A", "B", "B"],
                "role": ["control", "control", "experimental", "experimental"],
                "color": ["#111111", "#111111", "#E52521", "#E52521"],
                "voltage_mV": [0.0, 100.0, 0.0, 100.0],
                "endpoint_value": [5.0, 10.0, 4.0, 8.0],
            }
        )

        normalized = iv.normalize_cell_data(data, 100.0, 1.0)
        np.testing.assert_allclose(normalized["normalized_value"], [0.5, 1.0, 0.5, 1.0])

        data.loc[data["sample"] == "c2", "endpoint_value"] = [4.0, 0.1]
        with self.assertRaisesRegex(ValueError, "denominators"):
            iv.normalize_cell_data(data, 100.0, 1.0)


class RoutingTests(unittest.TestCase):
    def test_overlapping_group_prefixes_are_rejected(self):
        groups = [
            {"folder_prefix": "c", "name": "A"},
            {"folder_prefix": "con", "name": "B"},
        ]
        with self.assertRaisesRegex(ValueError, "matched 2"):
            iv.group_for_sample("con1", groups)


if __name__ == "__main__":
    unittest.main()
