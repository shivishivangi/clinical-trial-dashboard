import unittest
import sqlite3
import pandas as pd
import sys
import os
import plotly.graph_objects as go

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from analysis import (
    get_all_samples, compute_frequencies,
    get_melanoma_miraclib, run_mannwhitney, create_boxplot,
    get_baseline_melanoma, compute_subset_summary, compute_avg_cell,
)


class TestFrequencySummary(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.conn = sqlite3.connect("clinical_trial.db")
        cls.df_samples = get_all_samples(cls.conn)
        cls.df_summary = compute_frequencies(cls.df_samples)

    @classmethod
    def tearDownClass(cls):
        cls.conn.close()

    def test_row_count(self):
        self.assertEqual(len(self.df_samples), 10500)

    def test_columns(self):
        expected = {"sample_id", "b_cell", "cd8_t_cell", "cd4_t_cell",
                    "nk_cell", "monocyte", "condition", "project_id"}
        self.assertTrue(expected.issubset(set(self.df_samples.columns)))

    def test_frequency_row_count(self):
        self.assertEqual(len(self.df_summary), 52500)

    def test_frequency_columns(self):
        self.assertEqual(
            list(self.df_summary.columns),
            ["sample", "total_count", "population", "count", "percentage"]
        )

    def test_five_populations_per_sample(self):
        pop_counts = self.df_summary.groupby("sample")["population"].count()
        self.assertTrue((pop_counts == 5).all())

    def test_percentages_in_range(self):
        self.assertTrue(self.df_summary["percentage"].between(0, 100).all())

    def test_project_filter(self):
        df_prj1 = compute_frequencies(self.df_samples, project_id="prj1")
        expected = self.df_samples[self.df_samples["project_id"] == "prj1"]["sample_id"].nunique()
        self.assertEqual(df_prj1["sample"].nunique(), expected)

    def test_condition_filter(self):
        df_mel = compute_frequencies(self.df_samples, condition="melanoma")
        expected = self.df_samples[self.df_samples["condition"] == "melanoma"]["sample_id"].nunique()
        self.assertEqual(df_mel["sample"].nunique(), expected)

    def test_spot_check_sample00000(self):
        s = self.df_summary[self.df_summary["sample"] == "sample00000"]
        self.assertEqual(len(s), 5)
        bcell = s[s["population"] == "b_cell"].iloc[0]
        self.assertEqual(bcell["count"], 10908)
        self.assertEqual(bcell["total_count"], 93214)


class TestStatisticalAnalysis(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.conn = sqlite3.connect("clinical_trial.db")
        cls.df_mel = get_melanoma_miraclib(cls.conn)
        cls.df_freq = compute_frequencies(cls.df_mel, extra_id_vars=["response"])
        cls.df_stats = run_mannwhitney(cls.df_freq)
        cls.fig = create_boxplot(cls.df_freq)

    @classmethod
    def tearDownClass(cls):
        cls.conn.close()

    def test_melanoma_miraclib_row_count(self):
        self.assertEqual(len(self.df_mel), 1968)

    def test_only_melanoma_condition(self):
        self.assertTrue((self.df_mel["condition"] == "melanoma").all())

    def test_only_miraclib_treatment(self):
        self.assertTrue((self.df_mel["treatment"] == "miraclib").all())

    def test_only_pbmc_samples(self):
        self.assertTrue((self.df_mel["sample_type"] == "PBMC").all())

    def test_only_yes_no_response(self):
        self.assertSetEqual(set(self.df_mel["response"].unique()), {"yes", "no"})

    def test_response_column_preserved(self):
        self.assertIn("response", self.df_freq.columns)

    def test_frequency_row_count_with_response(self):
        self.assertEqual(len(self.df_freq), 9840)

    def test_frequency_columns_with_response(self):
        self.assertEqual(
            list(self.df_freq.columns),
            ["sample", "total_count", "response", "population", "count", "percentage"]
        )

    def test_mannwhitney_row_count(self):
        self.assertEqual(len(self.df_stats), 5)

    def test_mannwhitney_columns(self):
        self.assertEqual(
            list(self.df_stats.columns),
            ["population", "u_statistic", "p_value", "significant_p05", "significant_bonferroni"]
        )

    def test_pvalues_in_range(self):
        self.assertTrue(self.df_stats["p_value"].between(0, 1).all())

    def test_significant_columns_values(self):
        self.assertTrue(set(self.df_stats["significant_p05"].unique()).issubset({"yes", "no"}))
        self.assertTrue(set(self.df_stats["significant_bonferroni"].unique()).issubset({"yes", "no"}))

    def test_results_sorted_by_pvalue(self):
        self.assertTrue(self.df_stats["p_value"].is_monotonic_increasing)

    def test_boxplot_returns_figure(self):
        self.assertIsInstance(self.fig, go.Figure)

    def test_boxplot_two_traces(self):
        self.assertEqual(len(self.fig.data), 2)


class TestSubsetAnalysis(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.conn = sqlite3.connect("clinical_trial.db")
        cls.df_baseline = get_baseline_melanoma(cls.conn)
        (cls.samples_per_project,
         cls.response_counts,
         cls.sex_counts) = compute_subset_summary(cls.df_baseline)
        cls.avg_bcell = compute_avg_cell(cls.df_baseline, cell_col="b_cell", sex="M", response="yes")

    @classmethod
    def tearDownClass(cls):
        cls.conn.close()

    def test_baseline_row_count(self):
        self.assertEqual(len(self.df_baseline), 656)

    def test_only_time_zero(self):
        self.assertTrue((self.df_baseline["time_from_treatment_start"] == 0).all())

    def test_only_melanoma_miraclib(self):
        self.assertTrue((self.df_baseline["condition"] == "melanoma").all())
        self.assertTrue((self.df_baseline["treatment"] == "miraclib").all())

    def test_samples_per_project(self):
        prj = self.samples_per_project.set_index("project_id")["sample_count"]
        self.assertEqual(prj["prj1"], 384)
        self.assertEqual(prj["prj3"], 272)

    def test_response_counts(self):
        resp = self.response_counts.set_index("response")["subject_count"]
        self.assertEqual(resp["yes"], 331)
        self.assertEqual(resp["no"], 325)

    def test_sex_counts(self):
        sex = self.sex_counts.set_index("sex")["subject_count"]
        self.assertEqual(sex["M"], 344)
        self.assertEqual(sex["F"], 312)

    def test_avg_bcell_correct_value(self):
        self.assertAlmostEqual(self.avg_bcell, 10401.28, places=2)

    def test_avg_bcell_varies_by_sex(self):
        avg_female = compute_avg_cell(self.df_baseline, cell_col="b_cell", sex="F", response="yes")
        self.assertNotEqual(avg_female, self.avg_bcell)
        self.assertIsInstance(avg_female, float)

    def test_empty_result_for_invalid_timepoint(self):
        df_empty = get_baseline_melanoma(self.conn, time_point=999)
        self.assertEqual(len(df_empty), 0)


if __name__ == "__main__":
    unittest.main()
