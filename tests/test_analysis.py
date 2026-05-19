import pandas as pd
import sqlite3
import sys
import os
import plotly.graph_objects as go

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from analysis import (
    get_all_samples, compute_frequencies,
    get_melanoma_miraclib, run_mannwhitney, create_boxplot
)

conn = sqlite3.connect("clinical_trial.db")

print("=" * 50)
print("TESTING analysis.py - Part 2")
print("=" * 50)

print("\nTest 1: get_all_samples returns correct row count")
df_samples = get_all_samples(conn)
assert len(df_samples) == 10500, f"Expected 10500, got {len(df_samples)}"
print("    PASS")

print("\nTest 2: get_all_samples has correct columns")
expected_cols = {"sample_id", "b_cell", "cd8_t_cell", "cd4_t_cell", 
                 "nk_cell", "monocyte", "condition", "project_id"}
assert expected_cols.issubset(set(df_samples.columns)), \
    f"Missing columns: {expected_cols - set(df_samples.columns)}"
print("    PASS")

print("\nTest 3: compute_frequencies returns 52500 rows")
df_summary = compute_frequencies(df_samples)
assert len(df_summary) == 52500, f"Expected 52500, got {len(df_summary)}"
print("    PASS")

print("\nTest 4: correct columns in output")
expected = ["sample", "total_count", "population", "count", "percentage"]
assert list(df_summary.columns) == expected, \
    f"Column mismatch: {list(df_summary.columns)}"
print("    PASS")

print("\nTest 5: exactly 5 populations per sample")
pop_counts = df_summary.groupby("sample")["population"].count()
assert (pop_counts == 5).all(), "Some samples don't have exactly 5 populations"
print("    PASS")

print("\nTest 6: percentages between 0 and 100")
assert df_summary["percentage"].between(0, 100).all(), \
    "Some percentages are out of range"
print("    PASS")

print("\nTest 7: project_id filter works")
df_prj1 = compute_frequencies(df_samples, project_id="prj1")
assert df_prj1["sample"].nunique() == \
    df_samples[df_samples["project_id"] == "prj1"]["sample_id"].nunique()
print("    PASS")

print("\nTest 8: condition filter works")
df_melanoma = compute_frequencies(df_samples, condition="melanoma")
assert df_melanoma["sample"].nunique() == \
    df_samples[df_samples["condition"] == "melanoma"]["sample_id"].nunique()
print("    PASS")

print("\nTest 9: verify frequencies for sample00000")
s = df_summary[df_summary["sample"] == "sample00000"]
assert len(s) == 5, "sample00000 should have 5 rows"
bcell = s[s["population"] == "b_cell"].iloc[0]
assert bcell["count"] == 10908, f"Expected 10908, got {bcell['count']}"
assert bcell["total_count"] == 93214, f"Expected 93214, got {bcell['total_count']}"
print("    PASS")

print("\n" + "=" * 50)
print("TESTING analysis.py - Part 3")
print("=" * 50)

df_mel = get_melanoma_miraclib(conn)

print("\nTest 10: get_melanoma_miraclib returns correct row count")
assert len(df_mel) == 1968, f"Expected 1968, got {len(df_mel)}"
print("    PASS")

print("\nTest 11: get_melanoma_miraclib only contains melanoma condition")
assert (df_mel["condition"] == "melanoma").all(), "Non-melanoma rows found"
print("    PASS")

print("\nTest 12: get_melanoma_miraclib only contains miraclib treatment")
assert (df_mel["treatment"] == "miraclib").all(), "Non-miraclib rows found"
print("    PASS")

print("\nTest 13: get_melanoma_miraclib only contains PBMC samples")
assert (df_mel["sample_type"] == "PBMC").all(), "Non-PBMC rows found"
print("    PASS")

print("\nTest 14: get_melanoma_miraclib only contains yes/no response values")
assert set(df_mel["response"].unique()) == {"yes", "no"}, \
    f"Unexpected response values: {df_mel['response'].unique()}"
print("    PASS")

df_freq = compute_frequencies(df_mel, extra_id_vars=["response"])

print("\nTest 15: compute_frequencies with extra_id_vars preserves response column")
assert "response" in df_freq.columns, "response column missing from output"
print("    PASS")

print("\nTest 16: compute_frequencies with extra_id_vars returns correct row count")
assert len(df_freq) == 9840, f"Expected 9840, got {len(df_freq)}"
print("    PASS")

print("\nTest 17: compute_frequencies with extra_id_vars has correct columns")
expected = ["sample", "total_count", "response", "population", "count", "percentage"]
assert list(df_freq.columns) == expected, f"Column mismatch: {list(df_freq.columns)}"
print("    PASS")

df_stats = run_mannwhitney(df_freq)

print("\nTest 18: run_mannwhitney returns one row per cell population")
assert len(df_stats) == 5, f"Expected 5, got {len(df_stats)}"
print("    PASS")

print("\nTest 19: run_mannwhitney has correct columns")
expected = ["population", "u_statistic", "p_value", "significant"]
assert list(df_stats.columns) == expected, f"Column mismatch: {list(df_stats.columns)}"
print("    PASS")

print("\nTest 20: run_mannwhitney p_values are between 0 and 1")
assert df_stats["p_value"].between(0, 1).all(), "p_values out of range"
print("    PASS")

print("\nTest 21: run_mannwhitney significant column only contains yes/no")
assert set(df_stats["significant"].unique()).issubset({"yes", "no"}), \
    f"Unexpected significant values: {df_stats['significant'].unique()}"
print("    PASS")

print("\nTest 22: run_mannwhitney results sorted by p_value ascending")
assert df_stats["p_value"].is_monotonic_increasing, "Results not sorted by p_value"
print("    PASS")

fig = create_boxplot(df_freq)

print("\nTest 23: create_boxplot returns a plotly Figure")
assert isinstance(fig, go.Figure), f"Expected go.Figure, got {type(fig)}"
print("    PASS")

print("\nTest 24: create_boxplot has two traces (responders and non-responders)")
assert len(fig.data) == 2, f"Expected 2 traces, got {len(fig.data)}"
print("    PASS")

print("\n" + "=" * 50)
print("ALL TESTS PASSED")
print("=" * 50)

conn.close()