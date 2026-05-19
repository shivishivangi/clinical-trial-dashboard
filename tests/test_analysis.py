import pandas as pd
import sqlite3
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from analysis import get_all_samples, compute_frequencies

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
print("ALL TESTS PASSED")
print("=" * 50)

conn.close()