import sqlite3
import pandas as pd

conn = sqlite3.connect("clinical_trial.db")
df_raw = pd.read_csv("data/cell-count.csv")

print("=" * 50)
print("TESTING load_data.py")
print("=" * 50)

print("\nTest 1: All 3 tables exist")
tables = pd.read_sql("SELECT name FROM sqlite_master WHERE type='table'", conn)
expected = {"projects", "subjects", "samples"}
actual = set(tables["name"].tolist())
assert expected == actual, f"Missing tables: {expected - actual}"
print("    PASS")

print("\nTest 2: Row counts match original CSV")
sample_count = pd.read_sql("SELECT COUNT(*) as n FROM samples", conn).iloc[0]["n"]
assert sample_count == len(df_raw), f"Expected {len(df_raw)} samples, got {sample_count}"
print(f"    PASS ({sample_count} samples)")

print("\nTest 3: Correct number of unique projects")
project_count = pd.read_sql("SELECT COUNT(*) as n FROM projects", conn).iloc[0]["n"]
expected_projects = df_raw["project"].nunique()
assert project_count == expected_projects, f"Expected {expected_projects}, got {project_count}"
print(f"    PASS ({project_count} projects)")

print("\nTest 4: Correct number of unique subjects")
subject_count = pd.read_sql("SELECT COUNT(*) as n FROM subjects", conn).iloc[0]["n"]
expected_subjects = df_raw["subject"].nunique()
assert subject_count == expected_subjects, f"Expected {expected_subjects}, got {subject_count}"
print(f"    PASS ({subject_count} subjects)")

print("\nTest 5: No null sample IDs")
nulls = pd.read_sql("SELECT COUNT(*) as n FROM samples WHERE sample_id IS NULL", conn).iloc[0]["n"]
assert nulls == 0, f"Found {nulls} null sample IDs"
print("    PASS")

print("\nTest 6: Foreign keys are intact (check for dangling samples)")
orphans = pd.read_sql("""
    SELECT COUNT(*) as n FROM samples 
    WHERE subject_id NOT IN (SELECT subject_id FROM subjects)
""", conn).iloc[0]["n"]
assert orphans == 0, f"Found {orphans} orphaned samples"
print("    PASS")

print("\nTest 7: Cell counts are all positive")
negatives = pd.read_sql("""
    SELECT COUNT(*) as n FROM samples
    WHERE b_cell < 0 OR cd8_t_cell < 0 OR cd4_t_cell < 0 
    OR nk_cell < 0 OR monocyte < 0
""", conn).iloc[0]["n"]
assert negatives == 0, f"Found {negatives} negative cell counts"
print("    PASS")

print("\nTest 8: time_from_treatment_start only has expected values")
timepoints = pd.read_sql("SELECT DISTINCT time_from_treatment_start FROM samples", conn)
valid = {0, 7, 14}
actual_timepoints = set(timepoints["time_from_treatment_start"].tolist())
assert actual_timepoints.issubset(valid), f"Unexpected timepoints: {actual_timepoints - valid}"
print(f"    PASS (timepoints: {sorted(actual_timepoints)})")

print("\nTest 9: Verify specific row from CSV exists in DB")
result = pd.read_sql("""
    SELECT s.sample_id, s.b_cell, s.cd8_t_cell, s.cd4_t_cell, 
           s.nk_cell, s.monocyte, s.time_from_treatment_start,
           sub.condition, sub.sex, sub.treatment, sub.response,
           p.project_id
    FROM samples s
    JOIN subjects sub ON s.subject_id = sub.subject_id
    JOIN projects p ON sub.project_id = p.project_id
    WHERE s.sample_id = 'sample00000'
""", conn)

assert len(result) == 1, "sample00000 not found in DB"
row = result.iloc[0]
assert row["project_id"] == "prj1"
assert row["condition"] == "melanoma"
assert row["sex"] == "M"
assert row["treatment"] == "miraclib"
assert row["response"] == "no"
assert row["b_cell"] == 10908
assert row["cd8_t_cell"] == 24440
assert row["cd4_t_cell"] == 20491
assert row["nk_cell"] == 13864
assert row["monocyte"] == 23511
assert row["time_from_treatment_start"] == 0
print("    PASS (sample00000 verified)")

print("\n" + "=" * 50)
print("ALL TESTS PASSED")
print("=" * 50)

conn.close()