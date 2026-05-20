import unittest
import sqlite3
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


class TestLoadData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.conn = sqlite3.connect("clinical_trial.db")
        cls.df_raw = pd.read_csv("data/cell-count.csv")

    @classmethod
    def tearDownClass(cls):
        cls.conn.close()

    def test_all_tables_exist(self):
        tables = pd.read_sql("SELECT name FROM sqlite_master WHERE type='table'", self.conn)
        self.assertSetEqual({"projects", "subjects", "samples"}, set(tables["name"].tolist()))

    def test_sample_row_count_matches_csv(self):
        count = pd.read_sql("SELECT COUNT(*) as n FROM samples", self.conn).iloc[0]["n"]
        self.assertEqual(count, len(self.df_raw))

    def test_project_count_matches_csv(self):
        count = pd.read_sql("SELECT COUNT(*) as n FROM projects", self.conn).iloc[0]["n"]
        self.assertEqual(count, self.df_raw["project"].nunique())

    def test_subject_count_matches_csv(self):
        count = pd.read_sql("SELECT COUNT(*) as n FROM subjects", self.conn).iloc[0]["n"]
        self.assertEqual(count, self.df_raw["subject"].nunique())

    def test_no_null_sample_ids(self):
        nulls = pd.read_sql(
            "SELECT COUNT(*) as n FROM samples WHERE sample_id IS NULL", self.conn
        ).iloc[0]["n"]
        self.assertEqual(nulls, 0)

    def test_no_orphaned_samples(self):
        orphans = pd.read_sql("""
            SELECT COUNT(*) as n FROM samples
            WHERE subject_id NOT IN (SELECT subject_id FROM subjects)
        """, self.conn).iloc[0]["n"]
        self.assertEqual(orphans, 0)

    def test_cell_counts_non_negative(self):
        negatives = pd.read_sql("""
            SELECT COUNT(*) as n FROM samples
            WHERE b_cell < 0 OR cd8_t_cell < 0 OR cd4_t_cell < 0
            OR nk_cell < 0 OR monocyte < 0
        """, self.conn).iloc[0]["n"]
        self.assertEqual(negatives, 0)

    def test_timepoints_are_valid(self):
        timepoints = pd.read_sql(
            "SELECT DISTINCT time_from_treatment_start FROM samples", self.conn
        )
        actual = set(timepoints["time_from_treatment_start"].tolist())
        self.assertTrue(actual.issubset({0, 7, 14}))

    def test_specific_row_values(self):
        result = pd.read_sql("""
            SELECT s.sample_id, s.b_cell, s.cd8_t_cell, s.cd4_t_cell,
                   s.nk_cell, s.monocyte, s.time_from_treatment_start,
                   sub.condition, sub.sex, sub.treatment, sub.response,
                   p.project_id
            FROM samples s
            JOIN subjects sub ON s.subject_id = sub.subject_id
            JOIN projects p ON sub.project_id = p.project_id
            WHERE s.sample_id = 'sample00000'
        """, self.conn)
        self.assertEqual(len(result), 1)
        row = result.iloc[0]
        self.assertEqual(row["project_id"], "prj1")
        self.assertEqual(row["condition"], "melanoma")
        self.assertEqual(row["sex"], "M")
        self.assertEqual(row["treatment"], "miraclib")
        self.assertEqual(row["response"], "no")
        self.assertEqual(row["b_cell"], 10908)
        self.assertEqual(row["cd8_t_cell"], 24440)
        self.assertEqual(row["cd4_t_cell"], 20491)
        self.assertEqual(row["nk_cell"], 13864)
        self.assertEqual(row["monocyte"], 23511)
        self.assertEqual(row["time_from_treatment_start"], 0)

    def test_csv_crossvalidation_cell_counts(self):
        for _, csv_row in self.df_raw.head(10).iterrows():
            db_row = pd.read_sql(
                "SELECT b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte FROM samples WHERE sample_id = ?",
                self.conn, params=[csv_row["sample"]]
            ).iloc[0]
            self.assertEqual(db_row["b_cell"], csv_row["b_cell"])
            self.assertEqual(db_row["cd8_t_cell"], csv_row["cd8_t_cell"])
            self.assertEqual(db_row["cd4_t_cell"], csv_row["cd4_t_cell"])
            self.assertEqual(db_row["nk_cell"], csv_row["nk_cell"])
            self.assertEqual(db_row["monocyte"], csv_row["monocyte"])

    def test_csv_crossvalidation_conditions(self):
        db_conditions = set(
            pd.read_sql("SELECT DISTINCT condition FROM subjects", self.conn)["condition"]
        )
        csv_conditions = set(self.df_raw["condition"].unique())
        self.assertSetEqual(db_conditions, csv_conditions)


if __name__ == "__main__":
    unittest.main()
