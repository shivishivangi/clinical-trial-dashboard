import pandas as pd
import sqlite3
import os

os.makedirs("outputs", exist_ok=True)

DB_PATH = "clinical_trial.db"
CELL_COLS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def get_all_samples(conn):
    """
    Fetch all samples with cell counts, condition, and project
    from the database via a 3-table join.
    """
    return pd.read_sql("""
        SELECT 
            s.sample_id,
            s.b_cell, s.cd8_t_cell, s.cd4_t_cell, s.nk_cell, s.monocyte,
            sub.condition,
            p.project_id
        FROM samples s
        JOIN subjects sub ON s.subject_id = sub.subject_id
        JOIN projects p ON sub.project_id = p.project_id
    """, conn)


def compute_frequencies(df, project_id=None, condition=None):
    """
    Compute relative frequency of each cell population per sample.
    Optional filters: project_id, condition.
    Returns long-format dataframe with columns:
    sample, total_count, population, count, percentage
    """
    df = df.copy()

    if project_id:
        df = df[df["project_id"] == project_id]
    if condition:
        df = df[df["condition"] == condition]

    df["total_count"] = df[CELL_COLS].sum(axis=1)
    if (df["total_count"] == 0).any():
        print("WARNING: samples with zero total count detected")

    df_long = df.melt(
        id_vars=["sample_id", "total_count"],
        value_vars=CELL_COLS,
        var_name="population",
        value_name="count"
    )

    df_long["percentage"] = (
        df_long["count"] / df_long["total_count"] * 100
    ).round(4)

    df_long = df_long.rename(columns={"sample_id": "sample"})
    df_long = df_long[["sample", "total_count", "population", "count", "percentage"]]
    df_long = df_long.sort_values(["sample", "population"]).reset_index(drop=True)

    return df_long


def main():
    conn = sqlite3.connect(DB_PATH)

    print("Running Part 2: Frequency Summary Table")
    df_samples = get_all_samples(conn)
    df_summary = compute_frequencies(df_samples)
    df_summary.to_csv("outputs/summary_table.csv", index=False)
    print(f"  Successfully saved {len(df_summary)} rows to outputs/summary_table.csv")

    conn.close()


if __name__ == "__main__":
    main()