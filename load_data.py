import pandas as pd
import sqlite3

DB_PATH = "clinical_trial.db"

CREATE_TABLES_SQL = [
    """
    CREATE TABLE IF NOT EXISTS projects (
        project_id TEXT PRIMARY KEY
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS subjects (
        subject_id  TEXT PRIMARY KEY,
        project_id  TEXT NOT NULL,
        condition   TEXT,
        age         INTEGER,
        sex         TEXT,
        treatment   TEXT,
        response    TEXT,
        FOREIGN KEY (project_id) REFERENCES projects (project_id)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS samples (
        sample_id                   TEXT PRIMARY KEY,
        subject_id                  TEXT NOT NULL,
        sample_type                 TEXT,
        time_from_treatment_start   INTEGER,
        b_cell                      INTEGER,
        cd8_t_cell                  INTEGER,
        cd4_t_cell                  INTEGER,
        nk_cell                     INTEGER,
        monocyte                    INTEGER,
        FOREIGN KEY (subject_id) REFERENCES subjects (subject_id)
    )
    """
]

def create_tables(conn):
    for sql in CREATE_TABLES_SQL:
        conn.execute(sql)
    conn.commit()

def load_csv(file_path):
    df = pd.read_csv(file_path)
    return df

def insert_rows(conn, table, columns, rows):
    placeholders = ", ".join(["?"] * len(columns))
    column_list = ", ".join(columns)
    sql = f"INSERT OR IGNORE INTO {table} ({column_list}) VALUES ({placeholders})"
    conn.executemany(sql, rows)

def populate_tables(conn, df):

    # projects table
    # unique project IDs 
    projects = df[["project"]].drop_duplicates().rename(columns={"project": "project_id"})

    # subjects table
    subjects = df[[
        "subject", "project", "condition", "age", "sex", "treatment", "response"
    ]].drop_duplicates(subset=["subject"]).rename(columns={
        "subject": "subject_id",
        "project": "project_id"
    })

    # samples table
    # one row per sample, with subject_id as foreign key
    samples = df[[
        "sample", "subject", "sample_type", "time_from_treatment_start",
        "b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"
    ]].drop_duplicates(subset=["sample"]).rename(columns={
        "sample": "sample_id",
        "subject": "subject_id"
    })    
    
    # insert into tables
    insert_rows(conn, "projects", ["project_id"], projects.values.tolist())
    insert_rows(conn, "subjects", [
        "subject_id", "project_id", "condition", "age", "sex", "treatment", "response"
    ], subjects.values.tolist())
    insert_rows(conn, "samples", [
        "sample_id", "subject_id", "sample_type", "time_from_treatment_start",
        "b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"
    ], samples.values.tolist())

    conn.commit()

def main():
    df = load_csv("data/cell-count.csv")

    with sqlite3.connect(DB_PATH) as conn:
        create_tables(conn)
        populate_tables(conn, df)

    print("\nDatabase created and data loaded successfully.")
    print(f"  Projects: {df['project'].nunique()}")
    print(f"  Subjects: {df['subject'].nunique()}")
    print(f"  Samples:  {len(df)}")

if __name__ == "__main__":
    main()