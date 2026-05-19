import pandas as pd
from sqlalchemy import create_engine, text

engine = create_engine("sqlite:///clinical_trial.db")

with engine.connect() as conn:
    conn.execute(text("""
        CREATE TABLE IF NOT EXISTS projects (
            project_id TEXT PRIMARY KEY
        )
    """))

    conn.execute(text("""
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
    """))

    conn.execute(text("""
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
    """))
    conn.commit()

df = pd.read_csv("data/cell-count.csv")

with engine.connect() as conn:

    # projects table
    # unique project IDs 
    projects = df[["project"]].drop_duplicates()
    projects.columns = ["project_id"]
    projects.to_sql("projects", conn, if_exists="append", index=False)

    # subjects table
    subjects = df[[
        "subject", "project", "condition",
        "age", "sex", "treatment", "response"
    ]].drop_duplicates(subset=["subject"])
    subjects.columns = [
        "subject_id", "project_id", "condition",
        "age", "sex", "treatment", "response"
    ]
    subjects.to_sql("subjects", conn, if_exists="append", index=False)

    # samples table
    # one row per sample, with subject_id as foreign key
    samples = df[[
        "sample", "subject", "sample_type",
        "time_from_treatment_start",
        "b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"
    ]]
    samples.columns = [
        "sample_id", "subject_id", "sample_type",
        "time_from_treatment_start",
        "b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"
    ]
    samples.to_sql("samples", conn, if_exists="append", index=False)

    conn.commit()

print("Database created and data loaded successfully.")
print(f"  Projects: {df['project'].nunique()}")
print(f"  Subjects: {df['subject'].nunique()}")
print(f"  Samples:  {len(df)}")