# Clinical Trial Analysis Report

## Project Plan
- Part 1: Schema design + load_data.py
- Part 2: Frequency summary table
- Part 3: Statistical analysis + boxplot (melanoma, miraclib, PBMC)
- Part 4: Subset queries (baseline samples, project counts, responder breakdown)
- Dashboard: Plotly Dash with 4 sections

## Database Schema

Three-table relational schema using SQLite:

**projects**
- project_id (PK)

**subjects**
- subject_id (PK)
- project_id (FK - projects)
- condition
- age
- sex
- treatment
- response

**samples**
- sample_id (PK)
- subject_id (FK - subjects)
- sample_type
- time_from_treatment_start
- b_cell
- cd8_t_cell
- cd4_t_cell
- nk_cell
- monocyte

## Code Structure

- `load_data.py` handles data ingestion. Initializes the SQLite database schema and loads all rows from `data/cell-count.csv` into the three-table relational database.

- `analysis.py` performs all computation and saves all outputs to `outputs/`.

- `dashboard.py` involves display only for this project scale (see design decisions). Loads pre-computed outputs and renders an interactive Plotly Dash dashboard.

## Design Decisions

### Database Schema

I chose a 3 table relational schema since projects have a one-to-many relationship with subjects, and subjects have a one-to-many relationship with samples. If response status were stored per sample row, it would repeat identically across all 3 timepoints per patient and risk inconsistency.

Cell counts are stored as columns rather than rows in the samples table because the dataset has exactly 5 fixed populations. This simplifies frequency calculations and maps directly to pandas without requiring any reshaping from the database.

This decision allows for scalability as the number of projects, subjects per project, and samples per patient increase in the long-term. It prevents unnecessary duplication and keeps groups separated but connected through foreign keys. 

### Raw sqlite3 over SQLAlchemy for data insertion

Initially used SQLAlchemy because of its ORM convenience, however, the final implementation uses raw sqlite3 with a reusable `insert_rows()` helper. This allows for direct control over insert behavior (INSERT OR IGNORE) and removes an unnecessary abstraction layer for a project of this scale (see Challenges). Also, executemany() with raw sqlite3 is significantly faster than SQLAlchemy's to_sql() because it batches all inserts in one operation instead of row by row.

### Percentage Rounding

Relative frequencies are rounded to 4 decimal places. Due to floating point arithmetic, percentages for a given sample may sum to 99.9999% rather than exactly 100%. This is expected behavior and standard in clinical and scientific reporting.

### Plotly Dash for dashboard

Dash produces a professional, multi-section interactive dashboard that runs locally via a single Python script. 

**Dashboard Architecture (Display-Only)**
The dashboard reads pre-computed outputs from `outputs/` rather than querying the database directly; separating computation from presentation.

Interactive features (dropdowns, filters) operate on the loaded dataframe in memory, which is fast given the current dataset size (52,500 rows).

For future scale (millions of rows, real-time data updates, or multiple concurrent users), the dashboard will need to query the database directly with lazy loading, use a caching layer like Redis, or move to a backend API architecture such as FastAPI and Dash.

## Key Findings

**Average B cell count for melanoma male responders at time=0: 10206.72**

Calculated using this SQL query:

```sql
SELECT ROUND(AVG(s.b_cell), 2)
FROM samples s
JOIN subjects sub ON s.subject_id = sub.subject_id
WHERE sub.condition = 'melanoma'
AND sub.sex = 'M'
AND sub.response = 'yes'
AND s.time_from_treatment_start = 0
AND s.sample_type = 'PBMC'
```

## Challenges

### Avoiding duplicate inserts on rerun

The initial implementation of `load_data.py` used `to_sql(if_exists="append")` which would duplicate data if script was run more than once. Refactored to use raw sqlite3 with `INSERT OR IGNORE`, which silently skips rows that already exist based on primary key. This makes the script safe to rerun and easy to append new CSV files without duplicating existing data in the future.

## Testing

Tests are located in `tests/` and can be run with:
```bash
python tests/test_load_data.py
python tests/test_analysis.py
```

Tests validate database integrity (row counts, foreign keys, null checks) and analysis correctness (output shape, column names, filter behavior) using the provided `cell-count.csv` as reference. These are sanity checks rather than a full test suite given the scope of this project. The goal is to catch data pipeline errors early before they propagate into the analysis.