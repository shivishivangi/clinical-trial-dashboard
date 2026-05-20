# Clinical Trial Analysis Report

## Project Plan
- Part 1: Schema design + load_data.py
- Part 2: Frequency summary table
- Part 3: Statistical analysis + boxplot (melanoma, miraclib, PBMC)
- Part 4: Subset queries (baseline samples, project counts, responder breakdown)
- Dashboard: Plotly Dash with 3 sections

## Database Schema

Three-table relational schema using SQLite:

**projects**
- project_id (PK)

**subjects**
- subject_id (PK)
- project_id (FK -> projects)
- condition
- age
- sex
- treatment
- response

**samples**
- sample_id (PK)
- subject_id (FK -> subjects)
- sample_type
- time_from_treatment_start
- b_cell
- cd8_t_cell
- cd4_t_cell
- nk_cell
- monocyte

## Code Structure

- `load_data.py` handles data ingestion. 
  Initializes the SQLite database schema and loads all rows from `data/cell-count.csv` into the three-table relational database.

- `analysis.py` performs all computation and saves all outputs to `outputs/`.

- `dashboard.py` Renders an interactive Plotly Dash dashboard with 3 tabs: Frequency Summary (Part 2), Statistical Analysis (Part 3), and Subset Explorer (Part 4).
  
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

**Dashboard Architecture**
Part 2 data is loaded from the database at startup via `get_all_samples()` and filtered reactively in callbacks via `compute_frequencies()`. Part 3 is computed at startup directly from the database via `get_melanoma_miraclib()` and `create_boxplot()`. Part 4 queries the database directly on each user interaction via `get_baseline_melanoma()`.

For future scale (millions of rows, real-time data updates, or multiple concurrent users), the dashboard will need a caching layer like Redis, or move to a backend API architecture such as FastAPI and Dash.

### SQL Injection Prevention
`get_baseline_melanoma()` uses parameterized queries (`?` placeholders) rather than string interpolation for all user-supplied filter values, preventing SQL injection in interactive callbacks.

## Key Findings

### Part 3: Statistical Analysis

#### Findings
Comparing melanoma patients treated with miraclib (PBMC samples only for all timepoints):

| Population | U-Statistic | P-Value | p < 0.05 | Bonferroni (p < 0.01) |
|---|---|---|---|---|
| cd4_t_cell | 515277.5 | 0.0133 | Yes | No |
| b_cell | 459968.0 | 0.0557 | No | No |
| nk_cell | 464546.5 | 0.1211 | No | No |
| monocyte | 466509.0 | 0.1631 | No | No |
| cd8_t_cell | 478175.5 | 0.6391 | No | No |

CD4 T cells show a statistically significant difference in relative frequency between responders and non-responders (p = 0.0133, Mann-Whitney U test). No other populations reached significance at the p < 0.05 threshold. This suggests CD4 T cell frequency may be a potential predictor of response to miraclib in melanoma patients. 

However, cd4_t_cell does not survive Bonferroni correction for multiple comparisons (adjusted threshold p < 0.01). This result should be treated as a hypothesis worth investigating rather than a definitive biomarker claim.

#### Methodology

**Mann-Whitney U over t-test**
Cell frequency data does not follow a normal distribution in either response group, confirmed by Shapiro-Wilk normality tests (responders: p < 0.0001, non-responders: p = 0.0054). Mann-Whitney U is the appropriate non-parametric alternative because it makes no normality assumption and is valid for the group sizes observed (responders n=993, non-responders n=975).

**No multiple testing correction**
Five simultaneous hypothesis tests (one per population) were run against p < 0.05. Technically, a Bonferroni or Benjamini-Hochberg correction should be applied when testing multiple hypotheses simultaneously. For this exploratory analysis, uncorrected p-values are acceptable, but the Bonferroni-adjusted threshold (p < 0.01) is also reported. cd4_t_cell does not survive this stricter threshold, which Bob should acknowledge when presenting findings to Yah.

**Relative frequencies over raw counts**
Total cell counts vary across samples due to biological variability. A sample with 100,000 total cells and 10,000 B cells is immunologically equivalent to one with 50,000 total cells and 5,000 B cells since both are 10% B cells. Normalizing to relative frequency makes samples comparable regardless of total count, which is the correct unit for comparing immune composition between groups.

**All timepoints included in Part 3**
Part 3 does not filter to baseline and uses all timepoints for melanoma miraclib PBMC patients. Part 3's comparison captures the full treatment trajectory, while Part 4 is a baseline snapshot (time_from_treatment_start = 0). 

### Part 4: Subset Analysis
Baseline melanoma PBMC samples treated with miraclib (time_from_treatment_start = 0):

| Project | Sample Count |
|---|---|
| prj1 | 384 |
| prj3 | 272 |

| Response | Subject Count |
|---|---|
| Responders | 331 |
| Non-responders | 325 |

| Sex | Subject Count |
|---|---|
| Male | 344 |
| Female | 312 |

Note: prj2 has no melanoma miraclib baseline PBMC samples.

**Average B cell count for melanoma male responders at time=0: 10401.28**

Calculated using this SQL query:

```sql
SELECT ROUND(AVG(s.b_cell), 2)
FROM samples s
JOIN subjects sub ON s.subject_id = sub.subject_id
WHERE sub.condition = 'melanoma'
AND sub.treatment = 'miraclib'
AND sub.sex = 'M'
AND sub.response = 'yes'
AND s.time_from_treatment_start = 0
AND s.sample_type = 'PBMC'
```

## Challenges

### Avoiding duplicate inserts on rerun

The initial implementation of `load_data.py` used `to_sql(if_exists="append")` which would duplicate data if script was run more than once. Refactored to use raw sqlite3 with `INSERT OR IGNORE`, which silently skips rows that already exist based on primary key. This makes the script safe to rerun and easy to append new CSV files without duplicating existing data in the future.

### Percentage Rounding Precision
Initially rounded percentages to 2 decimal places in `analysis.py`, causing per-sample sums of 99.99% instead of 100%. Switched to 4 decimal places which reduces floating point error to 99.9999% (acceptable for clinical reporting).

### Incorrect B Cell Average (Caught and Corrected)
The initial SQL query used to verify the average B cell count for melanoma male responders at time=0 was missing `AND sub.treatment = 'miraclib'`, returning 10206.72 instead of the correct 10401.28. The spec states "among these samples" meaning the miraclib filter needs to be included. Caught this mistake by cross-checking against the implemented part 4 analysis.py output.

### dcc.Graph over iframe for Boxplot
Originally planned to embed `outputs/boxplot.html` as an iframe in the dashboard. Switched to calling `create_boxplot()` directly and passing the figure to `dcc.Graph` for cleaner Dash integration, consistent dark theme styling, and no file dependency at runtime.

### SQLite Connection Management in Dashboard Callbacks
Static data is loaded once at startup with a single connection that is closed before the app starts. Interactive callbacks open a fresh connection per call and close it immediately after. This avoids thread-safety issues with SQLite's single-writer model without requiring connection pooling.

## Testing

Tests are located in `tests/` and can be run with:
```bash
python tests/test_load_data.py
python tests/test_analysis.py
```

Tests validate database integrity (row counts, foreign keys, null checks) and analysis correctness (output shape, column names, filter behavior) using the provided `cell-count.csv` as reference. These are sanity checks rather than a full test suite given the scope of this project. The goal is to catch data pipeline errors early before they propagate into the analysis.