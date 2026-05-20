# Clinical Trial Immune Cell Dashboard

## Overview
A Python-based analysis pipeline and interactive dashboard for exploring immune cell population data from a clinical trial. Analyzes immune cell population dynamics across melanoma and carcinoma patients treated with miraclib or phauximab, and healthy untreated controls.

## How to Run

```bash
git clone https://github.com/shivishivangi/clinical-trial-dashboard.git
cd clinical-trial-dashboard
make setup
make pipeline
make dashboard
```

## Dashboard

Live: https://clinical-trial-dashboard-g1d7.onrender.com

Local: http://localhost:8050 after `make dashboard`

**Tab 1 - Frequency Summary**: filterable table of relative cell population frequencies across all samples. Filter by condition, project, and population.

**Tab 2 - Statistical Analysis**: boxplot comparing responders vs non-responders for each cell population, with Mann-Whitney U test results and significance highlighting.

**Tab 3 - Subset Explorer**: interactive subset analysis. Filter by condition, treatment, sample type, and timepoint to explore sample counts, subject breakdowns, and average cell counts.

## Code Structure

Input data is in `data/cell-count.csv` 
All generated outputs (tables, plots, results) are saved to `outputs/` after running `make pipeline`.

- `load_data.py` initializes the SQLite database schema and loads all rows from `data/cell-count.csv` into a three-table relational database.

- `analysis.py` performs all computation for Parts 2-4. Generates the frequency summary table, runs Mann-Whitney U statistical tests, creates the boxplot, and computes subset queries. Saves all outputs to `outputs/`.

- `dashboard.py` renders an interactive Plotly Dash dashboard with 3 tabs: Frequency Summary (Part 2), Statistical Analysis (Part 3), and Subset Explorer (Part 4).

## Database Schema

Three-table relational schema using SQLite:

**projects** - one row per clinical trial project
- project_id (PK)

**subjects** - one row per patient
- subject_id (PK)
- project_id (FK -> projects)
- condition, age, sex, treatment, response

**samples** - one row per biological sample
- sample_id (PK)
- subject_id (FK -> subjects)
- sample_type, time_from_treatment_start
- b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte

### Rationale

Subject-level attributes (condition, sex, treatment, response) are stored in `subjects` rather than repeated in every sample row. This eliminates redundancy since a patient's response status is the same across all their samples and should only be stored once.

Cell counts are stored as columns in `samples` because the dataset has exactly 5 fixed populations. This simplifies frequency calculations and maps directly to pandas without reshaping.

### Scalability

| Scenario | How schema handles it |
|---|---|
| Hundreds of projects | New rows in `projects` - no schema change |
| Thousands of samples | New rows in `samples` - no schema change |
| New timepoints | New rows in `samples` - no schema change |
| New cell populations | Add columns to `samples` for fixed sets; migrate to a long-format `cell_counts` table if populations become dynamic |
| New sample types | `sample_type` column already handles it |
| Analytics by project/condition/treatment | Foreign key joins across 3 tables |

For hundreds of projects and thousands of samples, adding indexes on frequently filtered columns (`condition`, `treatment`, `response`, `time_from_treatment_start`) would maintain query performance. For millions of rows or real-time analytics, migrating to PostgreSQL would be the best decision.

## Running Tests

```bash
make test
```

Unit and integration tests covering database integrity, analysis correctness, and parameterized query behavior.