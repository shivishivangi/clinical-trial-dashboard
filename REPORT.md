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

## Design Decisions

**Database Schema**
I chose a 3 table relational schema since projects have a one-to-many relationship with subjects, and subjects have a one-to-many relationship with samples. If response status were stored per sample row, it would repeat identically across all 3 timepoints per patient and risk inconsistency.

Cell counts are stored as columns rather than rows in the samples table because the dataset has exactly 5 fixed populations. This simplifies frequency calculations and maps directly to pandas without requiring any reshaping from the database.

This decision allows for scalability as the number of projects, subjects per project, and samples per patient increase in the long-term. It prevents unnecessary duplication and keeps groups separated but connected through foreign keys. 

**Plotly Dash for dashboard**
Dash produces a professional, multi-section interactive dashboard that runs locally via a single Python script. 

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
