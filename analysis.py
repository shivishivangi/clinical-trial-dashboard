import pandas as pd
import sqlite3
import os
from scipy import stats
import plotly.graph_objects as go

os.makedirs("outputs", exist_ok=True)

DB_PATH = "clinical_trial.db"
CELL_COLS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

POPULATION_LABELS = {
    "b_cell": "B Cell",
    "cd8_t_cell": "CD8 T Cell",
    "cd4_t_cell": "CD4 T Cell",
    "nk_cell": "NK Cell",
    "monocyte": "Monocyte",
}


def get_distinct_values(conn, table, column):
    """
    Fetch sorted distinct values for a column from a given table.
    Used by the dashboard to populate dropdowns dynamically.
    table and column must be trusted internal constants — not user input.
    """
    result = pd.read_sql(
        f"SELECT DISTINCT {column} FROM {table} ORDER BY {column}",
        conn
    )
    return result[column].tolist()


def get_all_samples(conn):
    """
    Fetch all samples with cell counts, condition, and project via a 3-table join. 
    Returns one row per sample across all cohorts.
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


def compute_frequencies(df, project_id=None, condition=None, extra_id_vars=None):
    """
    Compute relative frequency of each cell population per sample.
    Optional filters: project_id, condition.
    extra_id_vars: additional columns to preserve through the melt (e.g. ["response"]).
    Returns long-format dataframe with columns:
    sample, total_count, population, count, percentage
    """
    if extra_id_vars is None:
        extra_id_vars = []

    df = df.copy()

    if project_id:
        df = df[df["project_id"] == project_id]
    if condition:
        df = df[df["condition"] == condition]

    df["total_count"] = df[CELL_COLS].sum(axis=1)
    if (df["total_count"] == 0).any():
        print("WARNING: samples with zero total count detected")

    df_long = df.melt(
        id_vars=["sample_id", "total_count"] + extra_id_vars,
        value_vars=CELL_COLS,
        var_name="population",
        value_name="count"
    )

    df_long["percentage"] = (
        df_long["count"] / df_long["total_count"] * 100
    ).round(4)

    df_long = df_long.rename(columns={"sample_id": "sample"})
    base_cols = ["sample", "total_count"] + extra_id_vars + ["population", "count", "percentage"]
    df_long = df_long[base_cols]
    df_long = df_long.sort_values(["sample", "population"]).reset_index(drop=True)

    return df_long


def get_melanoma_miraclib(conn):
    return pd.read_sql("""
        SELECT
            s.sample_id,
            s.b_cell, s.cd8_t_cell, s.cd4_t_cell, s.nk_cell, s.monocyte,
            s.sample_type,
            sub.response,
            sub.condition,
            sub.treatment
        FROM samples s
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.condition = 'melanoma'
        AND sub.treatment = 'miraclib'
        AND s.sample_type = 'PBMC'
        AND sub.response IN ('yes', 'no')
    """, conn)


def run_mannwhitney(df):
    results = []

    for population in CELL_COLS:
        pop_data = df[df["population"] == population]
        responders = pop_data[pop_data["response"] == "yes"]["percentage"]
        non_responders = pop_data[pop_data["response"] == "no"]["percentage"]

        u_stat, p_value = stats.mannwhitneyu(
            responders,
            non_responders,
            alternative="two-sided"
        )

        results.append({
            "population": population,
            "u_statistic": round(u_stat, 4),
            "p_value": round(p_value, 4),
            "significant_p05": "yes" if p_value < 0.05 else "no",
            "significant_bonferroni": "yes" if p_value < (0.05 / len(CELL_COLS)) else "no"
        })

    return pd.DataFrame(results).sort_values("p_value")


def create_boxplot(df):
    fig = go.Figure()

    colors = {"yes": "#2ecc71", "no": "#e74c3c"}
    labels = {"yes": "Responder", "no": "Non-Responder"}

    for response in ["yes", "no"]:
        df_response = df[df["response"] == response]

        fig.add_trace(go.Box(
            x=df_response["population"].map(POPULATION_LABELS),
            y=df_response["percentage"],
            name=labels[response],
            marker_color=colors[response],
            boxmean=True
        ))

    fig.update_layout(
        title="Cell Population Frequencies: Responders vs Non-Responders",
        xaxis_title="Cell Population",
        yaxis_title="Relative Frequency (%)",
        boxmode="group",
        plot_bgcolor="#1e2130",
        paper_bgcolor="#1e2130",
        font={"color": "#ffffff"},
        legend={"title": "Response", "bgcolor": "#1e2130"},
    )

    return fig


def get_baseline_melanoma(conn, condition='melanoma', treatment='miraclib',
                          sample_type='PBMC', time_point=0):
    filters = ["s.sample_type = ?"]
    params  = [sample_type]

    if condition is not None:
        filters.append("sub.condition = ?")
        params.append(condition)
    if treatment is not None:
        filters.append("sub.treatment = ?")
        params.append(treatment)
    if time_point is not None:
        filters.append("s.time_from_treatment_start = ?")
        params.append(time_point)

    where = " AND ".join(filters)

    return pd.read_sql(f"""
        SELECT
            s.sample_id,
            s.b_cell, s.cd8_t_cell, s.cd4_t_cell, s.nk_cell, s.monocyte,
            s.time_from_treatment_start,
            sub.subject_id,
            sub.response,
            sub.sex,
            sub.condition,
            sub.treatment,
            p.project_id
        FROM samples s
        JOIN subjects sub ON s.subject_id = sub.subject_id
        JOIN projects p ON sub.project_id = p.project_id
        WHERE {where}
    """, conn, params=params)


def compute_subset_summary(df):
    """
    Compute summary counts from a filtered sample dataframe.
    Returns three dataframes:
    - samples_per_project: sample count grouped by project
    - response_counts: unique subject count grouped by response (yes/no)
    - sex_counts: unique subject count grouped by sex (M/F)

    Expects df to have columns: project_id, subject_id, response, sex.
    Call get_baseline_melanoma() or any parameterized query first to filter
    the data, then pass the result here.

    Dashboard usage:
        df = get_baseline_melanoma(conn, condition=selected_condition,
                                        treatment=selected_treatment,
                                        time_point=selected_timepoint)
        samples_per_project, response_counts, sex_counts = compute_subset_summary(df)
    """
    samples_per_project = (
        df.groupby("project_id")
        .size()
        .reset_index(name="sample_count")
    )

    subjects = df.drop_duplicates("subject_id")

    response_counts = (
        subjects.groupby("response")
        .size()
        .reset_index(name="subject_count")
    )

    sex_counts = (
        subjects.groupby("sex")
        .size()
        .reset_index(name="subject_count")
    )

    return samples_per_project, response_counts, sex_counts


def compute_avg_cell(df, cell_col='b_cell', sex='M', response='yes'):
    mask = (df["sex"] == sex) & (df["response"] == response)
    return round(df[mask][cell_col].mean(), 2)


def main():
    conn = sqlite3.connect(DB_PATH)

    print("Running Part 2: Frequency Summary Table")
    df_samples = get_all_samples(conn)
    df_summary = compute_frequencies(df_samples)
    df_summary.to_csv("outputs/summary_table.csv", index=False)
    print(f"  Saved {len(df_summary)} rows to outputs/summary_table.csv")

    print("Running Part 3: Statistical Analysis")
    df_mel = get_melanoma_miraclib(conn)
    df_freq = compute_frequencies(df_mel, extra_id_vars=["response"])
    df_stats = run_mannwhitney(df_freq)
    df_stats.to_csv("outputs/mannwhitney_results.csv", index=False)
    fig = create_boxplot(df_freq)
    fig.write_image("outputs/boxplot.png")
    print(f"  Saved mannwhitney_results.csv and boxplot.png to outputs/")

    print("\nRunning Part 4: Subset Analysis")
    df_baseline = get_baseline_melanoma(conn)
    samples_per_project, response_counts, sex_counts = compute_subset_summary(df_baseline)
    avg_bcell = compute_avg_cell(df_baseline, cell_col='b_cell', sex='M', response='yes')

    with open("outputs/part4_results.txt", "w") as f:
        f.write("PART 4: Baseline Melanoma Miraclib PBMC Samples\n")
        f.write("=" * 50 + "\n\n")
        f.write("Samples per project:\n")
        f.write(samples_per_project.to_string(index=False))
        f.write("\n\nResponders vs Non-responders:\n")
        f.write(response_counts.to_string(index=False))
        f.write("\n\nMales vs Females:\n")
        f.write(sex_counts.to_string(index=False))
        f.write(f"\n\nAverage B cell count (melanoma male responders, time=0):\n{avg_bcell}\n")

    print(f"  Samples per project:\n{samples_per_project.to_string(index=False)}")
    print(f"\n  Responders vs Non-responders:\n{response_counts.to_string(index=False)}")
    print(f"\n  Males vs Females:\n{sex_counts.to_string(index=False)}")
    print(f"\n  Avg B cell (melanoma male responders, time=0): {avg_bcell}")

    conn.close()


if __name__ == "__main__":
    main()