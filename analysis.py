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
            "significant": "yes" if p_value < 0.05 else "no"
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
        plot_bgcolor="white",
        legend_title="Response"
    )

    return fig


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
    fig.write_html("outputs/boxplot.html")
    fig.write_image("outputs/boxplot.png")
    print(f"  Saved mannwhitney_results.csv and boxplot.html to outputs/")

    conn.close()


if __name__ == "__main__":
    main()