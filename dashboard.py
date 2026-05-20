import sqlite3
import pandas as pd
import dash
from dash import dcc, html, dash_table, Input, Output
import dash_bootstrap_components as dbc

from analysis import (
    DB_PATH, CELL_COLS, POPULATION_LABELS,
    get_distinct_values,
    get_all_samples, compute_frequencies,
    get_melanoma_miraclib, run_mannwhitney, create_boxplot,
    get_baseline_melanoma, compute_subset_summary, compute_avg_cell,
)

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])
server = app.server

# Load static data at startup to minimize repeated DB queries in callbacks
_conn = sqlite3.connect(DB_PATH)

conditions   = get_distinct_values(_conn, "subjects", "condition")
projects     = get_distinct_values(_conn, "projects", "project_id")
treatments   = get_distinct_values(_conn, "subjects", "treatment")
timepoints   = get_distinct_values(_conn, "samples", "time_from_treatment_start")
sample_types = get_distinct_values(_conn, "samples", "sample_type")
sexes        = get_distinct_values(_conn, "subjects", "sex")

df_all_samples = get_all_samples(_conn)

_df_mel        = get_melanoma_miraclib(_conn)
_df_freq_p3    = compute_frequencies(_df_mel, extra_id_vars=["response"])
_fig_box       = create_boxplot(_df_freq_p3)
_df_stats      = run_mannwhitney(_df_freq_p3).copy()
_df_stats["population"] = _df_stats["population"].map(POPULATION_LABELS)

_conn.close()

# Helpers
def opts(values):
    return [{"label": str(v), "value": v} for v in values]

def pop_opts():
    return [{"label": POPULATION_LABELS[c], "value": c} for c in CELL_COLS]

TABLE_STYLE = dict(
    style_cell={"textAlign": "left", "padding": "8px", "fontFamily": "sans-serif",
                "backgroundColor": "#1e2130", "color": "#ffffff", "border": "1px solid #2d3748"},
    style_header={"fontWeight": "bold", "backgroundColor": "#0f1117",
                  "color": "#00d4aa", "border": "1px solid #2d3748"},
    style_data={"border": "1px solid #2d3748"},
    style_filter={"backgroundColor": "#2d3748", "color": "#ffffff", "border": "1px solid #4a5568"},
)

# Layout
app.layout = dbc.Container(
    style={"backgroundColor": "#0f1117", "minHeight": "100vh"},
    fluid=True,
    children=[
    html.H2("Clinical Trial Dashboard", className="my-4 text-center", style={"color": "#00d4aa"}),

    dbc.Tabs([

        # Tab 1: Frequency Summary (Part 2)
        dbc.Tab(label="Frequency Summary", children=[
            dbc.Row([
                dbc.Col([
                    html.Label("Condition"),
                    dcc.Dropdown(
                        id="p2-condition",
                        options=[{"label": "All", "value": "all"}] + opts(conditions),
                        value="all", clearable=False, style={"color": "#000000", "backgroundColor": "#ffffff"},
                    ),
                ], width=3),
                dbc.Col([
                    html.Label("Project"),
                    dcc.Dropdown(
                        id="p2-project",
                        options=[{"label": "All", "value": "all"}] + opts(projects),
                        value="all", clearable=False, style={"color": "#000000", "backgroundColor": "#ffffff"},
                    ),
                ], width=3),
                dbc.Col([
                    html.Label("Population"),
                    dcc.Dropdown(
                        id="p2-population",
                        options=[{"label": "All", "value": "all"}] + pop_opts(),
                        value="all", clearable=False, style={"color": "#000000", "backgroundColor": "#ffffff"},
                    ),
                ], width=3),
            ], className="my-3"),
            dash_table.DataTable(
                id="p2-table",
                page_size=200,
                sort_action="native",
                style_table={"overflowX": "auto"},
                **TABLE_STYLE,
            ),
        ]),

        # Tab 2: Statistical Analysis (Part 3)
        dbc.Tab(label="Statistical Analysis", children=[
            html.P(
                "Melanoma patients treated with miraclib — PBMC samples, all timepoints.",
                className="text-muted mt-3",
            ),
            dbc.Row([
                dbc.Col(dcc.Graph(figure=_fig_box), width=8),
                dbc.Col([
                    html.H5("Mann-Whitney U Results", className="mt-2"),
                    dash_table.DataTable(
                        data=_df_stats.to_dict("records"),
                        columns=[
                            {"name": "Population",  "id": "population"},
                            {"name": "U Statistic", "id": "u_statistic"},
                            {"name": "P Value",     "id": "p_value"},
                            {"name": "p < 0.05",    "id": "significant_p05"},
                            {"name": "Bonferroni",  "id": "significant_bonferroni"},
                        ],
                        style_data_conditional=[{
                            "if": {"filter_query": '{significant_p05} = "yes"'},
                            "backgroundColor": "#1a4731",
                            "color": "#ffffff",
                        }],
                        **TABLE_STYLE,
                    ),
                ], width=4),
            ], className="my-3"),
        ]),

        # Tab 3: Subset Explorer (Part 4)
        dbc.Tab(label="Subset Explorer", children=[
            html.H5("Subset Filters", className="mt-3"),
            dbc.Row([
                dbc.Col([
                    html.Label("Condition"),
                    dcc.Dropdown(
                        id="p4-condition",
                        options=[{"label": "All", "value": "all"}] + opts(conditions),
                        value="melanoma", clearable=False, style={"color": "#000000", "backgroundColor": "#ffffff"},
                    ),
                ], width=3),
                dbc.Col([
                    html.Label("Treatment"),
                    dcc.Dropdown(
                        id="p4-treatment",
                        options=[{"label": "All", "value": "all"}] + opts(treatments),
                        value="miraclib", clearable=False, style={"color": "#000000", "backgroundColor": "#ffffff"},
                    ),
                ], width=3),
                dbc.Col([
                    html.Label("Sample Type"),
                    dcc.Dropdown(id="p4-sample-type", options=opts(sample_types), value="PBMC", clearable=False, style={"color": "#000000", "backgroundColor": "#ffffff"}),
                ], width=2),
                dbc.Col([
                    html.Label("Timepoint"),
                    dcc.Dropdown(
                        id="p4-timepoint",
                        options=[{"label": "All", "value": "all"}] + opts(timepoints),
                        value=0, clearable=False, style={"color": "#000000", "backgroundColor": "#ffffff"},
                    ),
                ], width=2),
            ], className="mb-3"),

            dbc.Row([
                dbc.Col([
                    html.H6("Samples per Project"),
                    dash_table.DataTable(id="p4-project-table", **TABLE_STYLE),
                ], width=4),
                dbc.Col([
                    html.H6("Subjects by Response"),
                    dash_table.DataTable(id="p4-response-table", **TABLE_STYLE),
                ], width=4),
                dbc.Col([
                    html.H6("Subjects by Sex"),
                    dash_table.DataTable(id="p4-sex-table", **TABLE_STYLE),
                ], width=4),
            ], className="mb-4"),

            html.Hr(),

            html.H5("Average Cell Count Calculator"),
            dbc.Row([
                dbc.Col([
                    html.Label("Cell Population"),
                    dcc.Dropdown(id="p4-cell-col", options=pop_opts(), value="b_cell", clearable=False, style={"color": "#000000", "backgroundColor": "#ffffff"}),
                ], width=3),
                dbc.Col([
                    html.Label("Sex"),
                    dcc.Dropdown(id="p4-sex", options=opts(sexes), value="M", clearable=False, style={"color": "#000000", "backgroundColor": "#ffffff"}),
                ], width=2),
                dbc.Col([
                    html.Label("Response"),
                    dcc.Dropdown(id="p4-response-avg", options=opts(["yes", "no"]), value="yes", clearable=False, style={"color": "#000000", "backgroundColor": "#ffffff"}),
                ], width=2),
            ], className="mb-3"),
            html.Div(id="p4-avg-result"),
        ]),

    ]),
])


# Callbacks

@app.callback(
    Output("p2-table", "data"),
    Output("p2-table", "columns"),
    Input("p2-condition", "value"),
    Input("p2-project", "value"),
    Input("p2-population", "value"),
)
def update_summary_table(condition, project, population):
    df = compute_frequencies(
        df_all_samples,
        project_id=None if project == "all" else project,
        condition=None if condition == "all" else condition,
    )
    if population != "all":
        df = df[df["population"] == population]
    cols = [{"name": c, "id": c} for c in df.columns]
    return df.to_dict("records"), cols


@app.callback(
    Output("p4-project-table", "data"),
    Output("p4-project-table", "columns"),
    Output("p4-response-table", "data"),
    Output("p4-response-table", "columns"),
    Output("p4-sex-table", "data"),
    Output("p4-sex-table", "columns"),
    Input("p4-condition", "value"),
    Input("p4-treatment", "value"),
    Input("p4-sample-type", "value"),
    Input("p4-timepoint", "value"),
)
def update_subset_tables(condition, treatment, sample_type, timepoint):
    conn = sqlite3.connect(DB_PATH)
    df = get_baseline_melanoma(
        conn,
        condition=None if condition == "all" else condition,
        treatment=None if treatment == "all" else treatment,
        sample_type=sample_type,
        time_point=None if timepoint == "all" else int(timepoint),
    )
    conn.close()

    if df.empty:
        empty = [{"result": "No data for selected filters."}]
        empty_cols = [{"name": "Result", "id": "result"}]
        return empty, empty_cols, empty, empty_cols, empty, empty_cols

    samples_per_project, response_counts, sex_counts = compute_subset_summary(df)

    def to_table(d):
        return d.to_dict("records"), [{"name": c, "id": c} for c in d.columns]

    pr, pc = to_table(samples_per_project)
    rr, rc = to_table(response_counts)
    sr, sc = to_table(sex_counts)
    return pr, pc, rr, rc, sr, sc


@app.callback(
    Output("p4-avg-result", "children"),
    Input("p4-condition", "value"),
    Input("p4-treatment", "value"),
    Input("p4-sample-type", "value"),
    Input("p4-timepoint", "value"),
    Input("p4-cell-col", "value"),
    Input("p4-sex", "value"),
    Input("p4-response-avg", "value"),
)
def update_avg(condition, treatment, sample_type, timepoint, cell_col, sex, response):
    conn = sqlite3.connect(DB_PATH)
    df = get_baseline_melanoma(
        conn,
        condition=None if condition == "all" else condition,
        treatment=None if treatment == "all" else treatment,
        sample_type=sample_type,
        time_point=None if timepoint == "all" else int(timepoint),
    )
    conn.close()

    avg = compute_avg_cell(df, cell_col=cell_col, sex=sex, response=response)

    if pd.isna(avg):
        return dbc.Alert("No data for selected filters.", color="warning")

    label = POPULATION_LABELS.get(cell_col, cell_col)
    return dbc.Alert(
        f"Average {label} count — {sex}, response={response}, timepoint={timepoint}: {avg:.2f}",
        style={"backgroundColor": "#00d4aa", "color": "#0f1117", "fontWeight": "bold"},
    )


if __name__ == "__main__":
    app.run(debug=False, port=8050)