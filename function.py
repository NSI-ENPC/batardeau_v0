from dash import Dash, html, dash_table, dcc, callback, Output, Input, State
import plotly.graph_objects as go
import dash_bootstrap_components as dbc

style_label = {}#{'font-size':18, 'font-family':"Arial"}
style_graph = {'width': '90vh', 'height': '80vh'}
font_style_graph = dict(family="Arial",size=18)

def create_onglets(**kwargs):
    lst_tabs_ = kwargs.get('lst_tabs_', [])
    lst_label_ = kwargs.get("lst_label_", [])
    id_ = kwargs.get('id_', '')

    nav = html.Div([
        html.Div(
            [
                dbc.RadioItems(
                    id=id_,
                    className="btn-group",
                    inputClassName="btn-check",
                    labelClassName="btn btn-outline-primary",
                    labelCheckedClassName="active",
                    options=[{"label": label_, "value": k} for k, label_ in enumerate(lst_label_)],
                    value=0,
                ),
                html.H2(),
            ],
            className="radio-group",
        ),
    ]+lst_tabs_)

    return nav

def calcul_button(**kwargs):
    name_app_ = kwargs.get('name_app_', '')
    button = html.Div(
        [
            dbc.Button("Calculer", color="primary", id='{}-calcul-button'.format(name_app_)),
        ],
        className="d-grid gap-2",
    )
    return button

def general_window(**kwargs):
    tabs_ = kwargs.get('tabs_', None)
    graphs_ = kwargs.get('graphs_', None)
    titre_ = kwargs.get('titre_', None)
    name_app_ = kwargs.get('name_app_', '')
    return html.Div([
        calcul_button(name_app_=name_app_),
        html.Br(),
        html.H2(titre_, style={'textAlign': 'center'}),
        html.Br(),
        html.Div([
                dbc.Row([dbc.Col(tabs_, width=5), dbc.Col(graphs_, width=6),], justify="center"),
            ]),
        ])

def create_graph(**kwargs):
    id_ = kwargs.get('id_',None)
    fig = go.Figure()
    fig.update_layout(plot_bgcolor='white',
                      font=font_style_graph,
                      height=700)

    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey',
        zeroline=True, 
        zerolinewidth=1,
        zerolinecolor='black'
    )
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey',
        zeroline=True, 
        zerolinewidth=1,
        zerolinecolor='black',
        scaleanchor = "x",
        scaleratio = 1,
    )
    return dcc.Graph(figure=fig, id=id_,style=style_graph)

def update_layout(**kwargs):
    fig = kwargs.get('figure', go.Figure())
    titre_x_ = kwargs.get('titre_x_', '')
    titre_y_ = kwargs.get('titre_y_', '')
    legende_ = kwargs.get('legende_', '')
    equal_ = kwargs.get('equal_', False)

    fig.update_layout(xaxis_title=titre_x_,
                      yaxis_title=titre_y_,
                      legend_title=legende_,
                      plot_bgcolor='white',
                      font=font_style_graph,
                      height=700)

    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey',
        zeroline=True, 
        zerolinewidth=1,
        zerolinecolor='black'
    )
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey',
        zeroline=True, 
        zerolinewidth=1,
        zerolinecolor='black',
    )
    if equal_:
        fig.update_yaxes(scaleanchor = "x", scaleratio = 1)
    

def button_number(**kwargs):
    id_ = kwargs.get('id_', None)
    value_ = kwargs.get('value_', None)
    label_ = kwargs.get('label_', '')
    hidden_ = kwargs.get('hidden_', False)

    style_ = {'display':'none'} if hidden_ else {'display':'block'}

    butt = dbc.Col(html.Div([
        dcc.Markdown(label_, style=style_label, mathjax=True),
        dbc.Input(id=id_, value=value_, type='number'),
        html.H2()
    ],id=id_+'-box',style=style_))

    return butt

def button_text(**kwargs):
    id_ = kwargs.get('id_', None)
    value_ = kwargs.get('value_', None)
    label_ = kwargs.get('label_', '')
    hidden_ = kwargs.get('hidden_', False)

    style_ = {'display':'none'} if hidden_ else {'display':'block'}

    butt = html.Div([
        dcc.Markdown(label_, style=style_label, mathjax=True),
        dbc.Input(id=id_, value=value_, type='text'),
        html.H2()
    ],id=id_+'-box',style=style_)

    return butt

def dropdown(**kwargs):
    id_ = kwargs.get('id_', None)
    options_ = kwargs.get('options_', [])
    value_ = kwargs.get('value_', None)
    label_ = kwargs.get('label_', '')
    hidden_ = kwargs.get('hidden_', False)

    style_ = {'display':'none'} if hidden_ else {'display':'block'}

    drop = html.Div([
        dcc.Markdown(label_, style=style_label, mathjax=True),
        dcc.Dropdown(id=id_, options=options_, value=value_),
        html.H2()
    ],id=id_+'-box',style=style_)

    return drop