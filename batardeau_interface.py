from dash import Dash, html, dash_table, dcc, callback, Output, Input, State, dash_table
import pandas as pd
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import numpy as np
import os
import time

from function import button_number, general_window, create_graph, update_layout, create_onglets, button_text, dropdown

name_app = 'batardeau'

app = Dash(__name__, external_stylesheets=[dbc.themes.JOURNAL])

titre_="Batardeau (v0.3)"

app.title=titre_

# Buttons & co.

sout_z_sup = button_number(id_='{}-sout-z-sup'.format(name_app),label_='Cote supérieure (mNGF)')
sout_z_inf = button_number(id_='{}-sout-z-inf'.format(name_app),label_='Cote inférieure (mNGF)')
sout_lv = button_number(id_='{}-sout-lv'.format(name_app),label_=r'Taille verticale caractéritique $l_v$ (m)')
sout_lh = button_number(id_='{}-sout-lh'.format(name_app),label_=r'Taille horizontale caractéristique $l_h$ (m)')
sout_E = button_number(id_='{}-sout-E'.format(name_app),label_=r"Module d'Young $E$ (MPa)")
sout_nu = button_number(id_='{}-sout-nu'.format(name_app),label_=r'Coefficient de Poisson $\nu$')
sout_t = button_number(id_='{}-sout-t'.format(name_app),label_=r'Épaisseur $t$ (m)')

post_defo = dropdown(id_='{}-post-defo'.format(name_app),options_=['Oui','Non'],label_='Afficher le maillage déformé :',value_='Oui')
post_colorize = dropdown(id_='{}-post-colorize'.format(name_app),options_=['ux','uy','uz','|u|','mhh','mvv','mhv'],label_="Colorer par :",value_='|u|')
scaleF = button_number(id_='{}-scaleFactor'.format(name_app),label_='Scale factor for deformation',value_=1)

# Tableaux

sout_column_name = ["X [m]", "Y [m]"]
sout_column_id = ["x", "y"]
sout_line_init = [0,0]

sout_column_dico = [{
            'name': sout_column_name[i],
            'id': sout_column_id[i],
            'deletable': False,
            'renamable': False,
            'type':'numeric'
        } for i in range(2)]
#SM_column_dico[2]['presentation'] = 'dropdown'

sout_tab = html.Div([
    html.H6("Coordonnées du contour (horaire)"),
    dash_table.DataTable(
        id='{}-sout-tab'.format(name_app),
        columns=sout_column_dico,
        data=[
            {sout_column_id[i]: sout_line_init[i] for i in range(2)}
            for j in range(1)
        ],
        editable=True,
        row_deletable=True
    ),

    dbc.Button('Add Row', id='{}-sout-addRow'.format(name_app), n_clicks=0),
], id='{}-div-SM'.format(name_app))

# Onglets

onglet_structure = html.Div([
    dbc.Accordion([
        dbc.AccordionItem([
            html.H6("Caractéritiques principales"),
            dbc.Row([sout_z_inf,sout_z_sup]),
            #html.Br(),
            dbc.Row([sout_lh,sout_lv]),
            html.Br(),
            dbc.Row([sout_E,sout_nu,sout_t]),
            html.Br(),
            sout_tab,
        ],title='Soutènement'),
        dbc.AccordionItem([
            dcc.Markdown('Cette section sera complétée lorsque le code permettra de prendre en compte les fonctionnalités associées.')
        ],title='Liernes'),
        dbc.AccordionItem([
            dcc.Markdown('Cette section sera complétée lorsque le code permettra de prendre en compte les fonctionnalités associées.')
        ],title='Butons'),
        dbc.AccordionItem([
            dcc.Markdown('Cette section sera complétée lorsque le code permettra de prendre en compte les fonctionnalités associées.')
        ],title='Tirants'),
    ])
],id='{}-div-structure'.format(name_app))
onglet_sol = html.Div([
    dcc.Markdown('Cette section sera complétée lorsque le code permettra de prendre en compte les fonctionnalités associées.')
],id='{}-div-sol'.format(name_app))
onglet_eau = html.Div([
    dcc.Markdown('Cette section sera complétée lorsque le code permettra de prendre en compte les fonctionnalités associées.')
],id='{}-div-eau'.format(name_app))
onglet_postproc = html.Div([dbc.Row([post_defo,post_colorize,scaleF])],id='{}-div-postproc'.format(name_app))

model_3D = create_graph(id_='{}-model3D'.format(name_app))
graph_3D = create_graph(id_='{}-graph3D'.format(name_app))
graph_2D = create_graph(id_='{}-graph2D'.format(name_app))

# 1) Partie Données

lst_lab_onglets_left = ["Structure", "Sol", "Eau", "Post-processing"]
lst_id_onglets_left = ['{}-div-structure'.format(name_app),'{}-div-sol'.format(name_app),'{}-div-eau'.format(name_app),'{}-div-postproc'.format(name_app)]

tabs_ = create_onglets(lst_tabs_=[onglet_structure,onglet_sol,onglet_eau,onglet_postproc], lst_label_=lst_lab_onglets_left, id_='{}-onglet_left'.format(name_app))

# 2) Partie graphes

lst_lab_onglets_right = ["Modèle", "Résultats 3D", "Résultats 2D"]
lst_id_onglets_right = ['{}-model3D'.format(name_app),'{}-graph3D'.format(name_app),'{}-graph2D'.format(name_app)]

graphs_ = create_onglets(lst_tabs_=[model_3D,graph_3D,graph_2D], lst_label_=lst_lab_onglets_right, id_='{}-onglet_right'.format(name_app))

# 3) Assemblage final

app.layout = general_window(tabs_=tabs_, graphs_=graphs_, titre_=titre_, name_app_=name_app)

# III - Callbacks généraux
# ------------------------

# Affichage navigation gauche
@callback(
        [Output(component_id=id_, component_property='style') for id_ in lst_id_onglets_left],
        Input('{}-onglet_left'.format(name_app), 'value')
)
def affiche_onglets(nav1):
    lst_st = []
    for k in range(len(lst_id_onglets_left)):
        if k==nav1:
            lst_st.append({'display':'block'})
        else:
            lst_st.append({'display':'none'})
    return tuple(lst_st)

# Affichage navigation droite
@callback(
        [Output(component_id=id_, component_property='style') for id_ in lst_id_onglets_right],
        Input('{}-onglet_right'.format(name_app), 'value')
)
def affiche_onglets(nav1):
    lst_st = []
    for k in range(len(lst_id_onglets_right)):
        if k==nav1:
            lst_st.append({'display':'block'})
        else:
            lst_st.append({'display':'none'})
    return tuple(lst_st)

# Add Line for structure
@callback(
    Output('{}-sout-tab'.format(name_app), 'data'),
    Input('{}-sout-addRow'.format(name_app), 'n_clicks'),
    State('{}-sout-tab'.format(name_app), 'data'),
    State('{}-sout-tab'.format(name_app), 'columns'))
def add_row(n_clicks, rows, columns):
    if n_clicks > 0:
        rows.append({c['id']: '' for c in columns})
    return rows

# IV - Callbacks particuliers
# ---------------------------

# Callback d'affichage de la structure (soutènement + éventuelles poutres, etc etc)
@callback(
    Output('{}-model3D'.format(name_app),'figure'),
    Input('{}-sout-z-sup'.format(name_app),'value'),
    Input('{}-sout-z-inf'.format(name_app),'value'),
    Input('{}-sout-tab'.format(name_app), 'data')
)
def aff_struct(zSup,zInf,tabPt):
    fig = go.Figure()

    print()

    try:
        '''        

        #fig.add_trace(go.Surface(x=[0,1],y=[0,1],z=[[0,0],[0,0]]))
        '''
        n = len(tabPt)
        k = 0
        while (k<n):
            if k==n-1:
                fig.add_trace(go.Mesh3d(x=[tabPt[k]['x'], tabPt[k]['x']+0.01, tabPt[0]['x']+0.01, tabPt[0]['x'], tabPt[k]['x'], tabPt[k]['x']+0.01, tabPt[0]['x']+0.01, tabPt[0]['x']],
                                        y=[tabPt[k]['y'], tabPt[k]['y']+0.01, tabPt[0]['y']+0.01, tabPt[0]['y'], tabPt[k]['y'], tabPt[k]['y']+0.01, tabPt[0]['y']+0.01, tabPt[0]['y']],
                                        z=[zInf,zInf,zInf,zInf,zSup,zSup,zSup,zSup],
                                        i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2],
                                        j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3],
                                        k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6],
                                        color='lightpink',contour=dict(color='black',show=True,width=3)))
            else:    
                fig.add_trace(go.Mesh3d(x=[tabPt[k]['x'], tabPt[k]['x']+0.01, tabPt[k+1]['x']+0.01, tabPt[k+1]['x'], tabPt[k]['x'], tabPt[k]['x']+0.01, tabPt[k+1]['x']+0.01, tabPt[k+1]['x']],
                                        y=[tabPt[k]['y'], tabPt[k]['y']+0.01, tabPt[k+1]['y']+0.01, tabPt[k+1]['y'], tabPt[k]['y'], tabPt[k]['y']+0.01, tabPt[k+1]['y']+0.01, tabPt[k+1]['y']],
                                        z=[zInf,zInf,zInf,zInf,zSup,zSup,zSup,zSup],
                                        i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2],
                                        j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3],
                                        k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6],
                                        color='lightpink',contour=dict(color='black',show=True,width=3)))
            k += 1
        #'''
    except:
        print("Problème")

    return fig

@callback(
    #Output('{}-graph2D'.format(name_app),'figure'),
    Output('{}-calcul-button'.format(name_app), 'children'),
    Input('{}-calcul-button'.format(name_app), 'n_clicks'),
    State('{}-sout-z-inf'.format(name_app),'value'),
    State('{}-sout-z-sup'.format(name_app),'value'),
    State('{}-sout-lv'.format(name_app),'value'),
    State('{}-sout-lh'.format(name_app),'value'),
    State('{}-sout-E'.format(name_app),'value'),
    State('{}-sout-nu'.format(name_app),'value'),
    State('{}-sout-t'.format(name_app),'value'),
    State('{}-sout-tab'.format(name_app),'data'),
    prevent_initial_call=True
)
def calculate(n_click,zInf,zSup,lv,lh,E,nu,t,tabPt):
    E = E*1e6
    #H = zSup-zInf
    lst_coins = [(tabPt[k]['x'],tabPt[k]['y']) for k in range(len(tabPt))]
    recalculate = True

    # Traitement
    # ----------

    fileIn = open('input_batardeau.txt','w')

    fileIn.write(str(lh)+"\n")
    fileIn.write(str(lv)+"\n")
    fileIn.write(str(zInf)+"\n")
    fileIn.write(str(zSup)+"\n")
    fileIn.write(str(E)+"\n")
    fileIn.write(str(nu)+"\n")
    fileIn.write(str(t)+"\n")
    fileIn.write(str(len(lst_coins))+"\n")
    for k in range(len(lst_coins)):
        fileIn.write(str(lst_coins[k][0])+"\n")
        fileIn.write(str(lst_coins[k][1])+"\n")

    fileIn.close()
    #'''
    if recalculate:
        try:
            os.remove('output_batardeau.csv')
        except:
            "Rien"
        t0 = time.time()
        a = os.system(os.getcwd()+r'/batardeau_engine')
        temps = str(round(time.time()-t0,3))
        print('[s] '+temps)
    return "Calculer (dernier calcul: "+temps+" s)"

@callback(
    Output('{}-graph3D'.format(name_app),'figure'),
    Input('{}-calcul-button'.format(name_app), 'children'),
    Input('{}-post-defo'.format(name_app),'value'),
    Input('{}-post-colorize'.format(name_app),'value'),
    Input('{}-scaleFactor'.format(name_app),'value'),
    prevent_initial_call=True
)
def returnGraph3D(calcb, defo, colori, fact):

    fig = go.Figure()

    if colori in ['ux','uy','uz','|u|']:

        df = pd.read_csv('output_batardeau.csv')

        if defo=='Non':
            fact = 0
        #if defo=='Oui':
        #    max_coord = np.max(np.abs(df[['x','y']]))
        #    max_displ = np.max(np.abs(df[['ux','uy','uz']]))
        #    fact = 0.2*max_coord/max_displ
        
        if (colori=='|u|'):
            df['|u|'] = [(ux**2+uy**2+uz**2)**0.5 for ux,uy,uz in zip(df['ux'].tolist(),df['uy'].tolist(),df['uz'].tolist())]
        
        fig.add_trace(go.Scatter3d(x=[x+fact*ux for x,ux in zip(df['x'].tolist(),df['ux'].tolist())],
                                y=[y+fact*uy for y,uy in zip(df['y'].tolist(),df['uy'].tolist())],
                                z=[z+fact*uz for z,uz in zip(df['z'].tolist(),df['uz'].tolist())],
                                mode='markers',
                                marker=dict(
                                                #size=3,
                                                color=df[colori],                # set color to an array/list of desired values
                                                colorscale='plasma',   # choose a colorscale
                                                opacity=0.9,
                                                colorbar=dict(title=colori+' (m)')
                                ),
                                customdata=df[colori],
                                hovertemplate="x: %{x:.3f}m"+"<br>"+"y: %{y:.3f}m"+"<br>"+"z: %{z:.3f}mNGF"+"<br>"+colori+": %{customdata:,.3f}m"))
    else:
        df = pd.read_csv('output_batardeau_eff.csv')

        fig.add_trace(go.Scatter3d(x=df['x'].tolist(),
                                y=df['y'].tolist(),
                                z=df['z'].tolist(),
                                mode='markers',
                                marker=dict(
                                                #size=3,
                                                color=[m/1000 for m in df[colori].tolist()],                # set color to an array/list of desired values
                                                colorscale='plasma',   # choose a colorscale
                                                opacity=0.9,
                                                colorbar=dict(title=colori+' (kNm/m)')
                                            ),
                                customdata=df[colori],
                                hovertemplate="x: %{x:.3f}m"+"<br>"+"y: %{y:.3f}m"+"<br>"+"z: %{z:.3f}mNGF"+"<br>"+colori+": %{customdata:,.0f}Nm/m"))
    fig['layout']['scene']=dict(camera=dict(eye=dict(x=1.15, y=1.15, z=0.8)), #the default values are 1.25, 1.25, 1.25
                                xaxis=dict(),
                                yaxis=dict(),
                                zaxis=dict(),
                                aspectmode='data', #this string can be 'data', 'cube', 'auto', 'manual'
                                #a custom aspectratio is defined as follows:
                                #aspectratio=dict(x=1, y=1, z=1)
                                )
    return fig

if __name__ == '__main__':
    app.run(debug=True)