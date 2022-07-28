### Basic static chloropleth map

from dash import Dash, html, dcc
import plotly.express as px
import pandas as pd
import json
import plotly.io as pio

# Show graph in a new browser
pio.renderers.default ="browser"

data = pd.read_csv("/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/Data/test chloro.csv")
geo_world = json.load(open('/Users/FM/Public/HDS/Summer Project/wastewater_r_estimation/Data/NHS_BUC.geojson', 'r'))

state_id_map = {}
for feature in geo_world['features']:
    feature['id'] = feature['properties']['nhser20cd']
    state_id_map[feature['properties']['nhser20nm']] = feature['id']

data['id'] = data['nhser20nm'].apply(lambda x: state_id_map[x])

# Create figure
fig = px.choropleth_mapbox(
    data,
    geojson=geo_world,
    locations='id',
    color='rest',
    hover_name='nhser20nm',
    hover_data=['rest'],
    mapbox_style="carto-positron",
    center={'lat':51, 'lon':0},
    zoom=5,
    opacity=0.8,
    color_continuous_scale=px.colors.diverging.Picnic,
    color_continuous_midpoint=1
)

fig.update_geos(fitbounds="locations", visible = False)

# Display figure
fig.show()