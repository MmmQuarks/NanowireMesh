import NanowireMesh as nwm
import plotly.graph_objects as go
import numpy as np

g = nwm.NanowireMesh(250, 250)

fig = go.Figure()

# add edge points
fig.add_trace(go.Scatter(
	x = [g.edges[edge]['x'] for edge in g.edges],
	y = [g.edges[edge]['y'] for edge in g.edges],
	mode = 'markers',
	marker = {'color' : 'blue'}
))

fig.show()