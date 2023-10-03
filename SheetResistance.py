import numpy as np
import sparse_linear_circuit_solver as circ
import pandas as pd
from scipy.interpolate import interp2d
import pdb

def _Correction(parLen, perpLen, probeSep):
	# in the language of smits from the bell tecnical manual:
	a = parLen #  length of sample along axis parallel to probes
	d = perpLen # length of sample along axis perpendicular to probes
	s = probeSep # separation between probes

	df = pd.read_csv('sheetResCorrectionFactors',
				index_col = 0)
	# columns of df are
	# d/s | a/d=1 | a/d=2 | a/3=3 | a/d>=4
	# first column shows d/s vals
	# other columns show the correction factor values for given a/d

	df = df.dropna(axis = 0) # remove any rows containing NaN 

	# we are taking data from a dataframe where the x, y, and z values are stored as shown below

	#	y axis label | x | x | x | x
	#	-----------------------------
	#	y            | z | z | z | z
	#	y            | z | z | z | z
	#	y            | z | z | z | z
	#	y            | z | z | z | z


	# x are the a/d values at which the correction factors are evaluated
	# y are the d/s values at which the correction factors are evaluated
	# z are the correction factors C(x,y)


	x = [float(str.replace(col,'a/d=','')) for col in df.columns]
	y = list(df.index)
	z = np.array(df) # use iloc so we don't get the index values in our np.array

	f = interp2d(x, y, z, bounds_error = True)

	try:
		return f(a/d, d/s)[0]
	except ValueError:
		pdb.set_trace()
		raise ValueError('Probe separation is too large for sample size')

def measure(g, probeSep, orientation = 'vertical'):

	# first measure with probe axis pointing from top electrode to bottom electrode
	if orientation == 'vertical':
		parLen = g.height
		perpLen = g.width
		theta = np.pi/2
	elif orientation == 'horizontal':
		parLen = g.width
		perpLen = g.height
		theta = 0
	else:
		raise ValueError("Orientation must be either 'vertical' or 'horizontal'")

	# make correction factor (and make sure that our dimensions are allowed)
	C = _Correction(parLen, perpLen, probeSep)

	# finding the points where the probe tips should go in the ideal case
	# should be evenly spaced about the center
	# if . is center and x is probe locations
	# x---x-.-x---x
	
	center = np.array( (g.width/2, g.height/2))

	# finding current injection node
	point = center + 1.5 * probeSep * np.array([np.cos(theta), np.sin(theta)])
	injectionNode = g.find_closest_node(*point)

	# finding current withdrawal node
	point = center - 1.5 * probeSep * np.array([np.cos(theta), np.sin(theta)])
	withdrawalNode = g.find_closest_node(*point)

	# finding first potential measurement node
	point = center + .5 * probeSep * np.array([np.cos(theta), np.sin(theta)])
	topVoltageNode = g.find_closest_node(*point)

	# finding second potential measurement node
	point = center - 0.5 * probeSep * np.array([np.cos(theta), np.sin(theta)])
	bottomVoltageNode = g.find_closest_node(*point)
	
	I = 1

	circ.solve_with_admittances(g,
				currentInjections = {injectionNode : I, withdrawalNode : -I},
				voltageName = 'voltage')

	V = g.nodes[topVoltageNode]['voltage'] - g.nodes[bottomVoltageNode]['voltage']
	return V / I * C
