import NanowireMesh as nwm
import pdb
import pandas as pd
import networkx as nx
import random
import gnuplotlib as gp
import numpy as np
from scipy import stats, interpolate, optimize
from copy import deepcopy
import ParallelCentrality as PC
import csv
import itertools
from sortedcontainers import SortedDict
import os
import sys
import seaborn as sns
from matplotlib import pyplot as plt

transparencyParams = dict(wavelength = 550E-9,
				mr = 0.055 + 3.32j)

# allowing for customization without changing this script or the submit script 
resultsFolder = sys.argv[1]
if resultsFolder[-1] != '/':
	resultsFolder = resultsFolder + '/'
runtimeOptionsDf = pd.read_csv(resultsFolder + 'runtimeOptions.csv')
# columns 
runtimeOptions = dict(zip(runtimeOptionsDf.iloc[:,0], runtimeOptionsDf.iloc[:,1]))
for key,val in runtimeOptions.items(): # converting strings into other data types if possible
	if val in ['True', 'False']:
		runtimeOptions[key] = val == 'True'
	else:
		try:
			runtimeOptions[key] = float(val)
		except ValueError:
			pass # this is not a number

srcNetworkFolder = runtimeOptions['srcNetworkFolder']

imgFormat = 'png'
# note this gets the effective perc multiple
def get_perc_multiple(g):
	assert g.width == g.height, "Width and height are not the same"
	lTot = sum([g.nodes[node]['length'] for node in g if node not in [g.topElectrode, g.bottomElectrode]])
	l = g.nwLength
	Ls = g.width
	return lTot * l / (5.63726 * Ls**2 + l * Ls + 5.5 * l**2)


def effective_perc_multiple_to_density(percMultiple, nwLength, height, width):
	assert height == width, "Height and width must be equal"
	nc_eff = 5.63726 / nwLength**2 + 1 / (nwLength * height) + 5.5 / (height**2)
	n_s = percMultiple * nc_eff
	return n_s

def fom(resistance, transparency):
	return 188.5 / (resistance * (transparency**(-1/2) - 1))


params = {'width' : None, 
		'height' : None,
		'nwLength' : None,
		'rcMean' : None,
		'rcSD' : None,
		'buffer' : 1,
		'percMultiple' : None,
		'nwDiam' : None,
		'buffer' : 2,
		'initialTemp' : 298.15,
		'addInternalResistance' : None}

# setting the params that are specific in runtimeOptions.csv
for key in params.keys():
	if params[key] == None:
		params[key] = runtimeOptions[key]

if 'evolveData.csv' in os.listdir(resultsFolder):
	evolveData = pd.read_csv(filepath_or_buffer = resultsFolder + 'evolveData.csv',
		index_col = False)

if 'initialData.csv' in os.listdir(resultsFolder):
	initialData = pd.read_csv(resultsFolder + 'initialData.csv')
else:
	initialData = []
	for fileName in os.listdir(resultsFolder + 'bc_paper_network_generator_001/'):
		if fileName[-2:] == '.p':
			# file names are like task_01_trial_05_initial.p
			# so if we split by '_' then the task number will be right after task in the list from the split
			fileNameSplit = fileName.split('_')
			task = int( fileNameSplit[ fileNameSplit.index('task') + 1] )
			trial = int( fileNameSplit[ fileNameSplit.index('trial') + 1] )
			row = {'task' : task,
				'trial' : trial}
			g = nwm.NanowireMesh(inPickle = resultsFolder + 'bc_paper_network_generator_001/' + fileName)
			row['resistance'] = g.sheetResistance
			row['percMultiple'] = get_perc_multiple(g)
			# here we make a row with the same date (task, trial, resistance, percMultiple) for each centrality.
			# this way when we merge this into evolveData we get the complete trajectory including the starting point
			# for each unique combo of task, trial, centrality
			for centrality in pd.unique(evolveData.centrality):
				row['centrality'] = centrality
				initialData.append(deepcopy(row))
	initialData = pd.DataFrame(initialData)
	initialData.to_csv(resultsFolder + 'initialData.csv',
				index = False)

# getting theoretical data predicted by percolation theory
if 'theoryData.csv' in os.listdir(resultsFolder):
	theoryData = pd.read_csv(resultsFolder + 'theoryData.csv')
else:
	group = evolveData.groupby(by = ['task', 'trial'])
	alpha0 = initialData.percMultiple.mean() # alpha0 should be the mean starting density
#	alpha0 = group.aggregate(func = 'max').percMultiple.mean()	 # finding mean starting percMultiple (mean of max pm)
#	R0 = group.aggregate(func = 'min').resistance.mean()	# finding mean starting resistance (mean of min res)
	R0 = initialData.resistance.mean() # R0 should be average starting resistance
	gamma = 1.29
	percMultipleArray = np.linspace(alpha0, 1.01, 150)
	theoryData = pd.DataFrame({'percMultiple' : percMultipleArray,
					'resistance' : [R0 * ( (pm - 1)/(alpha0 - 1) )**(-gamma) for pm in percMultipleArray]})
	theoryData.to_csv(resultsFolder + 'theoryData.csv',
				index = False)

	

# making density data 
f = lambda x : effective_perc_multiple_to_density(x.percMultiple, 
			nwLength = params['nwLength'] * 1E-6,  # because nw Len given in um in params
			height = params['height'] * 1E-6, # because height given in um in params
			width = params['width'] * 1E-6), # because width given in um in params
if 'density' not in evolveData.columns:
	print('Calculating simulated densities')
	evolveData['density'] = evolveData.apply(f, axis = 1)
	evolveData.to_csv(resultsFolder + 'evolveData.csv',
				index = False)
if 'density' not in theoryData.columns:
	print('Calculating theoretical densities')
	theoryData['density'] = theoryData.apply(f, axis = 1)
	theoryData.to_csv(resultsFolder + 'theoryData.csv',
				index = False)
if 'density' not in initialData.columns:
	print('Calculating initial densities')
	initialData['density'] = initialData.apply(f, axis = 1)
	initialData.to_csv(resultsFolder + 'initialData.csv',
				index = False)


# making transparency data
f = lambda x : nwm.transparency(radius = 0.15E-6 / 2,
				n_s = x.density,
				nwLength = params['nwLength'] * 1E-6,
				wavelength = transparencyParams['wavelength'],
				mr = transparencyParams['mr'])
if 'transparency' not in evolveData.columns:
	print('Calculating simulated transparencies')
	evolveData['transparency'] = evolveData.apply(f, axis = 1)
	evolveData.to_csv(resultsFolder + 'evolveData.csv',
				index = False)
if 'transparency' not in theoryData.columns:
	print('Calculating theoretical transparencies')
	theoryData['transparency'] = theoryData.apply(f, axis = 1)
	theoryData.to_csv(resultsFolder + 'theoryData.csv',
				index = False)
if 'transparency' not in initialData.columns:
	print('Calculating initial transparencies')
	initialData['transparency'] = initialData.apply(f, axis = 1)
	initialData.to_csv(resultsFolder + 'initialData.csv',
				index = False)

# making fom data
if 'fom' not in evolveData.columns: 
	print('Calculating simulated figures of merit')
	evolveData['fom'] = evolveData.apply(lambda x : fom(x.resistance, x.transparency), axis = 1)
	evolveData.to_csv(resultsFolder + 'evolveData.csv',
				index = False)
if 'fom' not in theoryData.columns:
	print('Calculating theoretical figures of merit')
	theoryData['fom'] = theoryData.apply(lambda x : fom(x.resistance, x.transparency), axis = 1)
	theoryData.to_csv(resultsFolder + 'theoryData.csv',
				index = False)
if 'fom' not in initialData.columns:
	print('Calculating initial figures of merit')
	initialData['fom'] = initialData.apply(lambda x : fom(x.resistance, x.transparency), axis = 1)
	initialData.to_csv(resultsFolder + 'initialData.csv',
				index = False)

# merging initial data with evolve data
evolveData = pd.concat([evolveData, initialData], ignore_index = True)



summaryData = pd.DataFrame(columns = ['centrality', 'doubling density', 'doubling density err', 'failure density', 'failure density err'])

trialStats = []
resInterps = {key : [] for key in pd.unique(evolveData.centrality)}
transInterps = {key : [] for key in pd.unique(evolveData.centrality)}
fomInterps = {key : [] for key in pd.unique(evolveData.centrality)}


for task in evolveData.task.unique():
	taskData = evolveData[evolveData.task == task]
	for centrality in taskData.centrality.unique():
		centralityData = taskData[taskData.centrality == centrality]
		for trial in centralityData.trial.unique():
			trialData = centralityData[centralityData.trial == trial]
			row = {'task' : task,
				'trial' : trial,
				'centrality' : centrality,
				'failure density' : trialData.percMultiple.min(),
				'start density' : trialData.percMultiple.max()}
			# calculating interpolants to calculate doubling density and for use later
			resInterp = interpolate.interp1d(trialData['percMultiple'],
							trialData['resistance'],
							#fill_value = 'extrapolate',
							bounds_error = False)
			transInterp = interpolate.interp1d(trialData['percMultiple'],
							trialData['transparency'],
							#fill_value = 'extrapolate',
							bounds_error = False)
			fomInterp = interpolate.interp1d(trialData['percMultiple'],
							trialData['fom'],
							#fill_value = 'extrapolate',
							bounds_error = False)

			# appending functions to lists
			resInterps[centrality].append(resInterp)
			transInterps[centrality].append(transInterp)
			fomInterps[centrality].append(fomInterp)
			dblFunc = lambda x : resInterp(x) - 2 * trialData['resistance'].min()
			rootResults = optimize.root_scalar(dblFunc,
								bracket = (trialData['percMultiple'].min(), trialData['percMultiple'].max())
							)
			row['doubling density'] = rootResults.root
			trialStats.append(row)
trialStats = pd.DataFrame(trialStats)

summaryData = trialStats.groupby('centrality').aggregate(func = ['mean','sem'])
summaryData.to_csv(resultsFolder + 'summaryData.csv', index = True)

startDensity = {centrality: trialStats[trialStats['centrality'] == centrality]['start density'].min() for centrality in pd.unique(trialStats.centrality)}
stopDensities = {centrality: trialStats[trialStats['centrality'] == centrality]['failure density'].max() for centrality in pd.unique(trialStats.centrality)}

# now we clean the evolve data to exclude all densities after the first networkfails
rowsToKeep = evolveData.apply(lambda x : x.percMultiple >= stopDensities[x.centrality] , axis = 1)
evolveData = evolveData[rowsToKeep]

#using interpolating functions to make our plot data
evolvePlotData = []
for centrality in pd.unique(evolveData.centrality):
	for pm in np.linspace(startDensity[centrality], stopDensities[centrality], 100):
		row = {'centrality' : centrality,
				'percMultiple': pm}
		resList = np.squeeze([f(pm) for f in resInterps[centrality]] )  # squeeze turns resList from [array(blah) ...] into [blah ...]		
		resList = [el for el in resList if el == el] # filtering out nan
		row['resistance'] = np.mean(resList)
		row['resistanceSEM'] = stats.sem(resList)

		transList = np.squeeze([f(pm) for f in transInterps[centrality] ])
		transList = [el for el in transList if el == el] 
		row['transparency'] = np.mean(transList)
		row['resistanceSEM'] = stats.sem(transList)
		
		fomList = np.squeeze([f(pm) for f in fomInterps[centrality]])
		fomList = [el for el in fomList if el == el]
		row['fom'] = np.mean(fomList)
		row['fomSEM'] = stats.sem(fomList)

		evolvePlotData.append(row)
evolvePlotData = pd.DataFrame(evolvePlotData)


# making theoretical optimal performance data
theoreticalOptDf = pd.DataFrame(columns = ['percMultiple', 'resistance', 'fom'])

# simplifying notation
L = params['height'] # in um
l = params['nwLength'] # in um
Rc = params['rcMean']

print('Generating dummy network to get access to methods')
g = nwm.NanowireMesh(width = 20, height = 20, nwDiam = params['nwDiam'], initialTemp = params['initialTemp'])
rho = g.calc_internal_resistivity(diam = params['nwDiam'], temp = params['initialTemp'])
Rw = rho * l / (np.pi * (params['nwDiam']/2)**2)

# resistance of single series of wires spanning the network
Rspan = (L/l + 1) * Rc + L/l * Rw

for alpha in np.linspace(evolveData.percMultiple.max() * 1.01,evolveData.percMultiple.min() * 0.99, 100):
	assert params['height'] == params['width']


	#density of wires in network in wires per m^2
	n_s = effective_perc_multiple_to_density(percMultiple = alpha,
							nwLength = l * 10**(-6),
							height = L * 10**(-6),
							width = L * 10**(-6))
	# number of wires in network including converting L into meters
	Nw = n_s * (L * 1E-6)**2

	# number of wire series spanning the whole network
	Nspan = Nw / (L/l)

	# resistance of network
	Rtot = 2 * Rspan / Nspan
#	Rtot = 2 * ( (L/l + 1) * Rc + L/l * Rw) / (n_s * L * l)

	# transparency calculation
	# note that we multiply n_s by two to account for the horizontal wires
	T = nwm.transparency(radius = params['nwDiam'] / 2 * 1E-6,
				n_s = n_s,
				nwLength = l * 1E-6,
				wavelength = transparencyParams['wavelength'],
				mr = transparencyParams['mr'])

	fomVal = fom(resistance = Rtot,
			transparency = T)
	theoreticalOptDf = theoreticalOptDf.append({'percMultiple' : alpha,
							'resistance' : Rtot,
							'fom' : fomVal},
							ignore_index = True)


theoreticalOptInterp  = {'resistance': interpolate.interp1d(theoreticalOptDf['percMultiple'], theoreticalOptDf['resistance']),
				'fom' : interpolate.interp1d(theoreticalOptDf['percMultiple'], theoreticalOptDf['fom'])}



# merging theoryData, and evolvePlotData
theoryData['centrality'] = 'Percolation Theory'
evolvePlotData = pd.concat([evolvePlotData, theoryData], ignore_index = True)


evolvePlotData['Scaled Resistance'] = evolvePlotData.apply(lambda x : x.resistance / theoreticalOptInterp['resistance'](x.percMultiple), axis = 1)
evolvePlotData['Scaled FOM'] = evolvePlotData.apply(lambda x : x.fom / theoreticalOptInterp['fom'](x.percMultiple), axis = 1)

# removing random plot data
evolvePlotData = evolvePlotData[(evolvePlotData.centrality != 'random') & (evolvePlotData.centrality != 'random_cleaned')]

# setting text settings for matplotlib
plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.sans-serif' : ['Computer Modern Roman'],
	'text.latex.preamble' : r'\usepackage{amsmath}'
	})
legendLabels = {
	'betweenness_centrality' : r"$\displaystyle C_B$",
	'unweighted_betweenness_centrality' : r"$\displaystyle C_{UWB}$",
	'electrode_centrality' : r"$\displaystyle C_\text{elec}$",
	'percolation_centrality' : r"$\displaystyle C_\text{perc}$",
	'current_weighted_centrality' : r"$\displaystyle C_\text{curr}$",
	'power_weighted_centrality' : r"$\displaystyle C_\text{pow}$",
	'Percolation Theory' : 'Percolation Theory'
	}


# plotting adapted from https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
# making plot grid
fig = plt.figure(
		figsize = (6.75, 7), # width, height in inches
		)
gridspec = fig.add_gridspec(len(pd.unique(evolvePlotData.centrality)), 2, hspace = 0, wspace = 0.4)
axs = gridspec.subplots(sharex = True)
fig.supxlabel('Density Relative To Effective Percolation Threshold', 
	fontsize = 'medium',
	y = 0.04)
#fig.supylabel('Resistance Relative to Optimal Network at Given Density')

colorIter = iter(plt.cm.brg(np.linspace(0,1, len(axs)+1)))
colorsDict = {centrality : next(colorIter) for centrality in pd.unique(evolvePlotData.centrality)}

# make percolation theory black
colorsDict['Percolation Theory'] = 'black'
lineStyleDict = {centrality : 'solid' if centrality != 'Percolation Theory' else 'dotted' for centrality in pd.unique(evolvePlotData.centrality)}

## making all columns in title case
#def title_case(name):
#	nameList = name.split(' ')
#	capsNameList = [el.capitalize() for el in nameList]
#	return ' '.join(capsNameList)
#evolvePlotData.rename(columns = {col : title_case(col) for col in evolvePlotData.columns},
#			inplace = True)

evolvePlotData['Relative Resistance'] = evolvePlotData['resistance'] / initialData.resistance.mean()
#def parse_centrality_names(x):
#	name = x.centrality.split('_')
#	name = [el.capitalize() for el in name]
#	name = ' '.join(name)
#	return name
#evolvePlotData['centrality'] = evolvePlotData.apply(parse_centrality_names, axis = 1)


for col in range(2):
	plotType = 'Relative Resistance' if col == 0 else 'Scaled FOM'
	for row in range(len(axs)):
		centrality = evolvePlotData.centrality.unique()[row]
		ax = axs[row, col]
		df = evolvePlotData[evolvePlotData.centrality == centrality]
		# making the background dim plots
		for backgroundCentrality in evolvePlotData.centrality.unique():
			if backgroundCentrality == centrality:
				pass # don't plot in background if this is the centrality we want to highlight
			else:
				bdf = evolvePlotData[evolvePlotData.centrality == backgroundCentrality]
				ax.plot(bdf['percMultiple'], bdf[plotType], 
					c = 'silver' if backgroundCentrality != 'Percolation Theory' else 'black',
					linestyle = lineStyleDict[backgroundCentrality],
					linewidth = 2 if backgroundCentrality == 'Percolation Theory' else 1)
		# plotting the resistance relative to start and the FOM relative to optimal network
		ax.plot(df['percMultiple'], df[plotType], 
				c = colorsDict[centrality],
				linestyle = lineStyleDict[centrality],
				label = centrality,
				linewidth = 2)
		handles, labels = ax.get_legend_handles_labels()
		labels = [legendLabels[el] for el in labels]
		leg = ax.legend(handles,
				labels,
				loc = 'upper right',
				fontsize = 'small',
				frameon = False,
				borderpad = 0.1)
		if row == len(axs) / 2:
			if col == 0:
				ax.set_ylabel('Resistance Relative to Start')
			else:
				ax.set_ylabel('FOM Relative to Optimal Network')
		if row == 0:
			if col == 0:
				title = 'Resistance Evolution'
			else:
				title = 'FOM Evolution'
			ax.set_title(title, 
				fontsize = 'medium',
				y = 1.1)
#		if col == 1:
#			ax.yaxis.set_label_position('right')
#			ax.yaxis.tick_right()
		ax.set_xlim(2.3, 0)
		if col == 0:
			ax.set_ylim(0,15)
		if centrality != 'Percolation Theory':
			nonpercDf = evolvePlotData[evolvePlotData.centrality != 'Percolation Theory']
			ymin = nonpercDf[plotType].min() * 0.9
			ymax = nonpercDf[plotType].max() * (1.1 if col == 1 else 1.35)

		else:
			ymin = evolvePlotData[plotType].min() * 0.9
			ymax = evolvePlotData[plotType].max() * 1.3
		ax.set_ylim(ymin, ymax)
		

plt.savefig(resultsFolder + 'evolution.pdf',
		pad_inches = 0)
plt.close()


# saving trial stats to csv
trialStats.to_csv(resultsFolder + 'trialStats.csv', index = False)


#reshaping the data to be easy for seaborn to handle
violinData = trialStats.melt(id_vars = 'centrality', value_vars = ['failure density', 'doubling density'],
	value_name = 'density',
	var_name = 'type')

violinData = violinData[(violinData.centrality != 'random') & (violinData.centrality != 'random_cleaned')] # filtering out random





violinData['centrality'] = violinData.apply(lambda x : legendLabels[x.centrality], axis = 1)
#violinData['centrality'] = violinData.apply(lambda x : x.centrality.replace(' ','\n'), axis = 1)
violinData['type'] = violinData.apply(lambda x : x['type'].capitalize().replace(' density', ''), axis = 1)

fig = plt.figure(figsize = (3.375, 3.75))
ax = sns.violinplot(data = violinData, 
			y = 'density', 
			x = 'centrality',
			hue = 'type', 
			split = True,
			inner = 'quartile',
			linewidth = 1,
			palette = {'Doubling' : 'white', 'Failure': 'lightgray'})
#ax.set_xlim(1.5,0)
ax.set_xlabel('Centrality', labelpad = 10,
		fontsize = 'medium')
ax.set_title('Failure and Doubling Distributions',
		fontsize = 'medium',
		pad = 10)
ax.set_ylabel('Density Relative to Effective Percolation Threshold', 
	fontsize = 'medium',
	labelpad = 10)
handles, labels = ax.get_legend_handles_labels()
# inverting the order that the labels appear
handles.reverse()
labels.reverse()
ax.legend(handles, labels, title = None,
		fontsize = 'small')
fig.tight_layout()

plt.savefig(resultsFolder + 'trialStats.pdf', pad_inches = 0)

evolvePlotData.to_csv(resultsFolder + 'evolvePlotData.csv', index = False)

pdb.set_trace()
