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


centralities = dict(rand = dict(name = 'random',
			color = 'black',
			func = lambda g : {node : np.random.rand() for node in g},
			clean = False,
			pointtype = 4,
			),
		randCleaned = dict(name = 'random_cleaned',
			color = 'black',
			func = lambda g : {node : np.random.rand() for node in g},
			clean = True,
			pointtype = 6
			),
		lowBC = dict(name = 'betweenness_centrality',
				color = 'dark-violet',
				func = lambda g : nx.betweenness_centrality(g, weight = 'resistance'),
				clean = True,
				pointtype = 8
				),
		lowEC = dict(name = 'electrode_centrality', 
				color = 'dark-green', 
				func = lambda g : nwm.electrode_centrality(g, potential = 'voltage', weight = 'resistance'),
				clean = True,
				pointtype = 12
				),
		lowPC = dict(name = 'percolation_centrality',
				color = 'red', 
				func = lambda g : nx.percolation_centrality(g, attribute = 'voltage', weight = 'resistance'),
				clean = True,
				pointtype = 17
				),
		lowCW = dict(name = 'current_weighted_centrality', 
				color = 'dark-orange', 
				func = lambda g : PC.current_weighted_centrality(g),
				clean = True,
				pointtype = 32
				),
		lowPW = dict(name = 'power_weighted_centrality',
				color = 'dark-blue', 
				func = lambda g : PC.power_weighted_centrality(g),
				clean = True,
				pointtype = 10
				),
		eig = dict(name = 'eigenvector_centrality',
				color = 'gray50',
				func = lambda g : nx.eigenvector_centrality(g, weight = 'admittance'),
				clean = True,
				pointtype = 90
				)
			)

if not runtimeOptions['useEigenvectorCentrality']:
	del centralities['eig']


evolveData = pd.read_csv(filepath_or_buffer = resultsFolder + 'evolveData.csv')
#except FileNotFoundError:
#	# to make evolveData csv if simulation terminated early (no longer necessary since code changed to write to csv intermittently rather than just at the end of sim)
#	files = os.listdir(resultsFolder)
#	files = sorted(files)
#	evolveData = pd.DataFrame(columns = ['trial', 'centrality', 'percMultiple', 'resistance'])
#	for fileName in files:
#		if '100x100' in fileName and 'initial' not in fileName:
#			termsFound = dict.fromkeys(evolveData.columns, False) # to make sure we only add data to the data frame if they have 
#			# all the values populated
#			row = dict.fromkeys(evolveData.columns)
#			splitName = fileName.split('_')
#			# getting trial number 
#			for term in splitName:
#				if 'trial' in term:
#					val = term.replace('trial','')
#					val = int(val)
#					termsFound['trial'] = True
#					row['trial'] = val
#			# finding centrality type
#			for key in centralities.keys():
#				cName = centralities[key]['name']
#				if cName in fileName:
#					termsFound['centrality'] = True
#					row['centrality'] = cName
#			g = nwm.NanowireMesh(inPickle = resultsFolder + fileName)
#			row['resistance'] = g.sheetResistance # getting sheet resistance
#			termsFound['resistance'] = True
#			row['percMultiple'] = get_perc_multiple(g)
#			termsFound['percMultiple'] = True
#			if all(termsFound.values()):
#				evolveData = evolveData.append(row, ignore_index = True)
#
#	evolveData.to_csv(path_or_buf = resultsFolder + 'evolveData.csv')
#
#				
#if makeGenData:
#	# making generated data
#	genData = pd.DataFrame(columns = ['percMultiple', 'resistance'])
#	percMultipleList = list(np.linspace(1.0, params['percMultiple'], 100))
#	percMultipleList = percMultipleList + 2 * [pm for pm in percMultipleList if pm < 1.4]
#	for percMultiple in percMultipleList:
#		for repetitions in range(2):
#			params['percMultiple'] = percMultiple
#			g = nwm.NanowireMesh(**params)
#			contactEdges = [edge for edge in g.edges if g.edges[edge]['resistanceType'] == 'cont']
#			if runtimeOptions['rcBellew']:
#				contactResistances = gen_bellew_resistances(k = len(contactEdges),
#									exclude_outliers = runtimeOptions['excludeOutliers'])
#
#			elif runtimeOptions['rcLognormal']:
#				mm = g.rcMean
#				sd = g.rcSD
#				mu = np.log(mm**2 / np.sqrt(mm**2 + sd**2))
#				sigma = np.sqrt(2 * np.log(np.sqrt(mm**2 + sd**2)/mm))
#				contactResistances = np.random.lognormal(mean = mu, sigma = sd, size = len(contactEdges))
#
#			nx.set_edge_attributes(g, 
#					dict(zip(contactEdges, contactResistances)),
#					'resistance')
#
#			
#			# really the comparison data should be uncleaned
#			# cleaning the network if requested
##			if simOptions['cleanNetwork'] == True:
##				clean_network(g)
#
#			g.solve_circuit_using_xyce(xycePath = xycePath,
#							netlistName = resultsFolder + 'netlist')
#			genData = genData.append(dict(percMultiple = get_perc_multiple(g),
#						resistance = g.sheetResistance),
#						ignore_index = True)
#	
#	genData.to_csv(path_or_buf = resultsFolder + 'genData.csv')

# making figure of merit data


performanceDf = pd.DataFrame(columns = ['centrality', 'doubling density', 'doubling density err', 'failure density', 'failure density err'])
allEvolveData = pd.DataFrame(columns = ['centrality',
					'percMultiple',
					'resistanceMean',
					'resistanceSem',
					'transparencyMean',
					'transparencySem',
					'fomMean',
					'fomSem'])
for key in centralities.keys():
	evolveData = pd.read_csv(resultsFolder + 'evolveData.csv')
	evolveData['taskTrial'] = evolveData.task.apply(func = str) + '_' + evolveData.task.apply(func = str)
	#genData = pd.read_csv(resultsFolder + 'genData.csv')
	
	# filtering genData to not go so near the percolation threshold
	#genData = genData[ genData['percMultiple'] > 1.1]
	
	# filtering evolveData to only have the data from the relevant removal strategy
	evolveData = evolveData[ evolveData['centrality'] == centralities[key]['name']]

	# making binned statistics
#	evolveResMeans, evolveResMeanBinEdges, evolveResMeanBinNumbers = stats.binned_statistic(evolveData['percMultiple'],
#													evolveData['resistance'],
#													statistic = 'mean',
#													bins = 10)
#	evolveResSems, evolveResSemBinEdges, evolveResSemBinNumbers = stats.binned_statistic(evolveData['percMultiple'],
#													evolveData['resistance'],
#													statistic = stats.sem,
#													bins = 10)
#	genResMeans, genResMeanBinEdges, genResMeanBinNumbers = stats.binned_statistic(genData['percMultiple'],
#													genData['resistance'],
#													statistic = 'mean',
#													bins = 15)
#	genResSems, genResSemBinEdges, genResSemBinNumbers = stats.binned_statistic(genData['percMultiple'],
#													genData['resistance'],
#													statistic = stats.sem,
#													bins = 15)
	# making interpolating functions for each method
	fInterp = {}
	failureDensity = []
	doublingDensity = []
	for trial in pd.unique(evolveData.taskTrial):
		trialData = evolveData[ evolveData['taskTrial'] == trial] # really this is the taskTrial rather than the trial but 
										# keeping this named trial saves rewriting a bunch of code
		fInterp[trial] = interpolate.interp1d(trialData['percMultiple'],
							trialData['resistance'],
							bounds_error = False)
		# getting average failure density
		failureDensity.append(trialData['percMultiple'].min())
		# getting average doubling density 
		dblFunc = lambda x : fInterp[trial](x) - 2 * trialData['resistance'].min()
		rootResults = optimize.root_scalar(dblFunc,
						bracket = (trialData['percMultiple'].min(), trialData['percMultiple'].max()))
		doublingDensity.append(rootResults.root)
			
	performanceDf = performanceDf.append({'centrality' : centralities[key]['name'],
						'doubling density' : np.average(doublingDensity),
						'doubling density err': stats.sem(doublingDensity),
						'failure density' : np.average(failureDensity),
						'failure density err' : stats.sem(failureDensity)},
						ignore_index = True)
	# formatting the evolution data
	evolveCurveOptions = {'with' : 'yerrorbars linestyle 7 lw 6 lc \"black\"',
				'tuplesize' : 3}
#	evolveCurvePercMultiples = [np.average([evolveResMeanBinEdges[n], evolveResMeanBinEdges[n+1]]) for n in range(len(evolveResMeanBinEdges) - 1)]
#	evolveCurveDf = pd.DataFrame(dict(percMultiple = evolveCurvePercMultiples,
#					resistanceMean = evolveResMeans,
#					resistanceSem = evolveResSems))



	evolveCurveDf = pd.DataFrame(columns = ['percMultiple', 
							'resistanceMean', 
							'resistanceSem', 
							'transparencyMean', 
							'transparencySem',
							'fomMean',
							'fomSem'])
	for count, percMultiple in enumerate(np.linspace(2.2,0, 40)):
		resArray = [np.float64(fInterp[trial](percMultiple)) for trial in pd.unique(evolveData.taskTrial)]
		# filtering out nan
		usefulResArray = [el for el in resArray if el == el]
		if len(usefulResArray) >= 10 or (count <= 5 and len(usefulResArray) > 0):
			# calculating density from effective percolation multiple
			# using interpolation we measure all networks at the same density so their
			# transparencies will all be the same
			resistanceMean = np.average(usefulResArray)
			resistanceSem = stats.sem(usefulResArray)
			n_s = effective_perc_multiple_to_density(percMultiple = percMultiple,
									nwLength = params['nwLength']*1E-6,
									height = params['height']*1E-6,
									width = params['width']*1E-6)
			transparency = nwm.transparency(radius = 0.15E-6 / 2,
								n_s = n_s,
								nwLength = params['nwLength']*1E-6, # multiplied by 1E-6 bc nwLength given in um
								wavelength = transparencyParams['wavelength'],
								mr = transparencyParams['mr'])
			transparencySem = 0
			fomList = [fom(resistance = resistance, transparency = transparency) for resistance in usefulResArray]
			fomMean = np.average(fomList)
			fomSem = stats.sem(fomList)
			evolveCurveDf = evolveCurveDf.append({'percMultiple' : percMultiple,
								'resistanceMean' : np.average(usefulResArray),
								'resistanceSem' : stats.sem(usefulResArray),
								'transparencyMean' : transparency,
								'transparencySem' : transparencySem,
								'fomMean' : fomMean,
								'fomSem' : fomSem},
								ignore_index = True)

	# appending data to allEvolveData
	evolveCurveDf['centrality'] = centralities[key]['name'].replace('_',' ')
	allEvolveData = allEvolveData.append(evolveCurveDf, ignore_index = True)
	
	evolveCurve = (evolveCurveDf['percMultiple'], 
				evolveCurveDf['resistanceMean'],
				evolveCurveDf['resistanceSem'],
				evolveCurveOptions)
	
	# formatting the gen data
	genCurveOptions = {'with' : 'filledcurves lc \"orange\"',
				'tuplesize' : 3}
	
#	genCurvePercMultiples = [np.average([genResMeanBinEdges[n], genResMeanBinEdges[n+1]]) for n in range(len(genResMeanBinEdges) - 1)]
#	genCurveDf = pd.DataFrame(dict(percMultiple = genCurvePercMultiples,
#				resistanceMean = genResMeans,
#				resistanceSem = genResSems))
#	genCurve = (genCurveDf['percMultiple'], 
#				genCurveDf['resistanceMean'] - genCurveDf['resistanceSem'],
#				genCurveDf['resistanceMean'] + genCurveDf['resistanceSem'],
#				genCurveOptions)
#	
#	genCurveMeanOptions = {'with' : 'lines lw 3 lc \"blue\"'}
#	genCurveMean = (genCurveDf['percMultiple'],
#				genCurveDf['resistanceMean'],
#				genCurveMeanOptions)

	# making theory data

	alpha0 = evolveCurveDf.loc[0]['percMultiple']
	R0 = evolveCurveDf.loc[0]['resistanceMean']
	gamma = 1.29
	percMultipleArray = np.linspace(2.2, 1.01, 150)
	theoryDf = pd.DataFrame({'percMultiple' : percMultipleArray,
					'resistance' : [R0 * ( (pm - 1)/(alpha0 - 1) )**(-gamma) for pm in percMultipleArray]})
	theoryCurveOptions = {'with' : 'lines lw 6 lc \"red\"'}
	theoryCurve = (theoryDf['percMultiple'],
				theoryDf['resistance'],
				theoryCurveOptions)
	
	# making fom data for theory curve FOM
	# we remake the theory data for each plot (because we are plotting only one centrality at a time here and I am stupid)
	theoryFOM = []
	for index, row in theoryDf.iterrows():
		n_s = effective_perc_multiple_to_density(percMultiple = row['percMultiple'],
								nwLength = params['nwLength']*1E-6, 
								height = params['height']*1E-6,
								width = params['width']*1E-6)
		# note the lengths given above are multiplied by 1E-6 because originally they
		# were given in microns
		transparency = nwm.transparency(radius = 0.15E-6 / 2,
							n_s = n_s,
							nwLength = params['nwLength']*1E-6,
							wavelength = transparencyParams['wavelength'],
							mr = transparencyParams['mr'])
		theoryFOM.append( fom(transparency = transparency,
					resistance = row['resistance']))
	theoryDf['fom'] = theoryFOM
	# adding the zero values for theory curve
#	theoryDf = theoryDf.append(dict(resistance = [np.inf, np.inf],
#						percMultiple = [1,0],
#						fom = [0,0]),
#					ignore_index = True)
	theoryFOMOptions = {'with' : 'lines lw 6 lc \"red\"'}
	theoryFOMCurve = (theoryDf['percMultiple'],
				theoryDf['fom'],
				theoryFOMOptions)

	# making fom data for evolve curve FOM
#	evolveFOM = []
#	for index, row in evolveCurveDf.iterrows():
#		evolveFOM.append( fom(resistance = row['resistanceMean'],
#					percMultiple = row['percMultiple']))
#	evolveCurveDf['fom'] = evolveFOM
#	evolveFOMOptions = {'with' : 'points linestyle 7 lw 6 lc \"black\"'}
	evolveFOMOptions = {'with' : 'yerrorbars linestyle 7 lw 6 lc \"black\"',
				'tuplesize' : 3}
	evolveFOMCurve = (evolveCurveDf['percMultiple'], 
				evolveCurveDf['fomMean'],
				evolveCurveDf['fomSem'],
				evolveFOMOptions)


	# aspect ratio is 4:3
	# resolution is 72dpi
	# we want triple that
	dpcm = int(72 * 2.54) # dots per centimeter
	imgWidth = int(8.6 * dpcm) # img should be 8.6 cm wide
	imgHeight = int(3/4 * imgWidth)
	fomPlotOptions = {'xrange' : '2.2:0.2',
		'yrange' : '0:150',
		'terminal' : 'png size ' + str(imgWidth) + ',' + str(int(0.5 * imgHeight)),
		'unset' : ['grid'],
		'set' : ['xtics font \",40\" offset 0,-2',
				'ytics font \",40\" offset 0,0 0,50,150',
				'bmargin 10',
				'lmargin 19',
				'tmargin 5',
				'xlabel \"Density Relative to Percolation Threshold\" font \",40\" offset 0,-5',
				'ylabel \"FOM [1/Ohms]\" font \",40\" offset -8,0'],
		'hardcopy' : resultsFolder + 'plot_' + centralities[key]['name'] + '_fom.png',
		'cmds' : []}
	gp.plot(theoryFOMCurve, evolveFOMCurve, **fomPlotOptions)

	
	# formatting plot options
	plotOptions = {'xrange' : '2.2:0.2',
			'yrange' : '0:350',
			'terminal' : 'png size ' + str(imgWidth) + ',' + str(int(0.5 * imgHeight)),
			'unset' : ['grid'],
			'set' : ['xtics font \",40\" offset 0,-2',
					'ytics font \",40\" offset 0,0 0,100,300',
					'bmargin 10',
					'lmargin 19',
					'tmargin 2',
					'xlabel \"Density Relative to Percolation Threshold\" font \",40\" offset 0,-5',
					'ylabel \"Resistance [Ohms]\" font \",40\" offset -8,-1'],
			'hardcopy' : resultsFolder + 'plot_' + centralities[key]['name'] + '_evolution.png',
			'cmds' : []}
	
	gp.plot(theoryCurve, evolveCurve, **plotOptions)

performanceDf.to_csv(path_or_buf = resultsFolder + 'summary.csv')


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

for alpha in np.linspace(2.2,0.2, 100):
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

# make all plots on same axes with error bands
for plotType in ['resistance', 'fom']:
	plotWidth = '1200'
	plotHeight =  '1000' #'1400' if plotType == 'resistance' else '900'
	ylabel = 'Resistance Relative to Start' if plotType == 'resistance' else 'FOM [1/Ohms]'
	xlabel = 'Density Relative to Percolation Threshold'
	yMax = '35' if plotType == 'resistance' else '115'
	yTicIncrement = '10' if plotType =='resistance' else '20'
	#keySetting = 'key inside left top font \",30\" reverse vertical Left samplen 1' if plotType == 'resistance' else 'key off'
	keySetting = ' '.join(['key inside left',
				'top' if plotType == 'resistance' else 'bottom',
				'font \",27\" reverse vertical Left samplen 1'])

	resScale = allEvolveData['resistanceMean'].min() if plotType == 'resistance' else 1 # factor to normalize all resistances to starting resistance
	resPlotOptions = {'xrange' : '2.1:0.3',
			'yrange' : '0:' + yMax,
#			'terminal' : 'png truecolor size ' + ','.join([plotWidth, plotHeight]),
			'set' : ['xtics axis font \",30\" offset 0,-2',
					'ytics border font \",30\" offset 0,-.7 ' +  ','.join(['0', yTicIncrement, yMax,]),
					'tmargin ' + ('0' if imgFormat == 'ps' else  '0.5'),
					'lmargin 13',
					'bmargin 10',
					'rmargin 0.5',
					keySetting,
					'xlabel \"' + xlabel + '\" font \",30\" offset 0,-5',
					'ylabel \"' + ylabel + '\" font \",30\" offset -3,-1',
					'style fill transparent solid 0.3'],
			'unset' : ['grid'],
#			'hardcopy' : resultsFolder + 'plot_all_' + plotType + '_bands_evolution.png',
			'cmds' : []}

	resPlotOptions['terminal'] = 'postscript eps size 10in,10in' if imgFormat == 'ps' else 'png truecolor size ' + ','.join([plotWidth, plotHeight])
	resPlotOptions['hardcopy'] = ''.join([resultsFolder,
						'plot_all_',
						plotType,
						'_bands_evolution.',
						'ps' if imgFormat == 'ps' else 'png'])

	theoryOptions = {'with' : 'lines lw 6 lc \"black\"',
				'legend' : 'Percolation Theory'}
	theoryCurve = (theoryDf['percMultiple'],
			theoryDf[plotType] / resScale,
			theoryOptions)

	allTuples = [theoryCurve] 
	
	# error bands
	for n, key in enumerate(centralities.keys()):
		name = centralities[key]['name'].replace('_',' ')
		if 'random' in name: # we are not including the randoms in the plots
			continue
		thisDf = allEvolveData[allEvolveData['centrality'] == name]

		curveOptions = {'with' : 'filledcurves lc \"' + centralities[key]['color'] + '\"',
					'tuplesize' : 3}
		resTuple = (thisDf['percMultiple'], 
				(thisDf[plotType + 'Mean'] - thisDf[plotType + 'Sem']) / resScale, 
				(thisDf[plotType + 'Mean'] + thisDf[plotType + 'Sem']) / resScale,
				curveOptions)
		allTuples.append(resTuple)
	# points 
	for n, key in enumerate(centralities.keys()):
		name = centralities[key]['name'].replace('_',' ')
		if 'random' in name: # we are not including the randoms in the plots
			continue
		thisDf = allEvolveData[allEvolveData['centrality'] == name]
		curveOptions = {'with' : 'points  lc \"' + centralities[key]['color'] + '\" lw 4 pointtype ' + str(centralities[key]['pointtype']) + ' pointsize 2.5',
					'legend' : ' '.join([w.capitalize() for w in name.split()])}
		resTuple = (thisDf['percMultiple'], 
				thisDf[plotType + 'Mean'] / resScale,
				curveOptions)
		allTuples.append(resTuple)
	
	gp.plot(*allTuples, **resPlotOptions)

# making the plots with the optimal values in there as well
for plotType in ['resistance', 'fom']:
	plotWidth = '1200'
	plotHeight =  '850' #'1400' if plotType == 'resistance' else '900'
	ylabel = ' '.join(['Resistance' if plotType == 'resistance' else 'FOM', 'Relative to Optimal'])
	xlabel = 'Density Relative to Percolation Threshold'
	yMax = '30' if plotType == 'resistance' else '0.37'
	yTicIncrement = '5' if plotType =='resistance' else '0.05'
	#keySetting = 'key inside left top font \",30\" reverse vertical Left samplen 1' if plotType == 'resistance' else 'key off'
	keySetting = ' '.join(['key inside left',
				'top' if plotType == 'resistance' else 'bottom',
				'font \",27\" reverse vertical Left samplen 1'])

	resScale = theoreticalOptDf[plotType].min()
	resPlotOptions = {'xrange' : '2.1:0.5',
			'yrange' : '0:' + yMax,
			'terminal' : 'png truecolor size ' + ','.join([plotWidth, plotHeight]),
			'set' : ['xtics axis font \",30\" offset 0,-2',
					'ytics border font \",30\" offset 0,-.7 ' +  ','.join(['0', yTicIncrement, yMax,]),
					'tmargin 3',
					'lmargin 15',
					'bmargin 10',
					'rmargin 0.5',
					keySetting,
					'xlabel \"' + xlabel + '\" font \",30\" offset 0,-5',
					'ylabel \"' + ylabel + '\" font \",30\" offset -4,-1',
					'style fill transparent solid 0.3'],
			'unset' : ['grid'],
			'hardcopy' : resultsFolder + 'plot_all_' + plotType + '_bands_evolution_with_optimal.png',
			'cmds' : []}

	theoryOptions = {'with' : 'lines lw 6 lc \"black\"',
				'legend' : 'Percolation Theory'}
	theoryCurve = (theoryDf['percMultiple'],
			theoryDf[plotType] / resScale,
			theoryOptions)

	allTuples = [theoryCurve] 
	theoreticalMinCurve = (theoreticalOptDf['percMultiple'],
				theoreticalOptDf[plotType] / resScale,
				{'with' : 'lines lw 6 lc \"red\"',
					'legend' : 'Theoretical Optimum'})
#	allTuples.append(theoreticalMinCurve)
	
	# error bands
	for n, key in enumerate(centralities.keys()):
		name = centralities[key]['name'].replace('_',' ')
		if 'random' in name: # we are not including the randoms in the plots
			continue
		thisDf = allEvolveData[allEvolveData['centrality'] == name]
		curveOptions = {'with' : 'filledcurves lc \"' + centralities[key]['color'] + '\"',
					'tuplesize' : 3}

		optimalDf = pd.DataFrame(columns = ['percMultiple', plotType])
		optimalDf['percMultiple'] = thisDf['percMultiple']
		optimalDf = optimalDf[optimalDf['percMultiple'] >= theoreticalOptDf['percMultiple'].min()]
		optimalDf[plotType] = np.array([theoreticalOptInterp[plotType](pm) for pm in optimalDf['percMultiple'] ])
		resTuple = (thisDf['percMultiple'], 
				(thisDf[plotType + 'Mean'] - thisDf[plotType + 'Sem']) / optimalDf[plotType], 
				(thisDf[plotType + 'Mean'] + thisDf[plotType + 'Sem']) / optimalDf[plotType],
				curveOptions)
		allTuples.append(resTuple)
	# points 
	for n, key in enumerate(centralities.keys()):
		name = centralities[key]['name'].replace('_',' ')
		if 'random' in name: # we are not including the randoms in the plots
			continue
		thisDf = allEvolveData[allEvolveData['centrality'] == name]
		optimalDf = pd.DataFrame(columns = ['percMultiple', plotType])
		optimalDf['percMultiple'] = thisDf['percMultiple']
		optimalDf = optimalDf[optimalDf['percMultiple'] >= theoreticalOptDf['percMultiple'].min()]
		optimalDf[plotType] = np.array([theoreticalOptInterp[plotType](pm) for pm in optimalDf['percMultiple']])
		curveOptions = {'with' : 'points  lc \"' + centralities[key]['color'] + '\" lw 4 pointtype ' + str(centralities[key]['pointtype']) + ' pointsize 2.5',
					'legend' : ' '.join([w.capitalize() for w in name.split()])}
		resTuple = (thisDf['percMultiple'], 
				thisDf[plotType + 'Mean'] / optimalDf[plotType],
				curveOptions)
		allTuples.append(resTuple)
	
	gp.plot(*allTuples, **resPlotOptions)

# trying to make multiplot	
#allTuplesDict = dict(resistance = [], fom = [])
#for plotType in ['resistance', 'fom']:
#	ylabel = 'Resistance [Ohms]' if plotType == 'resistance' else 'FOM [1/Ohms]'
#	yMax = '300' if plotType == 'resistance' else '180'
#	resPlotOptions = {'xrange' : '2.2:0.2',
#			'yrange' : '0:' + yMax,
#			'terminal' : 'png truecolor size 1200,900',
#			'multiplot' : 'layout 2,1',
#			'set' : ['xtics font \",40\" offset 0,-2',
#					'ytics font \",40\" offset 0,0 0,50,' + yMax,
#					'bmargin 10',
#					'lmargin 19',
#					'tmargin 2',
#					'key inside top left font \",25\" reverse Left',
#					'xlabel \"Density Relative to Percolation Threshold\" font \",40\" offset 0,-5',
#					'ylabel \"' + ylabel + '\" font \",40\" offset -8,-1',
#					'style fill transparent solid 0.3'],
#			'hardcopy' : resultsFolder + 'plot_all_' + plotType + '_bands_evolution.png',
#			'cmds' : []}
#
#	theoryOptions = {'with' : 'lines lw 6 lc \"black\"',
#				'legend' : 'Theory'}
#	theoryCurve = (theoryDf['percMultiple'],
#			theoryDf[plotType],
#			theoryOptions)
#
#	allTuplesDict[plotType].append(theoryCurve) 
#	
#	for n, key in enumerate(centralities.keys()):
#		name = centralities[key]['name'].replace('_',' ')
#		if 'random' in name: # we are not including the randoms in the plots
#			continue
#		thisDf = allEvolveData[allEvolveData['centrality'] == name]
#		curveOptions = {'with' : 'filledcurves lc \"' + centralities[key]['color'] + '\"',
#					'tuplesize' : 3}
#		resTuple = (thisDf['percMultiple'], 
#				thisDf[plotType + 'Mean'] - thisDf[plotType + 'Sem'], 
#				thisDf[plotType + 'Mean'] + thisDf[plotType + 'Sem'],
#				curveOptions)
#		allTuples.append(resTuple)
#	for n, key in enumerate(centralities.keys()):
#		name = centralities[key]['name'].replace('_',' ')
#		if 'random' in name: # we are not including the randoms in the plots
#			continue
#		thisDf = allEvolveData[allEvolveData['centrality'] == name]
#		curveOptions = {'with' : 'points  lc \"' + centralities[key]['color'] + '\" lw 4',
#					'legend' : ' '.join([w.capitalize() for w in name.split()])}
#		resTuple = (thisDf['percMultiple'], 
#				thisDf[plotType + 'Mean'],
#				curveOptions)
#		allTuples.append(resTuple)
#	
#gp.plot(*allTuplesDict['resistance'], *allTuplesDict['fom'],**resPlotOptions)
	
# writing to formatted table for latex
tString = [' & '.join(['centrality',
				'doubling density',
				'failure density \\\\']
			)]
centralityToTexMapping = {'random' : 'rand',
				'random_cleaned' : 'random cleaned',
				'betweenness_centrality' : '$C_B$',
				'electrode_centrality' : '$C_{\\text{elec}}$',
				'percolation_centrality' : '$C_{\\text{perc}}$',
				'current_weighted_centrality' : '$C_{\\text{curr}}$',
				'power_weighted_centrality' : '$C_{\\text{pow}}$',
				'eigenvector_centrality' : '$C_{\\text{eig}}$'}

for index, row in performanceDf.iterrows():
	cent = row['centrality']
	symbol = centralityToTexMapping[cent]
	tRow = [symbol]
	tRow.append( '$' + str(round(row['doubling density'],3)) + ' \\pm ' + str(round(row['doubling density err'],3)) + '$')
	tRow.append( '$' + str(round(row['failure density'],3)) + ' \\pm ' + str(round(row['failure density err'],3)) + '$ ' + '\\\\')
	tString.append(' & '.join(tRow))

for row in tString:
	print(row)

with open(resultsFolder + 'summary_formatted.txt', 'w') as f:
	f.write('\n'.join(tString))

# making snapshots of network to use in paper
#filesToImg = ['100x100_percMult2_nwLen10_Rc10_trial6_current_weighted_centrality_iter00000.p',
#			'100x100_percMult2_nwLen10_Rc10_trial6_current_weighted_centrality_iter00280.p',
#			'100x100_percMult2_nwLen10_Rc10_trial6_current_weighted_centrality_iter00520.p',
#			'100x100_percMult2_nwLen10_Rc10_trial6_current_weighted_centrality_iter00793.p']
#for fileName in filesToImg:
#	g = nwm.NanowireMesh(inPickle = resultsFolder + fileName)
#	g.to_img(outFile = resultsFolder + fileName.replace('.p',''),
#			title = '',
#			xLabel = '',
#			yLabel = '',
#			lineWidth = 10,
#			openImage = False,
#			showJunctions = True)
#	g.show_percolating(output = resultsFolder + fileName.replace('.p','_v2.png'))
