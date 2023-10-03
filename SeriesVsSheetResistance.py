import NanowireMesh as nwm
import pandas as pd
import SheetResistance as sr
import sparse_linear_circuit_solver as splcs
import numpy as np
import gnuplotlib as gp
import pickle
import pdb
import time
from os import system

columns = ['sampleSize', 'nwLen', 'percMultiple', 'R_series', 
		'R_sheet_vertical', 'R_sheet_horizontal', 'R_sheet_avg', 'probeSep',
		'num_nodes', 'xyce_solve_time', 'sp_solve_time']
sampleSize = 100
data = []

resultsFolder = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/Results/SeriesVsSheetResistance_01/'
makeData = False
makePlot = True

#simulation params
percMultipleSequence = np.linspace(1.5, 7, 4)
nwLengthSequence = [5, 10, 20]
if makeData:
	counter = 0
	for percMultiple in percMultipleSequence:
		for nwLength in nwLengthSequence:
			for probeSep in [30]:
				for trial in range(5):
					success = False
					print('percMultiple, nwLength, probeSep, trial', (percMultiple, nwLength, probeSep,trial))
					while not success:
						dataRow = dict.fromkeys(columns, None)
						g = nwm.NanowireMesh(width = sampleSize,
									height = sampleSize,
									percMultiple = percMultiple,
									nwLength = nwLength,
									rcMean = 'bellew')
						try:
							start = time.time()
							dataRow['R_sheet_vertical'] = sr.measure(g, probeSep, orientation = 'vertical')
							dataRow['R_sheet_horizontal'] = sr.measure(g, probeSep, orientation = 'horizontal')
							dataRow['sp_solve_time'] = (time.time() - start)/2
						except AssertionError:
							print('Retrying')
							time.sleep(1)
							continue
						dataRow['R_sheet_avg'] = (dataRow['R_sheet_vertical'] + dataRow['R_sheet_horizontal'])/2
						start = time.time()
						g.solve_circuit_using_xyce()
						dataRow['xyce_solve_time'] = time.time() - start
						dataRow.update(dict(sampleSize = sampleSize,
									nwLength = nwLength,
									percMultiple = percMultiple,
									R_series = g.sheetResistance,
									probeSep = probeSep,
									num_nodes = len(g)))
	
						data.append(dataRow)
						if counter % 10 == 0:
							pickle.dump(data, open(resultsFolder + 'partial.p', 'wb'))
						counter += 1
						success = True
	
	df = pd.DataFrame(data)
	df.to_csv(path_or_buf = resultsFolder + 'data.csv', index = False)
	
	pdb.set_trace()
else:
	df = pd.read_csv(resultsFolder + 'data.csv')

if makePlot:
	sf = df.groupby(by = ['nwLength', 'percMultiple']).aggregate(func = ['mean', 'sem'])
	plotTuples = []
	plotOptions = {'hardcopy' : resultsFolder + 'plot.png',
			'set' : ['key font \",20\"',
				'pointsize 2',
				'title \"Sheet Resistance and Series Resistance vs Density\" font \",20\"',
				'xlabel \"Multiple of Percolation Threshold\" font \",18\"',
				'ylabel \"Resistance or Sheet Resistance\" font \",18\"'],
			'unset' : ['grid']}
	thisDf = sf.loc[5]
	for col in ['R_series', 'R_sheet_avg']:
		curveOptions = {'with' : 'yerrorbars pointtype 7',
					'tuplesize' : 3,
					'legend' : col.replace('_avg','') }
		thisTuple = (thisDf.index,
				thisDf[col]['mean'],
				thisDf[col]['sem'],
				curveOptions)
		plotTuples.append(thisTuple)
	
	gp.plot(*plotTuples, **plotOptions)
	system('open ' + resultsFolder + 'plot.png')


	print('making second attempt at plot based on percentdifference')
	plotTuples = []
	df['% Difference'] = (df['R_series'] - df['R_sheet_avg']) / df['R_sheet_avg']
	sf = df.groupby(by = ['nwLength', 'percMultiple']).aggregate(func = ['mean', 'sem'])
	plotOptions['hardcopy'] = resultsFolder + 'plot2.png'
	plotOptions['set'][-1] = 'ylabel \"(R_series - R_sheet) / R_sheet\" font \",18\"'
	for nwLength in nwLengthSequence:
		thisDf = sf.loc[nwLength]
		curveOptions = {'with' : 'yerrorbars pointtype 7',
					'tuplesize' : 3,
					'legend' : str(nwLength) + ' um' }
		thisTuple = (thisDf.index,
				thisDf['% Difference']['mean'],
				thisDf['% Difference']['sem'],
				curveOptions)
		plotTuples.append(thisTuple)
	gp.plot(*plotTuples, **plotOptions)
	system('open ' + resultsFolder + 'plot2.png')

