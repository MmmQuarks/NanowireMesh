import gnuplotlib as gp
import numpy as np
import pandas as pd
import NanowireMesh as nwm
import pdb

resultsFolder = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/Results/bc_paper_fom_calculation_01/'
assert resultsFolder[-1] == '/', 'results folder is missing a closing slash'
R0 = 10 # initial resistance in ohms
gamma = 1.29
startPercMultiple = 2
nwLength = 10E-6 # in meters
wavelength = 550E-9
mr = 0.055 + 3.32j # refractive index of silver at given wavelength
radius = 0.075E-6 # radius of nanowires in meters



df = pd.DataFrame(columns = ['percMultiple', 'resistance', 'transparency', 'fom'])

for percMultiple in np.linspace(startPercMultiple, 1.01, 100):
	resistance = R0 * ( (percMultiple - 1) / (startPercMultiple - 1) )**(-gamma)
	n_s = percMultiple * 5.63726 / nwLength**2
	transparency = nwm.transparency(wavelength = wavelength,
						mr = mr,
						radius = radius,
						n_s = n_s,
						nwLength = nwLength)
	fom = 188.5 / (resistance * (np.sqrt(1/transparency) - 1))
	df = df.append(dict(percMultiple = percMultiple,
				resistance = resistance,
				transparency = transparency,
				fom = fom),
				ignore_index = True)
print(df)

df2 = pd.DataFrame(columns = ['percMultiple', 'resistance', 'transparency', 'fom'])
for percMultiple in np.linspace(startPercMultiple, 1.01, 100):
	resistance = R0 * ( (percMultiple - 1) / (startPercMultiple - 1) )**(-gamma)
	n_s = percMultiple * 5.63726 / nwLength**2 / 2
	transparency = nwm.transparency(wavelength = wavelength,
						mr = mr,
						radius = radius,
						n_s = n_s,
						nwLength = nwLength)
	fom = 188.5 / (resistance * (np.sqrt(1/transparency) - 1))
	df2 = df2.append(dict(percMultiple = percMultiple,
				resistance = resistance,
				transparency = transparency,
				fom = fom),
				ignore_index = True)



theoryCurveOptions = {'with' : 'lines lc \"black\" lw 5',
			'legend' : 'theoretical'}
theoryCurve = (df['percMultiple'], df['fom'], theoryCurveOptions)
expCurveOptions = {'with' : 'lines lc \"red\" lw 5',
			'legend' : 'optimized'}
expCurve = (df2['percMultiple'], df2['fom'], expCurveOptions)
plotOptions = {'xrange' : '2.2:1',
		'hardcopy' : resultsFolder + 'figure_of_merit.png',
		'cmds' : ['set xlabel \"Density Relative to Percolation Threshold\" font \",20\" offset 0,-2',
				'set xtics font \",20\" offset 0,-1',
				'set ylabel \"FOM (1/Ohm)\" font \",20\"',
				'set ytics font \",20\" offset 0,0',
				'set bmargin 6',
				'set lmargin 10',
				'set rmargin 3',
				'set tmargin 5'] }
gp.plot(theoryCurve, expCurve, **plotOptions)
