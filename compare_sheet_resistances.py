import NanowireMesh as nwm
import numpy as np
import os
import pickle
import time

size = 400
percolationMultiples = np.linspace(1.5, 4, 11)
percolationMultiples = list(percolationMultiples) * 6
percolationMultiples = sorted(percolationMultiples)
percolationMultiples = [1.5]
probeSeparation = 40 # in um - should vary this

# things to change
# 1. move around center of sample - done
# 2. vary density - done
# 3. change probe distances

data = {'trial' : [], 
	'size' : [], 
	'linearResistance' : [],
	'sheetResistance' : [],
	'sheetResistanceSD' : [],
	'probeSeparation' : [],
	'percolationMultiple' : [],
	'solveTimeForLinear' : []
	}

xycePath = 'mpirun -n 64 /home/amwt/TPV/2DNanowires/installations/build_Xyce-6.10_parallel/src/Xyce'

for trial, percolationMultiple in enumerate(percolationMultiples):
	g = nwm.NanowireMesh(width = size, height = size, percMultiple=  percolationMultiple)
	# data.update({size : {}})
	g.measure_sheet_resistance(probeSeparation, trials = 5, xycePath= xycePath)
	g.to_netlist()
	t0 = time.time()
	os.system(xycePath + ' netlist')
	solveTimeForLinear = time.time() - t0
	g.update_with_xyce_output('netlist.csv')
	# data[size].update({'linear' : g.sheetResistance, 'measured' : g.sheetResistanceMeasured, 'measuredSD' : g.sheetResistanceMeasuredSD})

	# write all info to data dictionary
	data['trial'].append(trial)
	data['size'].append(size)
	data['linearResistance'].append(g.sheetResistance)
	data['sheetResistance'].append(g.sheetResistanceMeasured)
	data['sheetResistanceSD'].append(g.sheetResistanceMeasuredSD)
	data['probeSeparation'].append(probeSeparation)
	data['percolationMultiple'].append(percolationMultiple)
	data['solveTimeForLinear'].append(solveTimeForLinear)


# write to file
pickle.dump(data, open('sheetResData.p', 'wb'))
