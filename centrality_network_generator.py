import NanowireMesh as nwm
import numpy as np
import networkx as nx
import sys
import os

width = height = 100
nwLengthSD = 0
percMultiple = 2
nwDiam = .15
rcMean = 10
rcSD = 0

nwLengthVals = [3, 5, 7, 9]

for nwLength in nwLengthVals:
	outputFolder = sys.argv[1] 
	outputFolder = outputFolder + '/' +  '_'.join([str(width) + 'x' + str(height), 'perc' + str(percMultiple), 'nwLen' + str(nwLength), 'rc' + str(10)])
	if not os.path.exists(outputFolder):
		os.mkdir(outputFolder)
		print('Made folder to hold results')

	for n in range(100):
		g = nwm.NanowireMesh(width = width, height = height, nwLength = nwLength, nwLengthSD = nwLengthSD, percMultiple = percMultiple, nwDiam = nwDiam, nwDiamSD = 0, removeNonPercolating = True, rcMean = rcMean, rcSD = rcSD)
		g.to_pickle("/".join([outputFolder, 'net' + str(n).zfill(3)]) + '.p')

