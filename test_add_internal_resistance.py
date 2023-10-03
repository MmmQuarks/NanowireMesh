import pdb
import pandas as pd
import itertools
import networkx as nx
import numpy as np
from copy import deepcopy
import NanowireMesh as nwm

g = nwm.NanowireMesh(addInternalResistance = False,
			width = 100,
			height = 100)
h = deepcopy(g)
h.add_internal_resistance_4()

gNonElectrodeNodes = [node for node in g if node not in [g.topElectrode, g.bottomElectrode]]
hNonElectrodeNodes = [node for node in h if node not in [h.topElectrode, h.bottomElectrode]]

# g is network without internal resistance
# h is network with internal resistance

def ramp(x):
	if x > 0:
		return x
	else:
		return 0

# number of internal resistors should be SUM_i ( ramp(c_i - 1) ) where c_i is the 
# number of contacts a given NONELECTRODE wire has. This should be summed over all
# non electrode wires

results = pd.DataFrame(columns = ['predicted', 'actual'])

# calculating total wire length
predLength = sum([g.nodes[node]['length'] for node in g])
actualLength = sum([h.nodes[node]['length'] for node in h])
results.loc['total wire length'] = [predLength, actualLength]

# calculating number of internal resistors
predictedInternalResistorCount = sum([ramp(len(list(g.neighbors(node))) - 1) for node in gNonElectrodeNodes])
actualInternalResistorCount = len([edge for edge in h.edges if h.edges[edge]['resistanceType'] == 'int'])
results.loc['internal resistor count'] = [predictedInternalResistorCount, actualInternalResistorCount]


# calculating number of nodes
predNodes = 2 + sum([max([len(list(g.neighbors(node))), 1]) for node in gNonElectrodeNodes])
actualNodes = len(h.nodes)

# variant for add_internal_resistance_3()
#predNodes = 2 + 10 * len(gNonElectrodeNodes)
#actualNodes = len(h.nodes)
#results.loc['total nodes'] = [predNodes, actualNodes]

# calculating number isolates
predIsolates = len(list(nx.isolates(g)))
actualIsolates = len(list(nx.isolates(h)))
results.loc['num isolates'] = [predIsolates, actualIsolates]

# calculating sum of contact resistances
predTotalContactResistance = sum([g.edges[edge]['resistance'] for edge in g.edges if g.edges[edge]['resistanceType'] == 'cont'])
actualTotalContactResistance = sum([h.edges[edge]['resistance'] for edge in h.edges if h.edges[edge]['resistanceType'] == 'cont'])
results.loc['total contact res'] = [predTotalContactResistance, actualTotalContactResistance]

# calculating sum of internal resistances
# first calculate the total length bounded by contact resistors in the entire network
totalBoundedLength = 0
for node in gNonElectrodeNodes:
	# initialiing the length bounded
	boundedLength = 0
	# iterating through all possible pairs of nodes
	for edgePair in itertools.combinations(g.edges(node), 2):
		# getting coordinates
		p1 = np.array((g.edges[edgePair[0]]['x'], g.edges[edgePair[0]]['y']))
		p2 = np.array((g.edges[edgePair[1]]['x'], g.edges[edgePair[1]]['y']))

		# setting the bounded length to be the greatest distance between contact
		# edges of this node
		boundedLength = max(boundedLength, np.linalg.norm(p1 - p2))
	
	totalBoundedLength += boundedLength
# calculating sum of internal resistances
rho = g.calc_internal_resistivity(g.nwDiam, g.initialTemp)
A = np.pi * (g.nwDiam / 2)**2
predictedInternalResistance = rho * totalBoundedLength / A
actualInternalResistance = sum([h.edges[edge]['resistance'] for edge in h.edges if h.edges[edge]['resistanceType'] == 'int'])
results.loc['summed internal resistance'] = [predictedInternalResistance, actualInternalResistance]


print(results)
print('Matches?')
print(round(results['predicted'], 2) == round(results['actual'], 2))

g.to_img(outFile = 'g', showJunctions = True)
h.to_img(outFile = 'h', showJunctions = True)
