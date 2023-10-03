import NanowireMesh as nwm
from os import system
from matplotlib import pyplot as plt
import numpy as np
import math

g = nwm.NanowireMesh(height = 250, width = 250, nwLength = 7, nwLengthSD = 3, percMultiple = 1.5, removeNonPercolating=False, addInternalRes= True)
g.to_netlist()
system('xyce netlist')
g.update_with_xyce_output('netlist.csv')

powerMassRatio = [g.edges[edge]['power'] / g.edges[edge]['mass'] for edge in g.edges]

bins = np.linspace(math.ceil(min(powerMassRatio)), math.floor(max(powerMassRatio)), 500)
plt.hist(powerMassRatio, bins = bins)
plt.show()