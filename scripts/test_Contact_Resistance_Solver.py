import Nanowire_Mesh as nwm
import os
import networkx as nx
import numpy as np

percMultiple = 3
g = nwm.Nanowire_Mesh(width = 400, height = 400, percMultiple= percMultiple, nwLength = 7, nwLengthSD = 0)
netlistName = 'netlist'
g.toNetlist(netlistName= netlistName)
os.system('xyce ' + netlistName)
g.updateNetworkWithXyceOutput(netlistName + '.csv')

# get average contact resistance
rContact = []
lengths = []
for edge in g.edges:
    dict = g[edge[0]][edge[1]]
    if dict['resistanceType'] == 'cont':
        rContact.append(dict['resistance'])
    elif dict['resistanceType'] == 'int':
        pass
        #lengths.append(dict['length'])

diams = []
for node in g.nodes:
    diams.append(g.nodes[node]['diam'])

avgCont = np.average(rContact)
avgLen = 7 #np.average(lengths)
avgDiam = np.average(diams)

rho = nwm.calculateResistivity(avgDiam, 298.15)
Rw = rho * avgLen / (np.pi * (avgDiam/2)**2)
print('R_wire = ' + str(Rw))

Rc_estimate = nwm.solveForContactResistance(percMultiple, Rw, g.sheetResistance)
print('Average contact resistance is: ' + str(avgCont))
print('Calculated contact resistance is: ' + str(Rc_estimate))
