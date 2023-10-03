import Nanowire_Mesh as nwm
import networkx as nx
import numpy as np
import os

#inFolder = 'thomas_B2'
#inFolders = ['thomas_B1', 'thomas_B2', 'thomas_B3']
inFolders = ['thomas_B1', 'thomas_B2', 'thomas_B3']

expResults = {'B2' : {'resistance' : 329, 'percolationMultiple' : 1.52, 'tortuosity' : 1.063}}

# establishing trial parameters
for inFolder in inFolders:
        g = nwm.Nanowire_Mesh(inFolder = inFolder)
        simResults = {}

        rcMeanStart = 5
        rcSDStart = 0
        trialParams = {}

        for n in range(5):
                name = ''.join(['t',str(n)])
                rcMean = rcMeanStart * np.sqrt(10)**n
                rcSD = rcSDStart * np.sqrt(10)**n
                trialParams.update({name : {'rcMean' : rcMean, 'rcSD' : rcSD}})


        resistanceTypes = nx.get_edge_attributes(g, 'resistanceType')

        # find edges corresponding to contact resistances
        contacts = []
        for edge in resistanceTypes.keys():
                if resistanceTypes[edge] == 'cont':
                        contacts.append(edge)

        for trialKey in trialParams.keys():
                rcMean = trialParams[trialKey]['rcMean']
                rcSD = trialParams[trialKey]['rcSD']
                rcList = np.random.normal(rcMean, rcSD, len(contacts))
                rcDict = {}
                for counter, edge in enumerate(contacts):
                        rcDict.update({edge : {'resistance' : rcList[counter] } })
                nx.set_edge_attributes(g, rcDict)
                netlistName = '/'.join([inFolder,str(trialKey), 'netlist'])
                g.toNetlist(netlistName = netlistName)
                command = ' '.join(['xyce', netlistName])
                os.system(command)
                g.updateNetworkWithXyceOutput('.'.join([netlistName,'csv']))
                simResults.update({trialKey: {'resistance' : g.sheetResistance}})
                print(str(g.sheetResistance))
        print(trialParams)
        print(simResults)