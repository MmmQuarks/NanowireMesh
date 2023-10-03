import Nanowire_Mesh as nwm
import os
import math

# here we create a network, g, and find its properties with Xyce.
# we then write g to file and initialize a new network h from that file.
# we then compare all the properties of h and g to make sure that the reading and writing 

def findRowMatches(data1, data2):
    rowMatches = []
    rel_tol = 1e-5
    if type(data1) == type({}):
        for key in data1.keys():
            if type(data1[key]) == type([]) or type(data1[key]) == type(()):
                pass
                # list1 = data1[key]
                # list2 = data2[key].strip()
                # listTruths = []
                # for r in range(len(list1)):
                #     for c in range(len(list1[r])):
                #         print(list1[r][c])
                #         print(list2[r][c])
                #         print(nwm.isNumber(list2[r][c]))
                #         val1 = float(list1[r][c])
                #         val2 = float(list2[r][c])
                #         listTruths.append( math.isclose(val1, val2, rel_tol = rel_tol))
            elif nwm.isNumber(data1[key]):
                truthValue = math.isclose(data1[key], data2[key], rel_tol = rel_tol)
            else:
                truthValue = (data1[key] == data2[key])
            rowMatches.append(truthValue)
    return rowMatches

debugMode = False

g = nwm.Nanowire_Mesh()
folderName = 'testing_read_write'
g.toCSV(folderName = folderName)

netlistName = '/'.join([folderName, 'netlist'])
g.toNetlist(netlistName = netlistName)

os.system('xyce ' + netlistName)
g.updateNetworkWithXyceOutput(inFile = netlistName + '.csv')
g.toCSV(folderName = folderName)

h = nwm.Nanowire_Mesh(inFolder= folderName)

print('Checking that all g are in h and match')
matches = []
for edge in g.edges:
    gData = g[edge[0]][edge[1]]
    hData = h[edge[0]][edge[1]]
    matches.append(findRowMatches(gData,hData))
print('edges match = ' + str(all(matches)))

matches = []
for node in g.nodes:
    gData = g.nodes[node]
    hData = h.nodes[node]
    matches.append(findRowMatches(gData,hData))
print('nodes match = ' + str(all(matches)))

print('Checking that all h are in g and match')
matches = []
for edge in h.edges:
    gData = g[edge[0]][edge[1]]
    hData = h[edge[0]][edge[1]]
    matches.append(findRowMatches(hData,gData))
print('edges match = ' + str(all(matches)))

matches = []
for node in h.nodes:
    gData = g.nodes[node]
    hData = h.nodes[node]
    matches.append(findRowMatches(hData, gData))
print('nodes match = ' + str(all(matches)))

 print('note that list data is not being imported correctly, so no list data was verified')