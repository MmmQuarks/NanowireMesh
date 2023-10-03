import NanowireMesh as nwm
from os import system
import time
import csv

g = nwm.NanowireMesh(width = 1000, height = 1000, nwLength = 7, nwLengthSD = 4, percMultiple=1.75)

n = 0
folder = '/home/amwt/TPV/2DNanowires/Frames3'
g.to_pickle(folder + '/' + str(n) + '.p')
print('Nodes in network:', len(g.nodes))
oldNumNodes = len(g.nodes)
voltage = 0.05 * 2
while g.isPercolating:
	n+= 1
	g.to_netlist(voltage = voltage)
	g.voltage = voltage
	system('/home/amwt/TPV/2DNanowires/installations/build_Xyce-6.10_parallel/src netlist')
	g.update_with_xyce_output('netlist.csv')
	oldNumNodes = len(g.nodes)
	g.propagate_time(.01, option = 2)
	if oldNumNodes == len(g.nodes):
		voltage = 1.1 * voltage
	g.to_pickle(folder + '/' + str(n) + '.p')
	print('Nodes in network:', len(g.nodes))
	print('Is percolating:', g.isPercolating)
	print('voltage = ', voltage)
	# time.sleep(2)

# data = []	
# for nn in range(n+1):
# 	g = nwm.NanowireMesh(inPickle = 'Frames3/' + str(nn) + '.p')
# 	x1 = [g.nodes[node]['endpoints'][0][0] for node in g.nodes]
# 	x2 = [g.nodes[node]['endpoints'][1][0] for node in g.nodes]
# 	y1 = [g.nodes[node]['endpoints'][0][1] for node in g.nodes]
# 	y2 = [g.nodes[node]['endpoints'][1][1] for node in g.nodes]
# 	temps = [g.nodes[node]['temp'] for node in g.nodes]
# 	gData = [['x1', 'x2', 'y1', 'y2', 'temp']]
# 	gData = gData + [[x1[i], x2[i], y1[i], y2[i], temps[i]] for i in range(len(y2))]
# 	fileName = str(nn) + '_plotdata.csv'
# 	with open(fileName, mode = 'w') as file:
# 		writer = csv.writer(file, delimiter = ',')
# 		writer.writerows(gData)
