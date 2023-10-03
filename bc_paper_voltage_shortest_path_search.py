import NanowireMesh as nwm
import numpy as np
import networkx as nx

# load the example network
g = nwm.NanowireMesh(inPickle = '200x200_percMult1.5_nwLen10_Rc10.p')

voltage = nx.get_node_attributes(g, 'voltage')
#print(voltage)

def getSuitableNeighbors(g, thisNode, neighborsToExclude):
	neighbors = list(g.neighbors(thisNode))
	neighborsThatArentDanglingEnds = {node for node in neighbors if len(list(g.neighbors(node))) > 1}
	neighborsBelowVoltage = {node for node in neighbors if g.nodes[node]['voltage'] < g.nodes[thisNode]['voltage']}
	neighborsToExclude = set(neighborsToExclude[thisNode])
	# only want to return the intersection of these two sets (i.e. neighbors that are both down potential and not dangling)
	return (neighborsBelowVoltage & neighborsThatArentDanglingEnds) - neighborsToExclude

def getNextNode(g, thisNode, suitableNeighbors):
	# iterate through suitable neighbors and find the one with the with the highest  Voltage / vertical distance value
	bestNode = thisNode
	bestVoltageDrop = abs(g.nodes[g.topElectrode]['voltage'] - g.nodes[g.bottomElectrode]['voltage'])

	for node in suitableNeighbors:
		voltageDrop = abs(g.nodes[node]['voltage'] - g.nodes[thisNode]['voltage'] )
		if voltageDrop < bestVoltageDrop:
			bestVoltageDrop = voltageDrop
			bestNode = node
	return bestNode

pathComplete = False
path = [g.topElectrode]
# the below structure contains the nodes that should not be traveled to next from each node
neighborsToExclude = {node : [] for node in g.nodes}
while not pathComplete:
	print('Finding neighbors for node', path[-1])
	suitableNeighbors = getSuitableNeighbors(g, path[-1], neighborsToExclude)
	if len(suitableNeighbors) == 0:
		print('No suitable neighbors found. Moving back one node and removing present node from consideration as a next node') 
		previousNode = path[-2]
		thisNode = path[-1]
		del path[-1]
		neighborsToExclude[previousNode].append(thisNode)
	else:
		nextNode = getNextNode(g, path[-1], suitableNeighbors)
		path.append(nextNode)
		if nextNode == g.bottomElectrode:
			pathComplete = True
			print('Path is complete')

# visualize the thing
to_show_dict = {node : 1 if node in path else 0 for node in g.nodes}
g.to_img(z = to_show_dict, title = 'Shortest Paths by Minimizing Size of Voltage Drops', outFile = 'bc_paper_voltage_shortest_path_search.png')
