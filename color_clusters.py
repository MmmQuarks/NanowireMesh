import networkx as nx
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt


def color_clusters(img_mat):
	def add_edge_if_exists(img_mat,  g, p1, p2):
		if int(img_mat[p1[0], p1[1]]) == 0 and int(img_mat[p2[0], p2[1]]) == 0:
			g.add_edge(p1, p2)

	g = nx.Graph()

	print('Converting matrix to graph representation')
	for row in tqdm(range(len(img_mat))):
		for col in range(len(img_mat[row, :])):
			if img_mat[row, col] == 0:
				# add node to graph
				# g.add_node((row, col))
				p1 = (row, col)

				# testing the row above this for connection
				if row != 0:
					p2 = (row - 1, col)
					add_edge_if_exists(img_mat, g, p1, p2)
				
				# testing row below for connection
				if row != len(img_mat) - 1:
					p2 = (row + 1, col)
					add_edge_if_exists(img_mat, g, p1, p2)

				# testing column to the left for connection
				if col != 0:
					p2 = (row, col - 1)
					add_edge_if_exists(img_mat, g, p1, p2)

				# testing column to the right for connection
				if col != len(img_mat[row, :]) -1:
					p2 = (row, col + 1)
					add_edge_if_exists(img_mat, g, p1, p2)
	
	# getting all clusters
	print('Finding connected clusters')
	clusters = {}
	clusterNum = 0
	nodesToRemove = len(g.nodes)
	pbar = tqdm(total = nodesToRemove+1)
	
	while len(g.nodes) > 0:
		node = list(g.nodes)[0]
		T = nx.dfs_tree(g, node)
		connectedCluster = list(T.nodes)
		clusters.update({clusterNum : connectedCluster})
		pbar.update(len(T.nodes))
		g.remove_nodes_from(T.nodes)
		clusterNum += 1

	pixelsAfterClustering = sum([len(clusters[key]) for key in clusters.keys()])

	# # converting clusters back into numpy array
	img_output = np.zeros((len(img_mat), len(img_mat[0,:])))
	print('Converting graph back into image matrix')
	allClusterNums = [key for key in clusters.keys()]
	for clusterNum, nodes in tqdm(clusters.items()):
		if len(nodes) >= 10:
			randomIndex = np.random.randint(0, len(allClusterNums))
			for inds in nodes:
				img_output[inds[0], inds[1]] = allClusterNums[randomIndex]
			del allClusterNums[randomIndex]

	# masking parts of the image so they will be black

	img_output = np.ma.masked_where(img_output == 0, img_output)
	cmap = plt.cm.tab10
	cmap.set_bad(color = 'black')
	plt.imshow(img_output, cmap = cmap)
	plt.show()



	# graphs = nx.connected_component_subgraphs(g)
	# clusters = {}
	# for n, graph in enumerate(tqdm(graphs)):
	# 	clusters.update({n : list(graph.nodes)})

	return img_output