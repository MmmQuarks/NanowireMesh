import NanowireMesh as nwm
import networkx as nx

g = nwm.NanowireMesh()
thing = {node : {'random ting' : node, 'other ting' : int(node)*100 } for node in g.nodes}
nx.set_node_attributes(g, thing)
print(g.nodes[0])
g.to_pickle('testpickle')
h = nwm.NanowireMesh(inPickle = 'testpickle.p')
print(h.nodes[0])
