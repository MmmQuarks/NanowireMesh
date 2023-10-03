import csv
import numpy as np
import time
import NanowireMesh as nwm



percolationMultiples = np.arange(1.2, 2.5, .2)
sizes = np.arange(300,1700,200)

data = [['percolation multiple', 'side length (um)', 'time (s)']]
iter = 0
for pm in percolationMultiples:
    for s in sizes:
        iter +=1
        print('Iteration number', iter, 'of', len(percolationMultiples)*len(sizes))
        t = time.time()
        g = nwm.NanowireMesh(width = s, height = s, nwLength= 7, nwLengthSD= 3, percMultiple=pm, removeIsolates=False, removeNonPercolating=False)
        t = time.time() - t
        data.append([pm, s, t])

with open('test_scaling.csv', 'w') as file:
    writer = csv.writer(file)
    writer.writerows(data)

file.close()