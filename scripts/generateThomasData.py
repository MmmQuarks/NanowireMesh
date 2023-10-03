import NanowireMesh as nwm
import numpy as np

percMultiples = np.array([ 1.20863351,  1.52287823,  2.07884964])
width = height = 700

for n, percMultiple in enumerate(percMultiples,1):
    g = nwm.NanowireMesh(width = width, height = height, nwLength = 7, nwLengthSD = 3, percMultiple = percMultiple)
    failureCount = 0
    folderName = 'thomas_B' + str(n)
    g.toCSV(folderName)