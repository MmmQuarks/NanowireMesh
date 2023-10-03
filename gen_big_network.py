import NanowireMesh as nwm

percolationMultiples = [1.20863351, 1.52287823, 2.07884964, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]

for pm in percolationMultiples:
	print('current perc mult is', str(pm))

	# generate network
	g = nwm.NanowireMesh(width = 1000, height = 1000, nwLength = 7, nwLengthSD = 3, percMultiple = pm, removeNonPercolating = False, removeIsolates = False, addInternalRes = True)

	# make folder name
	folderName = '/home/amwt/TPV/2DNanowires/Data/1cm_perc_mult_' + str(pm)

	# write to csv
	print('writing to CSV')
	g.to_csv(folderName = folderName)
	print('writing to netlist')
	g.to_netlist(netlistName = folderName + '/netlist')

