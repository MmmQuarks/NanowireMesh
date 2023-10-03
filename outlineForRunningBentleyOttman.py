f = [[0,(-0.001,0.0),(700.001,0)]] 	#I included this format because I didn't have a network that offset the top and bottom electrodes so feel free to ignore it if you have that

for node in g.nodes(): 	    
	f.append( [node] +  [g.nodes('endpoints')[node][0], g.nodes('endpoints')[node][1]]) #read out the endpoint data in the way the bentleyOttman program likes
f.sort(key = lambda x: x[0]) #sort by height but we might not need to. We get the line endpoints themselves at the end anyway.
del f[1] #this is because i'm manually adding an offset
[row.pop(0) for row in f] #getting rid of the wire numbers. This wouldn't be here and the code would be even shorter if I don't include the ordering, but I'm doing it just in case.

import #whatever you call the function as BO

BO.isect_segments_include_segments(f) #outputs a list[ tuple( tuple(this is the intersction), list[ tuple(2x tuples endpts), tuple(2x tuples endpts) ] ) ]