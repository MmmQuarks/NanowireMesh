import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np
def main():
	def pdf(x,y):
		heavisideArg = (10**2 - (x - 50)**2 - (y - 50)**2)
		return 1./5 *(1. + 4 * np.heaviside(heavisideArg, 1))
	x, y  = draw(pdf, xBounds=[0, 100], yBounds=[0,100], numDraws = 4000)
	fig, ax = plt.subplots()
	fig = plt.scatter(x,y, s= .5)
	ax.set_aspect('equal')
	plt.ylim((0,100))
	plt.xlim((0,100))
	plt.show()

def draw(pdf, xBounds = [0,1], yBounds = [0,1], numDraws = 100):
	# note that pdf should be normalized such that it only varies between 0 and 1
#	procedure:
	#1. generate x,y,z vals
	#2. if z < p(x,y) add (x,y) to list of coordinates
	
	xSample = []
	ySample = []
	while len(xSample) < numDraws:
		testX = (xBounds[1] - xBounds[0]) * np.random.random() + xBounds[0]
		testY = (yBounds[1] - yBounds[0]) * np.random.random() + yBounds[0]
		testZ = np.random.random()
		if testZ < pdf(testX, testY):
			xSample.append(testX)
			ySample.append(testY)
	return xSample, ySample
if __name__ == "__main__":
	main()	
	
