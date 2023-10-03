import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

df = pd.read_csv('~/Documents/Research/TPV/2D_nanowires/3D_nanowires/Data/spice_files/2019_02_16_compiled_data',sep = '\|\|',engine = 'python')
df = df[df.junctionResistanceSD == 0]

#make nondimensionalized data
sigma = 5.63726 #constant in percolation calculations
l = df.nwLength
L = df.length #sample side length
alpha = df.percolationMultiple
#data = pd.DataFrame( {'nl2' : alpha * sigma + alpha * l/L + alpha * 5.5* l**2/L**2 ,
#						'Rs/Rj' : df.sheetResistance / df.junctionResistanceMean})
data = pd.DataFrame( {'nl2' : alpha * sigma, 'Rs/Rj' : df.sheetResistance / df.junctionResistanceMean})

plotData =  data.groupby(['nl2']).agg({'Rs/Rj' : ['mean', 'std']})

#fitting the data, assuming it has the form
# (Rs/Rj) = C (nl^2)^gamma
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(plotData.index), np.log(plotData.values[:,0]))

fit_x = np.linspace(plotData.index.min(), plotData.index.max(), 100)
fit_y = np.exp(intercept) * fit_x**slope

plt.figure(figsize = (8*1.2,6 * 1.2))
plt.errorbar(plotData.index, plotData.values[:,0], yerr = plotData.values[:,1],fmt = 'o')
plt.plot(fit_x, fit_y)
plt.xlabel('number density x length^2',fontsize = 'large')
plt.ylabel(r'$\frac{R_{sheet} }{R_{junction} }$', rotation = 'horizontal', fontsize = 'xx-large')
plt.title('Non-Dimensionalized Plot of Geometrical vs. Electrical Properties')
plt.legend(['Fit: y =' + str(round(np.exp(intercept),2)) + "x^" + str(round(slope,3)), 'Simulations'])


print("C = " + str(np.exp(intercept)))
print(" a = " + str(slope))


plt.show()