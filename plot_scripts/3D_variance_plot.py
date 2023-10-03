import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import stats

df = pd.read_csv('~/Documents/Research/TPV/2D_nanowires/3D_nanowires/Data/spice_files/2019_02_16_compiled_data',sep = '\|\|',engine = 'python')

#make nondimensionalized data
const = 5.63726 #constant in percolation calculations
l = df.nwLength
L = df.length #sample side length
alpha = df.percolationMultiple
#data = pd.DataFrame( {'nl2' : alpha * const + alpha * l/L + alpha * 5.5* l**2/L**2 ,
#						'Rs/Rj' : df.sheetResistance / df.junctionResistanceMean})
data = pd.DataFrame( {'nl2' : alpha * const, 'Rs/Rj' : df.sheetResistance / df.junctionResistanceMean, 'sj/Rj' : df.junctionResistanceSD / df.junctionResistanceMean})
data['Rs/Rj'] = round(data['Rs/Rj'],2)
m =  data.groupby(['nl2', 'sj/Rj']).mean()

#getting the values of the indices from the stupid multi index object that I don't understand


u = m.unstack()
rsrj = u.values # z axis
nl2 = u.index # x axis
sjrj = u['Rs/Rj'].columns # y axis

sjrj, nl2 = np.meshgrid(sjrj, nl2) #this order is important because it matches the dimensions of the array of R_sheet / R_junction






fig = plt.figure(figsize = (8*1.2,6 * 1.2))
ax = Axes3D(fig)

ax.set_xlabel(r'$D = n l^2$', fontsize = 'xx-large')
ax.set_ylabel(r'$\sigma_J/R_J$', fontsize = 'xx-large')
ax.set_zlabel(r'$R_s / R_J$', fontsize = 'xx-large',rotation = 'horizontal')
plt.title(r'$\frac{R_s}{R_J}$ vs $nl^2$ vs $\frac{\sigma_J}{R_J}$', fontsize = 'xx-large')




ax.view_init(10,-35)
ax.plot_surface(nl2, sjrj, u.values, cmap = cm.coolwarm, linewidth=0, antialiased=True)

plt.show()

