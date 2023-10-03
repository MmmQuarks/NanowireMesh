import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

df = pd.read_csv('~/Documents/Research/TPV/2D_nanowires/3D_nanowires/Data/spice_files/2019_02_16_compiled_data',sep = '\|\|',engine = 'python')

#make nondimensionalized data
sigma = 5.63726 #constant in percolation calculations
l = df.nwLength
L = df.length #sample side length
alpha = df.percolationMultiple

data = pd.DataFrame( {'nl2' : alpha * sigma, 'Rs/Rj' : df.sheetResistance / df.junctionResistanceMean, 'sj/Rj' : df.junctionResistanceSD / df.junctionResistanceMean})
data['sj/Rj'] = round(data['sj/Rj'],2)
thispdata = data[ data['sj/Rj']==0]
label = r'$R_s/R_J$ for $\sigma_J/R_J = 0$'
thispdata = thispdata.groupby(['nl2']).agg({'Rs/Rj' : ['mean', 'std']})
#fitting the data, assuming it has the form
# (Rs/Rj) = C (nl^2)^gamma
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(thispdata.index), np.log(thispdata.values[:,0]))
fit_x = np.linspace(thispdata.index.min(), thispdata.index.max(), 100)
fit_y = np.exp(intercept) * fit_x**slope
C = np.exp(intercept)
gamma = slope

plist = [[thispdata, label, fit_x, fit_y, C, gamma, r_value]]

plt.figure(figsize = (8*1.2,6 * 1.2))
plt.errorbar(thispdata.index, thispdata.values[:,0], yerr = thispdata.values[:,1],fmt = 'o',label = label)
plt.plot(fit_x, fit_y,label = label + " fit")
plt.xlabel('number density x length^2',fontsize = 'large')
plt.ylabel(r'$\frac{R_{sheet} }{R_{junction} }$', rotation = 'horizontal', fontsize = 'xx-large')
plt.title('Non-Dimensionalized Plot of Geometrical vs. Electrical Properties')

for i in range(0,5):
	thispdata = data[ data['sj/Rj'] > i*0.1]
	thispdata = thispdata[thispdata['sj/Rj'] <= (i+1)*0.1 ]
	label = r"$R_s/R_J$ for $\sigma_J/Rj$ in (" + str(round(i*0.1,2)) + ', ' + str( round((i + 1)*0.1,2)) + "]"

	thispdata = thispdata.groupby(['nl2']).agg({'Rs/Rj' : ['mean', 'std']})
	slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(thispdata.index), np.log(thispdata.values[:,0]))
	fit_x = np.linspace(thispdata.index.min(), thispdata.index.max(), 100)
	fit_y = np.exp(intercept) * fit_x**slope
	C = np.exp(intercept)
	gamma = slope

	plist.append([thispdata, label, fit_x, fit_y, C, gamma, r_value])

	plt.errorbar(thispdata.index, thispdata.values[:,0], yerr = thispdata.values[:,1],fmt = 'o', label = label)
	plt.plot(fit_x, fit_y, label = label + " fit")


#make legend for plot (note we must go in backwards order)
# leg = []
# for i in reversed(range(0,len(plist))):
# 	leg.append( plist[i][1] + ' fit')
# 	leg.append( plist[i][1])

# plt.legend(leg)
f = open('fit_summary.txt','w+')
#f.write('SD range,C,gamma,r_value')
print('SD range||C||gamma||r_value', file = f)
for i in range(len(plist)):
	print('%s||%f||%f||%f' % (plist[i][1], round(plist[i][4],2) , round(plist[i][5],2), round(plist[i][6],5) ) , file = f)
f.close()

plt.legend()

plt.show()


# plt.legend(['Fit: y =' + str(round(np.exp(intercept),2)) + "x^" + str(round(slope,3)), 'Simulations'])


# print("C = " + str(np.exp(intercept)))
# print(" a = " + str(slope))


# plt.show()