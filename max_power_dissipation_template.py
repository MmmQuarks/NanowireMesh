#!/home/amwt/mambaforge/envs/envnw/bin/python

#SBATCH --job-name=template
#SBATCH --output=template.out
#SBATCH --error=template.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=12000
#SBATCH --exclude=n16
#SBATCH --time=48:00:00


import sys
sys.path.insert(0, '/home/amwt/TPV/2DNanowires')
import pandas as pd
import numpy as np
from scipy.optimize import minimize_scalar
from multiprocessing import Pool, cpu_count
from pathlib import Path
import os
import NanowireMesh as nwm
import PardisoSolver
from copy import deepcopy
import SolverUtils
import networkx as nx

def power(g, edge, R):
	if R == 0:
		raise ValueError('R must be >= 0. Something wrong with edge {}'.format(edge))
	g.edges[edge]['resistance'] = R
	PardisoSolver.solve(g, voltage = 1, inplace = True)
	# we make the objective function negative to work with minimize_scalar
	return -g.edges[edge]['power']

def minimize_scalar_wrapper(g, edge):
	original_resistance = g.edges[edge]['resistance']
	ans = minimize_scalar(lambda R : power(g, edge, R), bracket = (1,1E7))
	max_power = g.edges[edge]['power']
	g.edges[edge]['resistance'] = original_resistance
	return (edge, ans, max_power)



if __name__ == '__main__':
	num_cores = int(os.environ['SLURM_CPUS_PER_TASK'])

	network_path = Path( '<NETWORK_PATH>' )
	g = nwm.NanowireMesh(
		makeEmpty = False,
		inPickle = network_path
		)
	SolverUtils.prune_network(g, pop = False)

	iterable = [(g, edge) for edge in g.edges]
	
	# parallel processing
	with Pool(num_cores) as pool:
		starmap_results = pool.starmap(
			minimize_scalar_wrapper,
			iterable
			)

		# converting results to csv
		df_list = []
		for (edge, ans, max_power) in starmap_results:
			df_list.append(
				dict(
					n1 = edge[0],
					n2 = edge[1],
					R_opt = ans['x'],
					P_max = max_power,
					success = ans['success'],
					message = ans['message']
				)
			)
		df_path = Path(
			network_path.parent,
			network_path.stem + 'maximal_power.csv'
			)
		pd.DataFrame(df_list).to_csv(df_path)
	print('Executed pool successfully. Closing Program.')
			

	

