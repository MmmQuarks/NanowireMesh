import NanowireMesh as nwm
import numpy as np
import argparse
import sys
from pathlib import Path
from copy import deepcopy
from tqdm import tqdm


def update_contact_resistance(g, Rc):
	contact_edges = filter(lambda e : g.edges[e]['resistanceType'] == 'cont', g.edges)
	for e in contact_edges:
		g.edges[e]['resistance'] = Rc

def impute(g, Rt, solver, disable_tqdm = False):
	'''
	Calculates an average effective contact resistance given a graph and a target
	network resistance

	args
		g = the input graph
		Rt = target resistance in Ohms
		solver = either 'numpy' or 'pardiso'
		disable_tqdm = if True, hides progress bars for tqdm. Default is False (shows progress bars)

	returns
		Rc_imputed = the imputed average effective contact resistance for the graph g that yields
					total network resistance of Rt
	'''
	g = deepcopy(g)
	if solver == 'numpy':
		import NumpySolver as SolverModule
	elif solver == 'pardiso':
		import PardisoSolver as SolverModule

	Rc_test = [Rt * 10**(tenpow) for tenpow in range(-3,5)]
	Rt_test = []
	
	# first pass
	pbar1 =  tqdm(Rc_test, desc = 'Pass 1', disable = disable_tqdm)
	for Rc in pbar1:
		update_contact_resistance(g, Rc)
		SolverModule.solve(
			g,
			voltage = 1,
			inplace = True
		)
		Rt_test.append(g.sheetResistance)
		if Rc == Rc_test[-1]:
			b = np.array(Rt_test)
			A = np.array(
				[(1, Rc) for Rc in Rc_test]
			)	
		
			(alpha, beta), residuals, rank, s = np.linalg.lstsq(A, b, rcond = None)
			Rc_imputed = (Rt - alpha) / beta
		
			# calculating error
			update_contact_resistance(g, Rc_imputed)
			SolverModule.solve(
				g,
				voltage = 1,
				inplace = True
			)
			err = abs(Rt - g.sheetResistance)
			rel_err = err / min(Rt, g.sheetResistance)
			trial_counter = 1
			pbar1.set_description('Pass 1, rel_err in Rs = {:.5f}'.format(rel_err))
	while rel_err > 0.001:
		trial_counter += 1
		Rc_test = Rc_imputed * np.linspace(1-rel_err, 1+rel_err ,10)
		Rt_test = []
		pbar2 = tqdm(
			Rc_test,
			desc = 'Pass {}'.format(
				trial_counter
			),
			disable = disable_tqdm
		)
		for Rc in pbar2:
			update_contact_resistance(g,Rc)
			SolverModule.solve(g, voltage = 1, inplace = True)
			Rt_test.append(g.sheetResistance)
			if Rc == Rc_test[-1]:
				b = np.array(Rt_test)
				A = np.array(
					[(1, Rc) for Rc in Rc_test]
				)
				(alpha, beta), residuals, rank, s = np.linalg.lstsq(A, b, rcond = None)
				Rc_imputed = (Rt - alpha) / beta
				# solving with imputed resistance and testing error
				update_contact_resistance(g, Rc_imputed)
				SolverModule.solve(
					g,
					voltage = 1,
					inplace = True
				)
				err = abs(Rt - g.sheetResistance)
				rel_err = err / min(Rt, g.sheetResistance)
				pbar2.set_description(
					'Pass {}, rel_err in Rs = {:.5f}'.format(
						trial_counter,
						rel_err
					)
				)


	if Rc_imputed <= 0:
		raise ValueError('Value of Rc_imputed is negative. This is unphysical.')
	else:
		return Rc_imputed

parser = argparse.ArgumentParser(description = 'Given some network G and a target resistance R, find the contact resistance Rc that gives G total resistance R')

parser.add_argument('G', type = str, help = 'path to the graph pickle')
parser.add_argument('Rt', type = float, help = 'The target resistance for the entire network')
parser.add_argument('solver', choices = ['numpy', 'pardiso'], help = 'Which solver to use ')


if __name__ == '__main__':
	args = parser.parse_args(sys.argv[1:])
	g = nwm.NanowireMesh(
		makeEmpty = False,
		inPickle = args.G
		)
	
	Rc_imputed = impute(
		g,
		args.R,
		args.solver
	)

	print(
		'\nEffective contact resistance = {:.4f}'.format(Rc_imputed)
	)