from mcts import mcts
import networkx as nx
import pdb
from MCTS_State import State
from pathlib import Path
import sys
import time
import argparse
import pickle
import NanowireMesh as nwm
import os

parser = argparse.ArgumentParser(description = "Execute Monte Carlo Tree Search to optimize Nanowire Networks")
parser.add_argument('--network-size', type = float, 
	help = 'the side length of the networks to generate. we assume height == width here')
parser.add_argument('--array-id', type = int, default = 0,
		    help = 'index of this job if using a job array')
parser.add_argument('--results-folder', type = str)
parser.add_argument('--wall-time', type = float, default = 0,
	help = 'The maximum time, in hours, we want this code to run for. Default is 0 which means no wall time. When code gets close to wall, will preemptively save safely and abort.'
	)
args = parser.parse_args(sys.argv[1:])


def main():
	global args
	run_start = time.time()
	best_states_filepath = Path(
		args.results_folder,
		'best_states{}.p'.format(str(args.array_id).zfill(2))
	)
	# if we are continuing a previous run
	if os.path.exists(best_states_filepath):
		with open(best_states_filepath, 'rb') as file:
			best_states = pickle.load(file)
	else: # otherwise make a new network and initialize a new state
		g = nwm.NanowireMesh(
			percMultiple = 2,
			height = args.network_size,
			width = args.network_size,
			rcMean = 'bellew',
			makeEmpty = False
		)
		print('Calculating Betweenness Centrality')
		bc = nx.betweenness_centrality(
			g,
			weight = 'resistance'
		)
		nx.set_node_attributes(
			g,
			bc,
			name = 'betweenness_centrality'
		)
		initial_state = State(g, parent = None)
		best_states = [initial_state]
	tree = mcts(iterationLimit = 2000)
	max_step_time = 0
	while not best_states[-1].isTerminal():
		step_start = time.time()
		best_states[-1].best_action = best_action = tree.search(
			initialState = best_states[-1]
		)
		step_time = time.time() - step_start
		max_step_time = max(max_step_time, step_time)

		# get next state
		next_state = tree.root.children[best_action].state
		best_states.append(next_state)
		# save progress
		with open(best_states_filepath, 'wb') as file:
			pickle.dump(best_states, file)
		run_time = time.time() - run_start
		wall_time = args.wall_time * 3600 # converting wall time from hours to seconds
		# if we have less than max_step_time between now and the wall time, terminate safely
		# only do this if wall_time > 0
		if wall_time - run_time < max_step_time and wall_time > 0:
			break
	
	if best_states[-1].isTerminal():
		print('\n\n\nCompleted Tree Search')

if __name__ == '__main__':
	main()
