#!/home/amwt/mambaforge/envs/envnw/bin/python

import pdb
from multiprocessing import Pool
import pandas as pd
from pathlib import Path
import os
import pickle
import sys
if str(Path.cwd()) not in sys.path:
	sys.path.insert(0, str(Path.cwd()))
import NanowireImageSim as nis

if __name__ == '__main__':
	jobdir = Path(os.environ['JOBDIR'])
	num_images = int(os.environ['NUM_IMAGES'])
	chunksize = int(os.environ['CHUNKSIZE'])
	args = pd.read_csv(
		Path(jobdir, 'args_NanowireImageSim.csv')
		)
	
	# drawing parameters
	draws = nis.param_draws(
		args,
		num_images
	)
	# saving draws
	draws_path = Path(jobdir, 'draws.p')
	with open(draws_path, 'wb') as file:
		pickle.dump(draws.to_dict(), file)
	
	draws_list = [row for idx, row in draws.iterrows()]
	sample_folder_len = len(str(num_images))
	numbers = [str(n).zfill(sample_folder_len) for n in range(num_images)]
	paths = [Path(jobdir, num) for num in numbers]
	iterable = zip(
		paths,
		draws_list
	)
	pool = Pool(
		int(os.environ['SLURM_CPUS_PER_TASK'])
	)
	pool.starmap(
		nis.make_sample,
		iterable,
		chunksize = chunksize
	)
	pool.close()
