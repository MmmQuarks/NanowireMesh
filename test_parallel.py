#!/home/amwt/mambaforge/envs/envnw/bin/python

#SBATCH --job-name=tp
#SBATCH --output=testpar.out
#SBATCH --error=testpar.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=12000
#SBATCH --exclude=n16
#SBATCH --time=48:00:00


import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import os
import time

def maker(n):
	print('Starting job {}'.format(n), flush = True)
	numrows = int(np.random.rand() * 1E7)
	data = [dict(cpu = n, row = i) for i in range(numrows)]
	df = pd.DataFrame(data)
	df.to_csv('test_par_job_{}.csv'.format(n))
	print('Done with job {}'.format(n), flush = True)



if __name__ == '__main__':
	par_start = time.time()
	num_cores = int(os.environ['SLURM_CPUS_PER_TASK'])
	with Pool(num_cores) as pool:
		pool.starmap(maker, [[int(n)] for n in range(num_cores)])
	print('Par Elapsed = {}'.format(time.time() - par_start))

	seq_start = time.time()
	for n in range(num_cores):
		maker(n)
	print('Seq Elapsed = {}'.format(time.time() - seq_start))
