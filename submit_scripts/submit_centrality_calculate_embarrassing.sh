#!/bin/bash
#
#SBATCH --job-name=centrality_calculate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=19200
#SBATCH -t 48:00:00
#SBATCH --output="/home/amwt/TPV/2DNanowires/parallel_centrality_embarrassing_calc.out"
#SBATCH --error="/home/amwt/TPV/2DNanowires/parallel_centrality_embarrassing_calc.err"
#SBATCH --array=1-10

source ~/.bashrc
module load intel
conda activate 2DNanowires_env
srun python -u parallel_centrality_calc_embarassing.py $SLURM_ARRAY_TASK_ID
