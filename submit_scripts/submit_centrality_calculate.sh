#!/bin/bash
#
#SBATCH --job-name=net_gen
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=32000
#SBATCH -t 48:00:00
#SBATCH --output="/home/amwt/TPV/2DNanowires/parallel_centrality_calc.out"
#SBATCH --error="/home/amwt/TPV/2DNanowires/parallel_centrality_calc.err"

source ~/.bashrc
module load intel
which mpirun
conda activate 2DNanowires_env
mpirun -n 6 python /home/amwt/TPV/2DNanowires/parallel_centrality_calc.py
