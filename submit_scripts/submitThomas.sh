#!/bin/bash
#
#SBATCH --job-name=make_NW_network
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-10
#SBATCH -t 48:00:00

srun ./matlab -nodisplay -nojvm -nosplash -nodesktop -r "run('/home/amwt/TPV/3D_nanowires/thomas_main.m')"
