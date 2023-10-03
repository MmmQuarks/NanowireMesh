#!/bin/bash
#
#SBATCH --job-name=node_removal_data
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=192000
#SBATCH -t 48:00:00
#SBATCH --output="/home/amwt/TPV/2DNanowires/Data/evolution_under_node_removal/bc_paper_node_removal_data_v4_009/slurm.out"
#SBATCH --error="/home/amwt/TPV/2DNanowires/Data/evolution_under_node_removal/bc_paper_node_removal_data_v4_009/slurm.err"

source /home/amwt/.bashrc
module load intel
conda activate envnw
mpirun python /home/amwt/TPV/2DNanowires/bc_paper_node_removal_data_v4.py
