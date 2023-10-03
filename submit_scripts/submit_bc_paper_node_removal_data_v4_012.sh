#!/bin/bash
#
#SBATCH --job-name=012
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=192000
#SBATCH -t 48:00:00
#SBATCH --output="/home/amwt/TPV/2DNanowires/Data/evolution_under_node_removal/bc_paper_node_removal_data_v4_012/slurm.out"
#SBATCH --error="/home/amwt/TPV/2DNanowires/Data/evolution_under_node_removal/bc_paper_node_removal_data_v4_012/slurm.err"

source /home/amwt/.bashrc
module load intel
conda activate envnw
mpirun python /home/amwt/TPV/2DNanowires/bc_paper_node_removal_data_v4.py True True /home/amwt/TPV/2DNanowires/Data/evolution_under_node_removal/bc_paper_node_removal_data_v4_012/
