#!/bin/bash
#
#SBATCH --job-name=1.8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=192000
#SBATCH -t 48:00:00
#SBATCH --output="/home/amwt/TPV/2DNanowires/bc_paper_node_removal_figure1.8.out"
#SBATCH --error="/home/amwt/TPV/2DNanowires/bc_paper_node_removal_figure1.8.err"

source /home/amwt/.bashrc
module load intel
conda activate envnw
mpirun python /home/amwt/TPV/2DNanowires/bc_paper_node_removal_figure.py 1.8
