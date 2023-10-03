#!/bin/bash
#
#SBATCH --job-name=node_removal_1.8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32000
#SBATCH -t 48:00:00
#SBATCH --output="/home/amwt/TPV/2DNanowires/bc_paper_node_removal_figure1.8.out"
#SBATCH --error="/home/amwt/TPV/2DNanowires/bc_paper_node_removal_figure1.8.err"

source /home/amwt/.bashrc
conda activate envnw
srun python /home/amwt/TPV/2DNanowires/bc_paper_node_removal_figure.py 1.8
