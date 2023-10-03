#!/bin/bash
#
#SBATCH --job-name=gen1
#SBATCH --output=/home/amwt/TPV/2DNanowires/Data/evolution_under_node_removal/bc_paper_network_generator_001/array_%a.out
#SBATCH --error=/home/amwt/TPV/2DNanowires/Data/evolution_under_node_removal/bc_paper_network_generator_001/array_%a.err
#SBATCH --array=1-12
#SBATCH --mem-per-cpu=90000
#SBATCH --time=48:00:00


source /home/amwt/.bashrc
module load intel
conda activate envnw
echo slurm job id - $SLURM_JOB_ID
python /home/amwt/TPV/2DNanowires/bc_paper_network_generator.py /home/amwt/TPV/2DNanowires/Data/evolution_under_node_removal/bc_paper_network_generator_001/ $SLURM_ARRAY_TASK_ID
