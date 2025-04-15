#!/bin/bash

#SBATCH -J cpk
#SBATCH -N 3
#SBATCH -n 27
#SBATCH --cpus-per-task=8
#SBATCH --mem=60000
#SBATCH -e errput_%A_%a_%J
#SBATCH -p q1 
#SBATCH -o output_%A_%a_%J
#SBATCH --array=1-27

date
echo "START----"

mapfile -t cell_types < subclass.lst 

current_cell=${cell_types[$SLURM_ARRAY_TASK_ID-1]}
echo ${current_cell}

Rscript 02-cpk_cluster.R -G ${current_cell}

echo "=========="
echo "FINISH WORK----"
date
