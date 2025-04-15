#!/bin/bash

#SBATCH -J getbed
#SBATCH -N 2
#SBATCH -n 4
#SBATCH --cpus-per-task=30
#SBATCH --mem=50000
#SBATCH -e errput_%A_%a_%J
#SBATCH -p q1 
#SBATCH -o output_%A_%a_%J
#SBATCH --array=1-7

date
echo "START----"

mapfile -t cell_types < rerun.lst 
#clu.lst

current_cell=${cell_types[$SLURM_ARRAY_TASK_ID-1]}
echo ${current_cell}

#Rscript testslurm.R -G $current_cell
Rscript get_bed.R -G ${current_cell}

echo "=========="
echo "FINISH WORK----"
date

