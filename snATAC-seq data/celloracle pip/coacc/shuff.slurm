#!/bin/bash

#SBATCH -J shuffcc
#SBATCH -N 2
#SBATCH -n 8
#SBATCH --mem=100000
#SBATCH -e errput_%A_%a_%J
#SBATCH -p q1 
#SBATCH -o output_%A_%a_%J
#SBATCH --array=1-5

date
echo "START----"

mapfile -t cell_types < layer.txt 

current_cell=${cell_types[$SLURM_ARRAY_TASK_ID-1]}
echo ${current_cell}

Rscript /cluster/share/atac_group/mafas5/chen_ws/SST/01.run_cicero.R \
-C ${current_cell} \
-P shuffle \
-O /cluster/share/atac_group/mafas5/chen_ws/SST/cicero/ 

echo "================"
echo "FINISH WORK----"
date

