#!/bin/bash

#SBATCH -J annda
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=500000
#SBATCH -e errput_annda_%J
#SBATCH -p q1
#SBATCH -o output_annda_%J

source activate celloracle_env

python /cluster/share/atac_group/mafas5/chen_ws/SST/celloracle/03.convertRNA.py

echo '=======================Finished============================'

