#!/bin/bash

#SBATCH -J mergetf
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=500000
#SBATCH -e errput_mergetf_%J
#SBATCH -p q1
#SBATCH -o output_mergetf_%J

source activate celloracle_env

python /cluster/share/atac_group/mafas5/chen_ws/SST/celloracle/02.mergeGimme.py

echo '=======================Finished============================'

