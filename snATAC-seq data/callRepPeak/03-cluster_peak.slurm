#!/bin/bash

#SBATCH -J repPeak
#SBATCH -N 3
#SBATCH -n 27
#SBATCH --cpus-per-task=8
#SBATCH --mem=50000
#SBATCH -e errput_%A_%a_%J
#SBATCH -p q1
#SBATCH -o output_%A_%a_%J
#SBATCH --array=1-27


date
echo "START----"

if [ ! -d "./peaks" ] || [ ! -d "./repPeaks" ]; then
  echo "Error: Directory ./peaks or ./repPeaks does not exist!"
  exit 1
fi

mapfile -t cell_types < subclass.lst

current_cell=${cell_types[$SLURM_ARRAY_TASK_ID-1]}
echo "Current cell type: ${current_cell}"

type=${current_cell//\//_}

bash 03-repPeak.sh ./peaks ./repPeaks "$type"

echo "=========="
echo "FINISH WORK----"
date
