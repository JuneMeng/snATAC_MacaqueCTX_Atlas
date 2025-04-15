#!/bin/bash

#SBATCH -J annoConn
#SBATCH -N 1
#SBATCH -n 10 
#SBATCH --mem=80000
#SBATCH -e errput_annoConn_%J
#SBATCH -p q1 
#SBATCH -o output_annoConn_%J

date
echo "START----"
for i in `cat rerun.txt`;
do
  echo $1
  bash /cluster/share/atac_group/mafas5/chen_ws/SST/03.addAnno.sh $i
  echo "finish #############################################"
done 

echo "=========="
echo "FINISH WORK----"
date

