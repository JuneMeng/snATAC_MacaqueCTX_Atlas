#!/bin/bash

#SBATCH -J GRN
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=400000
#SBATCH -e errput_grn_%J
#SBATCH -p q1
#SBATCH -o output_grn_%J

source activate celloracle_env
TMPDIR='/cluster/share/atac_group/mafas5/chen_ws/SST/celloracle/tmpdir/'
NEW_CACHE=$TMPDIR/cache
mkdir -p $NEW_CACHE
if [ -z $XDG_CACHE_HOME ]; then
    XDG_CACHE_HOME=$HOME/.cache
fi
cp -r $XDG_CACHE_HOME/gimmemotifs $NEW_CACHE/
export XDG_CACHE_HOME=$NEW_CACHE
echo 'Using $XDG_CACHE_HOME for cache'

for i in L2_3 L3 L3_5 L5 L5_6; do
  echo "Processing $i"  
  
  python /cluster/share/atac_group/mafas5/chen_ws/SST/celloracle/04.runGRN.py \
  --tfdir /cluster/share/atac_group/mafas5/chen_ws/SST/celloracle \
  --scRNAdir /cluster/share/atac_group/mafas5/chen_ws/SST/celloracle/RNA/layer \
  --outdir /cluster/share/atac_group/mafas5/chen_ws/SST/celloracle/grn \
  --group $i \
  --njobs 4 \
  --baseGRN_all 1

  echo '=======================Finished============================'
done

