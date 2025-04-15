#!/bin/bash

#SBATCH -J tfscan
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=400000
#SBATCH -e errput_tfscan_%J
#SBATCH -p q1
#SBATCH -o output_tfscan_%J

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

for i in L2_3 L3_5 L5_6; do
  echo "Processing $i"  
  
  pdcbedpe_file="/cluster/share/atac_group/mafas5/chen_ws/SST/cicero/${i}.pdc.bedpe"
  
  if [ ! -f "$pdcbedpe_file" ]; then
    echo "Error: File not found: $pdcbedpe_file"
    continue  # 跳过当前循环，处理下一个变量
  fi

  python /cluster/share/atac_group/mafas5/chen_ws/SST/celloracle/01.runGimmemotifs.py \
  --pdcbedpe ${pdcbedpe_file} \
  --outdir /cluster/share/atac_group/mafas5/chen_ws/SST/celloracle \
  --threshold 10

  echo '=======================Finished============================'
done
