#!/bin/bash

input_dir=$1
out_dir=$2
cl=$3

echo $cl

POOLED_PEAK="${input_dir}/${cl}-pool_peaks.narrowPeak"

function reuse_pooled_peak ()
{
    local r=$1
    if [ ! -f $1 ]; then
        r=$2
    fi
    echo "$r"
}

if [ ! -f ${POOLED_PEAK} ]; then
    exit 1
fi

REP1_PEAK=$(reuse_pooled_peak "${input_dir}/${cl}-rep1_peaks.narrowPeak" ${POOLED_PEAK})
REP2_PEAK=$(reuse_pooled_peak "${input_dir}/${cl}-rep2_peaks.narrowPeak" ${POOLED_PEAK})
PSEUDO1_PEAK=$(reuse_pooled_peak "${input_dir}/${cl}-pseudo1_peaks.narrowPeak" ${POOLED_PEAK})
PSEUDO2_PEAK=$(reuse_pooled_peak "${input_dir}/${cl}-pseudo2_peaks.narrowPeak" ${POOLED_PEAK})

PooledInRep1AndRep2="${out_dir}/${cl}.PooledInRep1AndRep2.narrowPeak.gz"
PooledInPsRep1AndPsRep2="${out_dir}/${cl}.PooledInPsRep1AndPsRep2.narrowPeak.gz"
naivePeakList="${out_dir}/${cl}.naivePeakList.narrowPeak.gz"

### Naive overlap
# Find pooled peaks that overlap Rep1 and Rep2
# where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
bedtools intersect -wo -a ${POOLED_PEAK} -b ${REP1_PEAK} -nonamecheck \
    | awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
    | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | bedtools intersect -wo -a stdin -b ${REP2_PEAK} \
    | awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
    | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | gzip -c > ${PooledInRep1AndRep2}

bedtools intersect -wo -a ${POOLED_PEAK} -b ${PSEUDO1_PEAK} -nonamecheck \
    | awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
    | cut -f 1-10 | sort -k1,1 -k2,2n | uniq \
    | bedtools intersect -wo -a stdin -b ${PSEUDO2_PEAK} -nonamecheck \
    | awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
    | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | gzip -c > ${PooledInPsRep1AndPsRep2}

# Combine peak lists
zcat ${PooledInRep1AndRep2} ${PooledInPsRep1AndPsRep2}  \
    | sort -k1,1 -k2,2n | uniq \
    | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
    | grep -P 'chr[0-9XY]+(?!_)' | gzip -c  > ${naivePeakList}

# Get summit 
POOLED_SUMMIT="${input_dir}/${cl}.summits.bed"
cat ${POOLED_PEAK} | awk 'BEGIN{FS=OFS="\t"}{s1=$2+$10; s2=$2+$10+1}{print $1,s1,s2,$4,$9}' > ${POOLED_SUMMIT}

naiveSummitList="${out_dir}/${cl}.naiveSummitList.bed"
join -1 1 -2 4 <(zcat ${naivePeakList} | cut -f 4 | sort) <(sort -k4,4 ${POOLED_SUMMIT}) -t$'\t' \
| awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$1,$5}' | sort -k1,1 -k2,2n > ${naiveSummitList} 
