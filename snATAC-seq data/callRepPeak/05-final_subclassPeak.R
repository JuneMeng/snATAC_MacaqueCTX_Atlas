rule intersect_mergepeak:
    input:
        f"{outd_sum_mergepeak}/mba.whole.union.peak.bed"
    output:
        touch(f"{fdir}/{{cl}}_unionpeak.done")
    threads: 1
    resources:
        walltime = 1,
        queue = "glean",
        mail = "a",
        email = "debug.pie@gmail.com"
    conda: "sa22"
    shell:
        """
        bash {script_dir}/shell/intersect_mergepeak.sh \
           {input} \
           {outd_iter_mergepeak}/{wildcards.cl}.filterNfixed.peakset \
           {outd_unionpeak}/{wildcards.cl}.union.peak.bed
        """


#!/bin/bash

# echo $1
# echo $2
# echo $3
# echo done

if [ ! -f $1 ]; then
    echo "$1 does not exist."
    exit 1
fi

if [ ! -f $2 ]; then
    echo "$2 does not exist."
    exit 1
fi

echo "intersect merge peak for $3 ."

intersectBed -wa -a $1 -b <(sed -e "1d" $2) -nonamecheck \
   | sort -k1,1 -k2,2n | uniq > $3

echo "done."