mm10ToHg38_CHAIN_FILE="/cluster/share/atac_group/PublishedData/liftover/mm10.hg38.all.chain.gz"
macFas5ToHg38_CHAIN_FILE="/cluster/share/atac_group/PublishedData/liftover/macFas5.hg38.all.chain.gz"
Hg38Tomm10_CHAIN_FILE="/cluster/share/atac_group/PublishedData/liftover/hg38ToMm10.over.chain.gz"
Hg38TomacFas5_CHAIN_FILE="/cluster/share/atac_group/PublishedData/liftover/hg38ToMacFas5.over.chain.gz"

HumanDir="/cluster/share/atac_group/PublishedData/renbin_AdultHumanBrain/CREs/SubClass/"
MacaqueDir="/cluster/share/atac_group/mafas5/creV2/FilterPeaks/SubClass/"
MouseDir="/cluster/share/atac_group/PublishedData/renbin_AdultMouseBrain2023/CREs/"
WorkDir="/cluster/share/atac_group/PublishedData/HumanMouseMacaque/CREV20241203"
cd $WorkDir
minMatch=0.3
mkdir OutPlot
#HumanName="ET";MacaqueName="L5";MouseName="L5_ET_CTX_Glut"; #examples
echo $HumanName,$MacaqueName,$MouseName
cd $WorkDir
HumanRawBed=${HumanDir}${HumanName}".bed"
MacaqueRawBed=${MacaqueDir}${MacaqueName}".UnionCRE.bed"
MouseRawBed=${MouseDir}${MouseName}".bed"
wc -l $HumanRawBed
wc -l $MacaqueRawBed
wc -l $MouseRawBed
awk -v OFS='\t' '{print $1,$2,$3,"Human:"$1":"$2"_"$3}' $HumanRawBed > HumanRaw.${HumanName}.bed
awk -v OFS='\t' '{print $1,$2,$3,"Macaque:"$1":"$2"_"$3}' $MacaqueRawBed > MacaqueRaw.${MacaqueName}.bed
awk -v OFS='\t' '{print $1,$2,$3,"Mouse:"$1":"$2"_"$3}' $MouseRawBed > MouseRaw.${MouseName}.bed
UNLIFT_FILE=MacaqueRaw.${MacaqueName}.bed
OUTPUT_FILE=macFas5ToHg38.${MacaqueName}.bed
UNMAP_FILE=macFas5ToHg38.${MacaqueName}.Unmapped.bed
time liftOver $UNLIFT_FILE $macFas5ToHg38_CHAIN_FILE $OUTPUT_FILE $UNMAP_FILE -minMatch=$minMatch
UNLIFT_FILE=MouseRaw.${MouseName}.bed
OUTPUT_FILE=mm10ToHg38.${MouseName}.bed
UNMAP_FILE=mm10ToHg38.${MouseName}.Unmapped.bed
time liftOver $UNLIFT_FILE $mm10ToHg38_CHAIN_FILE $OUTPUT_FILE $UNMAP_FILE -minMatch=$minMatch
UNLIFT_FILE=macFas5ToHg38.${MacaqueName}.bed
OUTPUT_FILE=macFas5ToHg38TomacFas5.${MacaqueName}.bed
UNMAP_FILE=macFas5ToHg38TomacFas5.${MacaqueName}.Unmapped.bed
time liftOver $UNLIFT_FILE $Hg38TomacFas5_CHAIN_FILE $OUTPUT_FILE $UNMAP_FILE -minMatch=$minMatch
wc -l $UNLIFT_FILE
wc -l $OUTPUT_FILE
UNLIFT_FILE=mm10ToHg38.${MouseName}.bed
OUTPUT_FILE=mm10ToHg38Tomm10.${MouseName}.bed
UNMAP_FILE=mm10ToHg38Tomm10.${MouseName}.Unmapped.bed
time liftOver $UNLIFT_FILE $Hg38Tomm10_CHAIN_FILE $OUTPUT_FILE $UNMAP_FILE -minMatch=$minMatch
wc -l $UNLIFT_FILE
wc -l $OUTPUT_FILE
UNLIFT_FILE=HumanRaw.${HumanName}.bed
OUTPUT_FILE=HumanRawTomacFas5.${HumanName}.bed
UNMAP_FILE=HumanRawTomacFas5.${HumanName}.Unmapped.bed
time liftOver $UNLIFT_FILE $Hg38TomacFas5_CHAIN_FILE $OUTPUT_FILE $UNMAP_FILE -minMatch=$minMatch
wc -l $UNLIFT_FILE
wc -l $OUTPUT_FILE
UNLIFT_FILE=HumanRaw.${HumanName}.bed
OUTPUT_FILE=HumanRawTomm10.${HumanName}.bed
UNMAP_FILE=HumanRawTomm10.${HumanName}.Unmapped.bed
time liftOver $UNLIFT_FILE $Hg38Tomm10_CHAIN_FILE $OUTPUT_FILE $UNMAP_FILE -minMatch=$minMatch
wc -l $UNLIFT_FILE
wc -l $OUTPUT_FILE
#keep Re-mapping ID
Infile=macFas5ToHg38.${MacaqueName}.bed
IDfile=macFas5ToHg38TomacFas5.${MacaqueName}.bed
Outfile=macFas5ToHg38TomacFas5.${MacaqueName}.qc.bed
awk 'NR==FNR{a[$4]++}NR>FNR&&($4 in a)' $IDfile <(cat $Infile) > $Outfile
Infile=mm10ToHg38.${MouseName}.bed
IDfile=mm10ToHg38Tomm10.${MouseName}.bed
Outfile=mm10ToHg38Tomm10.${MouseName}.qc.bed
awk 'NR==FNR{a[$4]++}NR>FNR&&($4 in a)' $IDfile <(cat $Infile) > $Outfile
#
awk -F"\t" '{print $4}' HumanRawTomacFas5.${HumanName}.bed > ${HumanName}.tmpfile
awk -F"\t" '{print $4}' HumanRawTomm10.${HumanName}.bed >> ${HumanName}.tmpfile
cat ${HumanName}.tmpfile | sort | uniq > ${HumanName}.tmpfile2
wc -l ${HumanName}.tmpfile2
Infile=HumanRaw.${HumanName}.bed
Outfile=Human.${HumanName}.qc.bed
IDfile=${HumanName}.tmpfile2
awk 'NR==FNR{a[$1]++}NR>FNR&&($4 in a)' $IDfile <(cat $Infile) > $Outfile
cat mm10ToHg38Tomm10.${MouseName}.qc.bed |grep -v 'Un\|random\|chrM\|_alt'| sort -k1,1 -k2,2n | awk -F"\t" '$3-$2<1000{print}' > mm10ToHg38Tomm10.${MouseName}.qc2.bed
cat macFas5ToHg38TomacFas5.${MacaqueName}.qc.bed |grep -v 'Un\|random\|chrM\|_alt'| sort -k1,1 -k2,2n | awk -F"\t" '$3-$2<1000{print}' > macFas5ToHg38TomacFas5.${MacaqueName}.qc2.bed
cat Human.${HumanName}.qc.bed |grep -v 'Un\|random\|chrM\|_alt'| sort -k1,1 -k2,2n | awk -F"\t" '$3-$2<1000{print}' > Human.${HumanName}.qc2.bed
mkdir OutPlot/${MacaqueName}
cp macFas5ToHg38TomacFas5.${MacaqueName}.qc2.bed OutPlot/${MacaqueName}/1_macFas5ToHg38TomacFas5.${MacaqueName}.qc2.bed
cp Human.${HumanName}.qc2.bed OutPlot/${MacaqueName}/2_Human.${HumanName}.qc2.bed
cp mm10ToHg38Tomm10.${MouseName}.qc2.bed OutPlot/${MacaqueName}/3_mm10ToHg38Tomm10.${MouseName}.qc2.bed
cd OutPlot/${MacaqueName}
wc -l *bed
#./*bed:
intervene venn -i ./*bed \
			   --names=Macaque,Human,Mouse \
			   --bedtools-options f=0.1 \
			   --colors='#272E6A','#D51F26','#208A42' \
			   --output ./venn --save-overlaps