library(ArchR)
library(Seurat)
library(tidyverse)
library(presto)
library(parallel)
library(chromVARmotifs)
library("GenomicRanges")
library("BSgenome.Mfascicularis.NCBI.5.0")
gene_file="/cluster/share/atac_group/mafas5/ref/MJ_creat_Gene_Macaca_fascicularis_5.0.91.2.2.V20230616.txt"
exon_file="/cluster/share/atac_group/mafas5/ref/MJ_creat_Exon_Macaca_fascicularis_5.0.91.2.2.V20230616.txt"
tss_file="/cluster/share/atac_group/mafas5/ref/MJ_creat_Transcript_Macaca_fascicularis_5.0.91.2.2.V20230616.txt"
backlist_file="/cluster/share/atac_group/mafas5/ref/macFas5.huada.backlist.bed"
gene=read.csv(gene_file,header = F,sep="\t",stringsAsFactors = FALSE)
colnames(gene)=c("chr","start","end","strand","gene_id","symbol")
gene$symbol=make.unique(gene$symbol) 
gene1=gene[-grep("MT",gene$chr),]
MF_gene=makeGRangesFromDataFrame(gene1,keep.extra.columns = T)
exon=read.csv(exon_file,header = F,sep="\t",stringsAsFactors = FALSE)
colnames(exon)=c("chr","start","end","strand","gene_id","symbol")
exon$symbol=make.unique(exon$symbol)
exon1=exon[-grep("MT",exon$chr),]
MF_exon=makeGRangesFromDataFrame(exon1,keep.extra.columns = T)
tss=read.csv(tss_file,header = F,sep="\t",stringsAsFactors = FALSE)
colnames(tss)=c("chr","start","end","strand","tx_id","tx_name")
tss1=tss[-grep("MT",tss$chr),]
MF_tss=makeGRangesFromDataFrame(tss1,keep.extra.columns = T)
MF_tss=resize(MF_tss, 1, fix = "start") 
backlist=read.table(backlist_file,header = F,stringsAsFactors = FALSE)
backlist=backlist[,1:3]
colnames(backlist)=c("chr","start","end")
MF_backlist=makeGRangesFromDataFrame(backlist,keep.extra.columns = T)
seqnames(BSgenome.Mfascicularis.NCBI.5.0)=gsub("MFA","chr",seqnames(BSgenome.Mfascicularis.NCBI.5.0))
seqnames(BSgenome.Mfascicularis.NCBI.5.0)
genomeAnnotation <- createGenomeAnnotation(
	BSgenome.Mfascicularis.NCBI.5.0,
	blacklist = MF_backlist,
	)
geneAnnotation <- createGeneAnnotation(
  TSS = MF_tss,
  exons = MF_exon,
  genes = MF_gene
)
geneAnnotation$genes
geneAnnotation$exons
geneAnnotation$TSS

directory="/cluster/share/atac_group/tmpdata/files/fragment/"
fragment_file=dir(directory,pattern="*gz$")
name=sapply(fragment_file,function(x) unlist(strsplit(x,"[.]"))[1])
fragment_file=paste0(directory,fragment_file)
names(fragment_file)=name
length(fragment_file)
head(fragment_file)
uesthreads=35
addArchRThreads(threads = uesthreads)
ArrowFiles <- createArrowFiles(
  inputFiles = fragment_file,
  sampleNames = names(fragment_file),
  filterTSS = 3,minFrags = 1000, 
  filterFrags = 1000, 
  promoterRegion = c(2000, 100),
  addTileMat = TRUE,
  gsubExpression= ":.*",
  excludeChr = c("chrM", "chrY"),
  addGeneScoreMat = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  offsetPlus = 0,
  offsetMinus = 0,
  force=T,
  nChunk = 5,threads = getArchRThreads(),
  subThreading=FALSE
)
ArrowFiles
arrow_path=ArrowFiles
doubScores <- addDoubletScores(
    input = arrow_path, 
    k = 10, 
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)
getwd()
arrow_path=dir("./","*.arrow$")
arrow_path %>%length