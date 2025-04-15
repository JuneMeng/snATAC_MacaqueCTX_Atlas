library(GenomicRanges)
library(BSgenome.Mfascicularis.NCBI.5.0)
library(ArchR)
library(parallel)
library(tidyverse)


geneAnnotation <- readRDS('macaca_geneAnnotation.Rds')

seqnames(BSgenome.Mfascicularis.NCBI.5.0)=gsub("MFA","chr",seqnames(BSgenome.Mfascicularis.NCBI.5.0))
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Mfascicularis.NCBI.5.0)
blacklist=read.table('/cluster/ref/macFas5.huada.backlist.bed',header = F,stringsAsFactors = FALSE)
blacklist=blacklist[,1:3]
colnames(blacklist)=c("chr","start","end")
MF_blacklist=makeGRangesFromDataFrame(blacklist,keep.extra.columns = T)
genomeAnnotation$blacklist=makeGRangesFromDataFrame(blacklist,keep.extra.columns = T)

