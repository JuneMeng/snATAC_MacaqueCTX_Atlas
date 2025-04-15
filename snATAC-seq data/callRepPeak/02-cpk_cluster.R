source('/cluster/home/chencheng/Mac_gaba/prepare.R')
library(ArchR)
library(parallel)
library(tidyverse)
library(GenomicRanges)
addArchRThreads(threads = 30)

library(argparse)
parser <- ArgumentParser(description='extarct bedfiles from ArchR project')
parser$add_argument('-G', '--group', help='specify a column in colData for pseudobulk')
args <- parser$parse_args()

macsPeaks <- function(input, outname, outdir){
  system(paste0("macs2 callpeak -t ", input, 
                " -f BED -n ", outname, 
                " --outdir ", outdir,
                " -g 2871842584 -q 0.01 -B --SPMR", 
                " --nomodel --keep-dup all --call-summits --nolambda --shift 75 --extsize 150"))
}

callpeak <- function(celltype, suffix){
  if(suffix %in% c('rep','pseudo')){
    if(suffix == 'rep') {directory='./bedfiles_rep/'} else{directory='./bedfiles_pseudo/'}
    
    bedfiles <- list.files(directory)
    bedfiles <- bedfiles[grepl(celltype, bedfiles)]
    
    bed1 <- bedfiles[grepl('1_blacklistrm.bed',bedfiles)]
    bed2 <- bedfiles[grepl('2_blacklistrm.bed',bedfiles)]
    
    macsPeaks(input = paste0(directory, bed1),outname = paste0(celltype,'-',suffix,'1'), outdir=paste0("./peaks/"))
    macsPeaks(input = paste0(directory, bed2),outname = paste0(celltype,'-',suffix,'2'), outdir=paste0("./peaks/"))
    
  }
  if(suffix == 'pool') {
    directory='./bedfiles/'
    bedfiles <- list.files(directory)
    bedfiles <- bedfiles[grepl(celltype, bedfiles)]
    macsPeaks(input = paste0(directory, bedfiles),outname = paste0(celltype,'-',suffix), outdir=paste0("./peaks/"))
  }
}


Groups_all <- readRDS('Groups_all_meta.Rds')
Groups_all <- Groups_all[Groups_all$Anno2st == args$group,]

runinput <- data.frame(celltype = rep(unique(Groups_all$CellType), each=3),
                       suffix = rep(c('pool','rep','pseudo'), length(unique(Groups_all$CellType))))

timestart<-Sys.time();

mcmapply(function(X,Y) callpeak(X,Y), 
         X=runinput$celltype, Y=runinput$suffix,  mc.cores=20)

timeend<-Sys.time()
runningtime<-timeend-timestart
message("finish peak calling ---- time consumed: ", runningtime)
