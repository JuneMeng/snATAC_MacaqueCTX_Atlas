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

#metafile <- readRDS('/cluster/share/atac_group/mafas5/5.merge2macaque/Save-projHeme-Tile-MergeAnno/metaAll.anno.RefineClusterV1205.rds')
proj_all <- loadArchRProject('/cluster/share/atac_group/mafas5/cre/Save_noGABA_project')
#proj_all <- proj_all[proj_all$cellNames %in% rownames(metafile)]
#proj_all$CellType <- metafile[proj_all$cellNames,'CellType']
proj <- proj_all[proj_all$Anno2st == args$group]
rm(proj_all);gc()

ArrowFiles <- getArrowFiles(proj)

blacklist <- blacklist[1:20,]
blacklist <- makeGRangesFromDataFrame(blacklist)

chromLengths <- getChromLengths(proj)
chromLengths <- chromLengths [names(chromLengths) %ni% "chrX"]

Groups_all <- getCellColData(ArchRProj = proj, select = c("CellType","batch"), drop = TRUE)
Groups_all$cells <- rownames(Groups_all)
sample.list <- split(Groups_all$cells,Groups_all$batch)

pseu_rep <- lapply(sample.list, function(x) {
  nn <- rep('pseu_rep1',length(x))
  names(nn) <- x
  mm <- sample(x,floor(length(x)/2))
  nn[mm] <- 'pseu_rep2'
  return(nn)
})
pseu_rep <- unlist(pseu_rep)
names(pseu_rep) <- gsub('^M[1|2]\\.','',names(pseu_rep))
pseu_rep <- pseu_rep[Groups_all$cells]
Groups_all$pseudo_rep <- pseu_rep

# get beds ----------------------------------------------------------------
prepare_bed <- function(suffix, uselable, exchr=blacklist){
  availableChr <- paste0(unique(seqnames(exchr)))
  for(j in unique(Groups_all$CellType)){
    Groups_type <- Groups_all[Groups_all$CellType==j,]
    
    for(cluster in unique(Groups_type[,uselable])){
      tmp <- Groups_type[,uselable] == cluster
      Groups <- Groups_type[tmp,]
      
      usearrow <- unique(gsub('#.*','',Groups$cells))
      
      bedfile <- list()
      for(i in usearrow){
        cells = grep(i,Groups$cells,value=T)
        ArrowFiles.sub <- ArrowFiles[i]
        
        bed <- lapply(availableChr,function(cc){
          ArchR:::.getFragsFromArrow(ArrowFiles.sub, chr = cc, out = "GRanges", cellNames = cells)
        })
        bed <- do.call(c, bed)
        strand(bed) <- '+'
        bedfile[[i]] <- sort(bed)
        message('running arrow ',i)
      }
      
      all_beds <- as(bedfile, "GRangesList")
      all_beds <- unlist(all_beds)
      all_beds <- subsetByOverlaps(all_beds, exchr, invert=TRUE)
      
      if(suffix == 'pool') {directory='./bedfiles/'}
      if(suffix == 'rep') {directory='./bedfiles_rep/'}
      if(suffix == 'pseudo') {directory='./bedfiles_pseudo/'}
      j_oname <- gsub("/", "_", j)
      cluster_oname <- gsub("/", "_", cluster)
      rtracklayer::export(all_beds, con=paste0(directory,j_oname,'-',suffix,'-',cluster_oname,"_blacklistrm.bed"), format="bed")
      message("finished ##", paste0(j,'-',cluster))
    }
  }
}

suffix_name=c('pool','rep','pseudo')
label_name=c('CellType','batch','pseudo_rep')

timestart<-Sys.time();

mcmapply(function(X,Y) prepare_bed(X,Y, exchr=blacklist),
    X=suffix_name, Y=label_name, mc.cores=30)

timeend<-Sys.time()

runningtime<-timeend-timestart
message("prepare bed for #pool #rep #pseudo_rep ---- time consumed: ", timeend)

