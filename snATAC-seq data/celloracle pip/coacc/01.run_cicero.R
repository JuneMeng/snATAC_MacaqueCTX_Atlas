library(Matrix)
library(monocle3)
library(cicero)
library(stringr)
library(argparse)
library(GenomicRanges)
library(BSgenome.Mfascicularis.NCBI.5.0)
library(ArchR)
library(parallel)
library(tidyverse)

parser <- ArgumentParser(description='extarct bedfiles from ArchR project')
#parser$add_argument('-I', '--input', help='input peak matrix')
parser$add_argument('-C', '--ctype', help='cell type info')
parser$add_argument('-P', '--ptype', help='cell peak info')
parser$add_argument('-O', '--output', help='output path')
args <- parser$parse_args()


#inputfile <- list.files('/cluster/share/atac_group/mafas5/chen_ws/cicero/inputMtx',full.names=T)
#inputfile <- inputfile[grepl(input,inputfile)]
inputfile <- paste0('/cluster/share/atac_group/mafas5/chen_ws/SST/SST-peakMtx.Rds')
metafile <- paste0('/cluster/share/atac_group/mafas5/chen_ws/SST/SST-meta.Rds')
stopifnot(file.exists(inputfile))

indata <- readRDS(inputfile)
meta <- readRDS(metafile)
meta <- meta[meta$layer == args$ctype,]
indata <- indata[,rownames(meta)]
indata@x[indata@x > 0] <- 1
#proj <- loadArchRProject('/cluster/share/atac_group/mafas5/chen_ws/save-All-cortex-proj/')
#proj <- proj[colnames(indata)]

#cellinfo <- getCellColData(proj,select=c('Anno2st','CellType'))
#colnames(cellinfo) <- c('SubClass','CellType')
cellinfo = meta
colnames(cellinfo) <- c('SubClass','CellType','Layer')


if(args$ptype == 'shuffle'){
  idx <- sample(nrow(indata),dim(indata)[1],replace=F)
  oriname <- rownames(indata)
  indata <- indata[idx,]
  rownames(indata) <- oriname
}

cellinfo <- cellinfo[colnames(indata),]

peakinfo <- str_split(rownames(indata),'_',simplify=T)
colnames(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo <- as.data.frame(peakinfo)
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
rownames(peakinfo) <- peakinfo$site_name

input_cds <-  suppressWarnings(new_cell_data_set(indata,cell_metadata = cellinfo,gene_metadata = peakinfo))

input_cds <- monocle3::detect_genes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP

cicero_cds <- make_cicero_cds(input_cds, k = 50, reduced_coordinates = umap_coords)
chromosome_length <- read.table('/cluster/share/atac_group/mafas5/chen_ws/CREs/macaca_chr.genome')
#'/cluster/share/atac_group/mafas5/chen_ws/CREs/shuff_bed/macaca_chr.genome')
conns <- run_cicero(cicero_cds, chromosome_length, window = 500000)

#subclass = unique(cellinfo$SubClass)
path = args$output
input = args$ctype
ptype = args$ptype
input = gsub('/','_',input)
saveRDS(conns,paste0(path,input,'-',ptype,"_cicero_connections.Rds"))

all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = paste0(path,input,'-',ptype,"_all_peaks.csv"))
#write.csv(x = conns, file = paste0(path,subclass,'-',ptype,"_cicero_connections.csv"))
