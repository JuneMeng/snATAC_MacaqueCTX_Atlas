# merge peak -------------------------------------------------------------------------
source('/cluster/home/chencheng/Mac_gaba/prepare.R')
library(future.apply)
library(argparse)
parser <- ArgumentParser(description='merge peaks to union peakSet')
parser$add_argument('-I', '--input', help='cluster peak information')
args <- parser$parse_args()

workers <- as.integer(20)
options(future.globals.maxSize = 10e9)
plan(multicore, workers = workers)

inF <- args$input

outDir <- './mergePeaks'
if (!dir.exists(outDir)) {
  message("create outdir: ", outDir)
  dir.create(outDir, recursive = TRUE)
}
blacklist <- blacklist[1:20,]
blacklist <- makeGRangesFromDataFrame(blacklist)
chromLengths <- getChromLengths(proj)
chromLengths <- chromLengths [names(chromLengths) %ni% "chrX"]

read2gr <- function(bedF, label) {
  df <- fread(bedF, sep = "\t", header = F)
  colnames(df) <- c("chr", "start", "end", "name", "score")
  df$label <- label
  gr <- GRanges(
    df$chr,
    IRanges(df$start, df$end)
  )
  mcols(gr)$score <- df$score
  mcols(gr)$name <- df$name
  mcols(gr)$label <- df$label
  return(gr)
}

# extend summit to 500 bp
extendSummit <- function(gr, size=501) {
  gr <- resize(gr, width=size, fix="center")
  return(gr)
}

# filter blacklist
filter4blacklist <- function(gr, black_list.gr = blacklist) {
  
  idx <- queryHits(findOverlaps(gr, black_list.gr));
  if(length(idx) > 0) { gr <- gr[-idx] }
  return(gr)
  
}

# filter non-chromosome
filter4chrom <- function(gr, chromL=chromLengths) {
  chrom.gr <- GRanges(
    names(chromLengths),
    IRanges(0, chromLengths)
  )
  idx <- queryHits(
    findOverlaps(gr, chrom.gr, type="within")
  )
  if(length(idx) > 0) {
    gr <- gr[idx]
  }
  return(gr)
}

# filter N containing regions
filter4N <- function(gr, genome=BSgenome.Mfascicularis.NCBI.5.0) {
  #genome <- getBSgenome(genome)
  nucFreq <- BSgenome::alphabetFrequency(getSeq(genome, gr))
  mcols(gr)$GC <- round(
    rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  mcols(gr)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  gr[which(mcols(gr)$N < 0.001)] #Remove N Containing Peaks
  return(gr)
}

nonOverlappingGR <- function(
    gr = NULL, 
    by = "score", 
    decreasing = TRUE, 
    verbose = FALSE
) {
  stopifnot(by %in% colnames(mcols(gr)))
  #-----------
  # Cluster GRanges into islands using reduce and then select based on input
  #-----------
  .clusterGRanges <- function(gr = NULL, filter = TRUE, by = "score", decreasing = TRUE) {
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r, ignore.strand = TRUE)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  
  if(verbose) { message("Converging", appendLF = FALSE) }
  i <-  0
  grConverge <- gr
  while(length(grConverge) > 0) {
    if(verbose){ message(".", appendLF = FALSE) }
    i <-  i + 1
    grSelect <- .clusterGRanges(
      gr = grConverge, 
      filter = TRUE, 
      by = by, 
      decreasing = decreasing)
    
    grConverge <- subsetByOverlaps(
      grConverge,
      grSelect, 
      invert=TRUE, 
      ignore.strand = TRUE) #blacklist selected gr
    
    if(i == 1){ #if i=1 then set gr_all to clustered
      grAll <- grSelect
      
    }else{
      grAll <- c(grAll, grSelect)
    } 
  }
  message(sprintf("Converged after %s iterations!", i))
  
  if(verbose){
    message("\nSelected ", length(grAll), " from ", length(gr))
  }
  grAll <- sort(sortSeqlevels(grAll))
  return(grAll)
}

# normlize to score per million
norm2spm <- function(gr, by = "score") {
  mlogp <- mcols(gr)[,by]
  normmlogp <- 10^6 * mlogp / sum(mlogp)
  mcols(gr)$spm <- normmlogp
  return(gr)
}

summitF <- read.table(inF, sep="\t", header=T)
type.lst <- as.character(summitF$SubClass)
label.lst <- as.character(summitF$label)
file.lst <- as.character(summitF$filepath)

peak.list <- future_lapply(
  seq_along(file.lst), function(i) {
    message("working on summit set for... ", label.lst[i])
    p.gr <- read2gr(file.lst[i], label = label.lst[i])
    p.gr <- extendSummit(p.gr, size = 501)
    p.gr <- filter4blacklist(p.gr, black_list.gr=blacklist)
    p.gr <- filter4chrom(p.gr)
    p.gr <- filter4N(p.gr)
    p.gr <- nonOverlappingGR(p.gr, by = "score", decreasing = TRUE)
    p.gr <- norm2spm(p.gr, by = "score")
    outPeak <- as.data.frame(p.gr)
    outFname <- file.path(outDir,
                          paste0(label.lst[i], ".filterNfixed.peakset", sep = ""))
    message("write fixed & filtered peak set to: ", outFname)
    fwrite(outPeak, file = outFname,
           sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    return(p.gr)
  })
names(peak.list) <- type.lst

message("Finish filtered peak for each type.")

#peak.list <- list()
#for(i in inF$V1){
#	tmp <- read.table(paste0('/cluster/share/atac_group/mafas5/cre/mergePeaks/',i,'.filterNfixed.peakset'),sep='\t',header=T)
#	peak.list[[i]] <- makeGRangesFromDataFrame(tmp, keep.extra.columns=T)
#}

message("Merge all the peaks.")
#merged.gr <- do.call(c, peak.list)
merged.gr <- Reduce(c, peak.list)
merged.gr <- nonOverlappingGR(merged.gr, by = "spm", decreasing = TRUE)

outUnion <- as.data.frame(merged.gr)
outF='macaque.whole'
outfname <- file.path(outDir,paste0(outF, ".filteredNfixed.union.peakSet"))
message("save union peak set to: ", outfname)
fwrite(outUnion, file = outfname, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
message("Done.")

## filter based on spm

final.pk <- subset(outUnion, spm > 2)


