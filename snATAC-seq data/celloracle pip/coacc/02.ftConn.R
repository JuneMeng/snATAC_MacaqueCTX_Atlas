library(Matrix)
library(fitdistrplus)
library(ggplot2)
library(argparse)
parser <- ArgumentParser(description='fit model for filtering peak link')
parser$add_argument('-T', '--type', help='input subclass type')
parser$add_argument('-O', '--path', help='output directory path')
args <- parser$parse_args()

ctype <- args$type
ctype <- gsub('/','_',ctype)

peak_in <- read.table(paste0('/cluster/share/atac_group/mafas5/chen_ws/CREs/subclass_cre/','SST.UnionCRE.bed'),sep='\t')[,1:3]
peakCord <- paste(peak_in[,1],peak_in[,2],peak_in[,3],sep="_")

check_conn <- function(raw.conn){
  idx <- intersect(
    which(raw.conn$Peak1 %in% peakCord),
    which(raw.conn$Peak2 %in% peakCord))
  conns.f <- raw.conn[idx,]
  peak.pairs <- with(conns.f, paste(Peak1, Peak2, sep = "."))
  conns.f <- conns.f[!duplicated(peak.pairs), ]
  return(conns.f)
}

raw.conns <- readRDS(paste0('/cluster/share/atac_group/mafas5/chen_ws/SST/cicero/',ctype,'-union_cicero_connections.Rds'))
raw.conns <- raw.conns[!is.na(raw.conns$coaccess), ]
conns <- check_conn(raw.conns)

raw.shuf.conns <- readRDS(paste0('/cluster/share/atac_group/mafas5/chen_ws/SST/cicero/',ctype,'-shuffle_cicero_connections.Rds'))
raw.shuf.conns <- raw.shuf.conns[!is.na(raw.shuf.conns$coaccess), ]
shuf.conns <- check_conn(raw.shuf.conns)

# * filter based on shuf signals
fitnorm.shuf <- fitdistrplus::fitdist(data = shuf.conns$coaccess, distr = "norm", keepdata = FALSE)
mean.shuf <- fitnorm.shuf$estimate[[1]]
std.shuf <- fitnorm.shuf$estimate[[2]]

opath='/cluster/share/atac_group/mafas5/chen_ws/SST/cicero/'
outFitDistFile=paste0(opath,ctype,'-fitConns.dist.pdf')
withr::with_pdf(new = outFitDistFile, code = {
  plot(fitnorm.shuf)
  hist(shuf.conns$coaccess, pch = 20, breaks = 25, prob = TRUE,
    main = paste(ctype, " shuffled co-accessibility scores"))
  dat <- rbind(data.frame(coaccess = shuf.conns$coaccess,class = rep("shuf", nrow(shuf.conns))),
               data.frame(coaccess = conns$coaccess,class = rep("conns", nrow(conns))))
  print(ggplot(dat, aes(x = coaccess, color = class)) + 
          geom_density() +
          theme_bw() +
          theme(axis.text.x = element_text(colour = "black", size = 8),
            axis.text.y = element_text(colour = "black", size = 8),
            axis.title.x = element_text(colour = "black", size = 8),
            axis.title.y = element_text(colour = "black", size = 8)))
})

para <- data.frame(
  group = ctype,
  metaCol = 'SubClass',
  meanShuf = as.numeric(mean.shuf),
  stdShuf = as.numeric(std.shuf)
)
write.table(para, paste0(opath,ctype,'-fitConns.para.txt'), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

cal_p <- function(x, mean, sd, lower = TRUE) {
  p <- 1 - pnorm(x, mean = mean, sd = sd, lower.tail = TRUE, log.p = FALSE)
  return(p)
}

testPval <- vapply(conns$coaccess, cal_p, mean = mean.shuf, sd = std.shuf, lower = TRUE,FUN.VALUE = 0.0)
testFDR <- p.adjust(testPval, method = "BH")
conns$nlog10p <- (-log10(testPval))
conns$FDR <- testFDR
conns$p <- testPval

write.table(conns,paste0(opath,ctype,'-fitConns.res.txt'),quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


