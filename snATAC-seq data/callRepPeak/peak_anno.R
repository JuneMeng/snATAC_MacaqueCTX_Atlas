nn.cre <- lapply(nn,\(x){read.table(paste0(path,x),sep='\t')})
nn.cre <- do.call(rbind,nn.cre)
nn.cre$name <- paste(nn.cre$V1,nn.cre$V2,nn.cre$V3,sep='_')
nn.cre <- nn.cre[!duplicated(nn.cre$name),]

gaba.cre <- lapply(gaba,\(x){read.table(paste0(path,x),sep='\t')})
gaba.cre <- do.call(rbind,gaba.cre)
gaba.cre$name <- paste(gaba.cre$V1,gaba.cre$V2,gaba.cre$V3,sep='_')
dim(gaba.cre)
gaba.cre <- gaba.cre[!duplicated(gaba.cre$name),]
dim(gaba.cre)

glu.cre <- lapply(glu,\(x){read.table(paste0(path,x),sep='\t')})
glu.cre <- do.call(rbind,glu.cre)
glu.cre$name <- paste(glu.cre$V1,glu.cre$V2,glu.cre$V3,sep='_')
dim(glu.cre)
glu.cre <- glu.cre[!duplicated(glu.cre$name),]
dim(glu.cre)

peak.list <- lapply(file,\(x){read.table(paste0(path,x),sep='\t',header=T)})
peak.list <- lapply(peak.list, \(x) {makeGRangesFromDataFrame(x, keep.extra.columns=T)})
merged.gr <- Reduce(c, peak.list)
merged.gr <- nonOverlappingGR(merged.gr, by = "spm", decreasing = TRUE)


library(ChIPseeker)
library(GenomicFeatures)

gtf_file="/cluster/share/atac_group/mafas5/ref/Macaca_fascicularis_5.0.91.2.2_huada.gtf"

gtf_file="/cluster/share/atac_group/mafas5/chen_ws/Macaca_mulatta.Mmul_10.112.gtf"

gtf_file="/cluster/share/atac_group/human/ref/Homo_sapiens.GRCh38.110_MJ.gtf"

gtf_file="/cluster/share/atac_group/mouse/ref/Mus_musculus.GRCm38.93_MJ.gtf"

txdb <- makeTxDbFromGFF(file = gtf_file, format="gtf", dataSource="Ensembl", organism="Homo sapiens")

load(file = 'SA_maca_class_CRE.Rdata')

colnames(glu.cre) <- c('seqname','start','end','name')
colnames(gaba.cre) <- c('seqname','start','end','name')
colnames(nn.cre) <- c('seqname','start','end','name')

data=makeGRangesFromDataFrame(glu.cre,keep.extra.columns = T)
peak_anno <- annotatePeak(data,
                          tssRegion = c(-1000, 1000),
                          TxDb = txdb,
                          assignGenomicAnnotation = TRUE,
                          genomicAnnotationPriority = c("Promoter","5UTR", "3UTR", "Exon","Intron","Intergenic"),
                          addFlankGeneInfo = TRUE,
                          flankDistance = 5000)
glu.anno <- peak_anno

data=makeGRangesFromDataFrame(gaba.cre,keep.extra.columns = T)
peak_anno <- annotatePeak(data,
                          tssRegion = c(-1000, 1000),
                          TxDb = txdb,
                          assignGenomicAnnotation = TRUE,
                          genomicAnnotationPriority = c("Promoter","5UTR", "3UTR", "Exon","Intron","Intergenic"),
                          addFlankGeneInfo = TRUE,
                          flankDistance = 5000)
gaba.anno <- peak_anno

data=makeGRangesFromDataFrame(nn.cre,keep.extra.columns = T)
peak_anno <- annotatePeak(data,
                          tssRegion = c(-1000, 1000),
                          TxDb = txdb,
                          assignGenomicAnnotation = TRUE,
                          genomicAnnotationPriority = c("Promoter","5UTR", "3UTR", "Exon","Intron","Intergenic"),
                          addFlankGeneInfo = TRUE,
                          flankDistance = 5000)
nn.anno <- peak_anno

save(glu.anno,gaba.anno,nn.anno, file ='SA_macaca_anno.Rdata')

df <- lapply(c(glu.anno,gaba.anno,nn.anno),\(x) {as.data.frame(getAnnoStat(x)$Frequency)}) 
df <- do.call(cbind,df)
colnames(df) <- paste0('MaM_SA.',c('Glu','GABA','NN'))
df <- rbind(df[1:3,],colSums(df[4:5,]),colSums(df[6:7,]),df[8:9,])

peak.ann <- cbind(peak.ann,df)

library(ggplot2)
library(ggsci)


df <- readRDS('all_anno_peak.Rds')
df$type <- rownames(df)

pl.df <- reshape2::melt(df,id.vars = 'type' ,value.name = 'frac')

pl.df$type <- factor(pl.df$type, levels = sort(unique(pl.df$type))[c(7,1,2,5,6,4,3)])

pdf('peak_anno.pdf',width=5,height=4)
ggplot(pl.df, aes
       (variable, weight = frac, fill = type)) +
  geom_bar(color = "black", width = .7, position = 'stack') +
  labs(y ='Percentage (%)', x = 'Class') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_fill_tron()
dev.off()
