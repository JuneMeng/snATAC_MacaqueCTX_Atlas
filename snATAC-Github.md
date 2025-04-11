## snATAC-seq data process

### Fragment2Arrow.R

~~~R
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
~~~

### CreateArchRProject.R

~~~R
uesthreads=33
projHeme1 <- ArchRProject(
  ArrowFiles = arrow_path, 
  outputDirectory = "Save-projHeme-Raw-V20240710",
  copyArrows = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  threads = getArchRThreads()  
)
projHeme1 #2909835
output_directory=getOutputDirectory(projHeme1)
output_directory
saveRDS(projHeme1,file=paste0(output_directory,"/projHeme.raw.rds"))
meta=projHeme1@cellColData %>% as.data.frame()
saveRDS(meta,file=paste0(output_directory,"/projHeme.raw.meta.rds"))
~~~



### QualityControl.R

~~~R
#By Library
lst=readRDS("/cluster/home/mengjuan/project/snATAC_snRNA/5.WB_mafas5/5.merge2macaque/snATAC.AnnoColor.V1105.rds")
MonkeyCols=lst$MonkeyCols
region_L0_cols=lst$region_L0_cols
outdir=getOutputDirectory(proj)
outdir
samples=proj$Sample %>% unique
FragmentSizes=list()
TSSEnrichment=list()
FeaturePlot=list()
for(i in 6:length(samples)){
sample=samples[i]
cells=rownames(proj@cellColData)[proj$Sample==sample]
print(paste0(sample," : ",length(cells)," cells"))
proj_sub=proj[cells,]
df1=plotFragmentSizes(proj_sub,groupBy = "SampleName",
	chromSizes = getChromSizes(proj),
	maxSize = 750,returnDF = TRUE,threads = getArchRThreads()
	)
df2=plotTSSEnrichment(proj_sub,groupBy = "SampleName",
	chromSizes = getChromSizes(proj),
	TSS = getTSS(proj),
    flank = 2000,norm = 100,smooth = 11,
	returnDF = TRUE,threads = getArchRThreads()
	)
df1$Sample=sample
df2$Sample=sample
FragmentSizes[[i]]=df1
TSSEnrichment[[i]]=df2
p1=plotFragmentSizes(proj_sub,groupBy = "SampleName",
	chromSizes = getChromSizes(proj),
	maxSize = 750,returnDF = FALSE,threads = getArchRThreads()
	)
p2=plotTSSEnrichment(proj_sub,groupBy = "SampleName",
	chromSizes = getChromSizes(proj),
	TSS = getTSS(proj),
    flank = 2000,norm = 100,smooth = 11,
	returnDF = FALSE,threads = getArchRThreads()
	)
p=cowplot::plot_grid(p1,p2)
FeaturePlot[[i]]=p
}
names(FragmentSizes)=samples
names(TSSEnrichment)=samples
saveRDS(FragmentSizes,paste0(outdir,"/Raw_Stat_BySample.FragmentSizes.rds"))
saveRDS(TSSEnrichment,paste0(outdir,"/Raw_Stat_BySample.TSSEnrichment.rds"))
saveRDS(FeaturePlot,paste0(outdir,"/Raw_Stat_BySample.FragSizes-TSSProfile.rds"))
p=cowplot::plot_grid(plotlist=FeaturePlot,ncol=1,labels=samples)
pdf(paste0(outdir,"/Raw_Stat_BySample-FragSizes-TSSProfile.pdf"),w=8,h=4*length(FeaturePlot))
print(p)
dev.off()

tmp=proj@cellColData %>% as.data.frame() %>%
	dplyr::select(Sample,SampleName,brain_area,MajorRegion,SubRegion,batch)

df1
df=merge(df1,tmp,id="Sample")
dim(df)
df1=FragmentSizes %>% purrr::reduce(rbind)
cols=ArchR::paletteDiscrete(sort(unique(df1$Sample)));cols
p1 <- ggplot(as.data.frame(df1), aes(fragmentSize, fragmentPercent,color = Sample)) + 
            geom_line(size = 1) + theme_ArchR() +
            xlab("ATAC-seq Fragment Size (bp)") + ylab("Percentage of Fragments") +
            scale_color_manual(values = cols) + 
			guides(color = "none")+
			scale_y_continuous(limits = c(0,max(df1$fragmentPercent) * 1.05), expand = c(0,0)) + 
			scale_x_continuous(limits = c(min(df1$fragmentSize),max(df1$fragmentSize)), expand = c(0, 0))
pdf(paste0(outdir,"/Raw_Stat_BySample.FragmentSizes.pdf"),h=5,w=5)
print(p1)
dev.off()

#By Monkey
lst=readRDS("/cluster/home/mengjuan/project/snATAC_snRNA/5.WB_mafas5/5.merge2macaque/snATAC.AnnoColor.V1105.rds")
MonkeyCols=lst$MonkeyCols
region_L0_cols=lst$region_L0_cols
proj=readRDS("./Save-projHeme-Raw-V202407/Save-ArchR-Project.rds")
outdir=getOutputDirectory(proj)
getChromSizes(proj)
df1=plotFragmentSizes(proj,groupBy = "batch",
	chromSizes = getChromSizes(proj),
	maxSize = 750,pal = MonkeyCols,
	returnDF = TRUE,threads = getArchRThreads()
	)
saveRDS(df1,paste0(outdir,"/Raw_Stat_ByMonkey.FragmentSizes.rds"))
p1 <- ggplot(as.data.frame(df1), aes(fragmentSize, fragmentPercent,
            color = group)) + geom_line(size = 1) + theme_ArchR() +
            xlab("ATAC-seq Fragment Size (bp)") + ylab("Percentage of Fragments") +
            scale_color_manual(values = MonkeyCols) + 
			scale_y_continuous(limits = c(0,max(df1$fragmentPercent) * 1.05), expand = c(0,0)) + 
			scale_x_continuous(limits = c(min(df1$fragmentSize),max(df1$fragmentSize)), expand = c(0, 0))
pdf(paste0(outdir,"/Raw_Stat_ByMonkey.FragmentSizes.pdf"),h=5,w=5)
print(p1)
dev.off()
tmp=split(rownames(proj@cellColData),proj$batch)
tmp=sapply(tmp,function(x) sample(x,size=50000,replace=FALSE))
cells=tmp%>% unlist %>% as.character
proj_sub=proj[cells,]
df2=plotTSSEnrichment(proj_sub,groupBy = "batch",
	chromSizes = getChromSizes(proj),
	TSS = getTSS(proj),
    flank = 2000,norm = 100,smooth = 11,
	returnDF = TRUE,threads = getArchRThreads()
	)
saveRDS(df2,paste0(outdir,"/Raw_Stat_ByMonkey.TSSEnrichment.rds"))
p2 <- ggplot(as.data.frame(df2), aes(x, smoothValue, color = group)) + 
			geom_line(size = 1) +theme_ArchR() + 
            xlab("Distance From Center (bp)") +ylab("Normalized Insertion Profile") +
            scale_color_manual(values = MonkeyCols) +
			guides(color = "none")+
            scale_y_continuous(limits = c(0, max(df2$smoothValue) *1.05), expand = c(0, 0)) + 
			scale_x_continuous(limits = c(min(df2$x), max(df2$x)), expand = c(0, 0))
pdf(paste0(outdir,"/Raw_Stat_ByMonkey.TSSEnrichment.pdf"),h=5,w=5)
print(p2)
dev.off()

#Raw_Stat-TSS-vs-Frags-density
df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
df %>%head
p=ggPoint(
    x = df[,1], y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.9999)) ) +
 geom_hline(yintercept = c(5,30), lty = "dashed") + 
 geom_vline(xintercept = 3, lty = "dashed")
# pdf(paste0(outdir,"/Raw_Stat-TSS-vs-Frags-density.pdf"),h=5,w=5)
# print(p)
# dev.off()
tiff(paste0(outdir,"/Raw_Stat-TSS-vs-Frags-density.tiff"),h=400,w=400)
print(p)
dev.off()

#QC metrics
projHeme1=proj
idxPass <- which(projHeme1$TSSEnrichment >= 5 & projHeme1$TSSEnrichment <= 30 &
projHeme1$PromoterRatio >= 0.05 & projHeme1$PromoterRatio <= 0.6 &
projHeme1$DoubletScore<=0.15 & projHeme1$log10nFrags>=3
)
idxPass %>% length 
cellsPass <- projHeme1$cellNames[idxPass]
length(cellsPass) 
projHeme2=projHeme1[cellsPass, ]
projHeme2 #1943019
table(projHeme2$batch)
 # M1      M2
 # 456001 1487018
~~~



### UnsupervisedClustering.R

~~~R
#Read project rds
projHeme=readRDS(proj_file)
projHeme$group=projHeme$SampleName
output_directory=getOutputDirectory(projHeme)

addArchRThreads(threads = 8)
options(future.globals.maxSize= 1000*1024^3)
projHeme <- addGroupCoverages(ArchRProj = projHeme, groupBy = "group",force = TRUE)
saveRDS(projHeme,file=paste0(output_directory,"/",name,".BySampleName.projHeme.rds"))
pathToMacs2 <- findMacs2() #~/.conda/envs/ArchR/bin/macs2
pathToMacs2
projHeme <- addReproduciblePeakSet(
    ArchRProj = projHeme, 
    groupBy = "group", 
	reproducibility = "2",
    pathToMacs2 = pathToMacs2,
	peaksPerCell = 500,maxPeaks = 500000,minCells = 25,
	geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
	promoterRegion = c(2000, 100),
	genomeSize = 2.8e9
)
peakset=getPeakSet(projHeme)
length(peakset)
peakset$peakType %>% table
head(peakset)
saveRDS(peakset,file=paste0(output_directory,"/",name,".BySample.peakset.rds"))
projHeme <- addPeakMatrix(projHeme) 
getAvailableMatrices(projHeme)
saveRDS(projHeme,file=paste0(output_directory,"/",name,".BySample.projHeme.rds"))

resolution=3
maxClusters=20
min_dist=0.3 
n_neighbors=30
varFeatures = 50000 
varFeatures/length(getPeakSet(projHeme))
addArchRThreads(threads = 10)
options(future.globals.maxSize= 100000*1024^2)
t1=date()
projHeme <- addIterativeLSI(
  ArchRProj = projHeme,
  useMatrix = "PeakMatrix",
  name = "IterativeLSI",
  iterations = 4,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(resolution),
    sampleCells = 10000,
    n.start = 15
  ),
  seed=123,
  sampleCellsPre=10000,
  nPlot = 10000,
  projectCellsPre=FALSE,sampleCellsFinal=50000,
  varFeatures = varFeatures, #25000 50000
  dimsToUse = 1:30,
  UMAPParams = list(n_neighbors = n_neighbors, min_dist = min_dist, metric = "cosine", 
		verbose =FALSE, fast_sgd = TRUE),
  force = TRUE
)
saveRDS(projHeme,file=paste0(output_directory,"/",name,".projHeme.IterativeLSI.rds"))
projHeme <- addHarmony(
    ArchRProj = projHeme,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "SampleName",
	force = TRUE
)
projHeme <- addUMAP(
    ArchRProj = projHeme, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = n_neighbors, 
    minDist = min_dist, 
    metric = "cosine",
	force = TRUE
)
options(future.globals.maxSize= 100000*1024^2)
projHeme <- addClusters(
    input = projHeme,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "HarmonyClustersSub",
	maxClusters = maxClusters,seed=1,
    resolution = resolution,
	force = TRUE
)
saveRDS(projHeme,file=paste0(output_directory,"/",name,".projHeme.HarmonySampleName.rds"))
p2 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name = "batch", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name = "HarmonyClustersSub", embedding = "UMAPHarmony")
p=cowplot::plot_grid(p2,p4,nrow=1,ncol=2)
ggsave(paste0(output_directory,"/",name,".UmapHarmony.cluster.png"),plot=p,h=600,w=1200,units='mm',limitsize = FALSE)
p1 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name = "brain_area", embedding = "UMAPHarmony")
ggsave(paste0(output_directory,"/",name,".UmapHarmony.brain_area.png"),plot=p1,h=1200,w=600,units='mm',limitsize = FALSE)
p1 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
ggsave(paste0(output_directory,"/",name,".UmapHarmony.Sample.png"),plot=p1,h=1200,w=600,units='mm',limitsize = FALSE)
~~~



