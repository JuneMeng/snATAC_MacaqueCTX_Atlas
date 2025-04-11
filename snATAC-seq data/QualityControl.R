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