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