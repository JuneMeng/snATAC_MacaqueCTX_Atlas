name="MergeGLU.V20241022"
outdir = "/cluster/share/atac_group/PublishedData/HumanMouseMacaque/GLUintegrated"
setwd(outdir)
homoGeneDir = "/cluster/share/atac_group/PublishedData"
oneToOneOrthGeneTb = readr::read_tsv(paste0(homoGeneDir, "/mart_export.humanMacaqeMouse.oneToOneOrth.ensembl91.20220428.txt"))
head(oneToOneOrthGeneTb)
HumanRawSeurat=readRDS("/cluster/share/atac_group/PublishedData/renbin_AdultHumanBrain/snapfile/GLUHuman.All.seuratSCT.rds")
MouseRawSeurat=readRDS("/cluster/share/atac_group/PublishedData/renbin_AdultMouseBrain2023/GSE246791_wmb_SnapATAC2_anndata/GLUMouse23.rename.seurat.rds")
MacaqueRawSeurat=readRDS("/cluster/home/mengjuan/project/snATAC_snRNA/5.WB_mafas5/5.merge2macaque/OutPlot/Anno.GLU.snATAC.rds")

#DownSampling
#The species with the largest number of clusters under a given subclass was allowed a maximum of 200 nuclei per cluster. The remaining species then split this theoretical maximum (200 nuclei multiplied by the maximum number of clusters under the subclass) evenly across their clusters. 
tmp=table(MacaqueRawSeurat$SubClass,MacaqueRawSeurat$CellType)
apply(tmp,1,function(x) sum(x>0)) %>% summary
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # 1.0(L5)     3.0     4.0     6.5     8.0    20.0
DownsamplingSeuratV2<-function(seurat,size=2000){
meta=seurat@meta.data
meta$X=rownames(meta)
tmplst=split(meta,meta$SubClass)
cells=c()
for(j in 1:length(tmplst)){
tmp=split(tmplst[[j]]$X,tmplst[[j]]$CellType)
Num=ceiling(size/length(tmp))
for(i in 1:length(tmp)){
tmp_dat=tmp[[i]]
if(length(tmp_dat)>Num){
tmp_dat=tmp_dat[sample(1:length(tmp_dat),size=Num)]
}
cells=c(cells,tmp_dat)
}
}
meta_sub=meta[cells,]
matrix_sub=seurat@assays$RNA@counts[,rownames(meta_sub)]
print(dim(matrix_sub))
seurat_sub=CreateSeuratObject(matrix_sub,project = "snATAC",meta.data=meta_sub)
return(seurat_sub)
}
macaqueSeurat=DownsamplingSeuratV2(MacaqueRawSeurat,size=2000) #24411
mouseSeurat=DownsamplingSeuratV2(MouseRawSeurat,size=2600) #24602
humanSeurat=DownsamplingSeuratV2(HumanRawSeurat,size=2800) #24589
macaqueSeurat$CellType %>% table 
mouseSeurat$CellType %>% table 
humanSeurat$CellType %>% table 
macaqueSeurat$SubClass %>% table 
mouseSeurat$SubClass %>% table 
humanSeurat$SubClass %>% table 
#
length(intersect(oneToOneOrthGeneTb$humanGene, rownames(humanSeurat)))
length(intersect(oneToOneOrthGeneTb$macaqueGene, rownames(macaqueSeurat)))
length(intersect(oneToOneOrthGeneTb$mouseGene, rownames(mouseSeurat)))
comGeneTb = oneToOneOrthGeneTb[c("humanGene", "macaqueGene", "mouseGene")]
comGeneTb$homoGeneSymbol = comGeneTb$macaqueGene #映射到猴子的同源基因上？
comGeneTb = subset(
    comGeneTb,
    humanGene %in% rownames(humanSeurat[["RNA"]]@counts) &
    macaqueGene %in% rownames(macaqueSeurat[["RNA"]]@counts) &
    mouseGene %in% rownames(mouseSeurat[["RNA"]]@counts)
)
comGeneTb #14382
prepareHomoSeurat = function(seurat, originalGeneName, homoGeneName, keepMetaCol, commonMetaCol, speciesPrefix="macaque_", clusterCol="anno3") {
    homoCountMx = seurat[["RNA"]]@counts[originalGeneName,]
    rownames(homoCountMx) = homoGeneName
    colnames(homoCountMx) = paste0(speciesPrefix, colnames(homoCountMx))
    homoSeurat = CreateSeuratObject(counts=homoCountMx)
    tmpDf = seurat@meta.data[sub(speciesPrefix, "", rownames(homoSeurat@meta.data)), keepMetaCol]
    tmpDf$species = sub("(_|-)", "", speciesPrefix)
    tmpDf$clusterName = paste0(speciesPrefix, tmpDf[[clusterCol]])
    commonMetaCol = c(commonMetaCol, "species", "clusterName")
    colnames(tmpDf) = ifelse(
        colnames(tmpDf) %in% commonMetaCol, colnames(tmpDf),
        paste0(speciesPrefix, colnames(tmpDf))
    )
    rownames(tmpDf) = paste0(speciesPrefix, rownames(tmpDf))
    homoSeurat = AddMetaData(homoSeurat, metadata=tmpDf)
    return(homoSeurat)
}
humanSeurat$specie="human"
mouseSeurat$specie="mouse"
macaqueSeurat$specie="macaque"
macaqueHomoSeurat = prepareHomoSeurat(
    macaqueSeurat, originalGeneName=comGeneTb$macaqueGene, homoGeneName=comGeneTb$homoGeneSymbol,
    keepMetaCol=c("orig.ident","nCount_RNA","nFeature_RNA","Sample","brain_area","batch","specie","Class","SubClass","CellType"),
    commonMetaCol=c("orig.ident","nCount_RNA","nFeature_RNA","Sample","brain_area","batch","specie","Class","SubClass","CellType"),
    speciesPrefix="macaque_", clusterCol="SubClass"
)
mouseHomoSeurat = prepareHomoSeurat(
    mouseSeurat, originalGeneName=comGeneTb$mouseGene, homoGeneName=comGeneTb$homoGeneSymbol, # Attention, seurat already use human gene name
    keepMetaCol=c("orig.ident","nCount_RNA","nFeature_RNA","Sample","brain_area","batch","specie","Class","SubClass","CellType"),
    commonMetaCol=c("orig.ident","nCount_RNA","nFeature_RNA","Sample","brain_area","batch","specie","Class","SubClass","CellType"),
    speciesPrefix="mouse_", clusterCol="SubClass"
)
humanHomoSeurat = prepareHomoSeurat(
    humanSeurat, originalGeneName=comGeneTb$humanGene, homoGeneName=comGeneTb$homoGeneSymbol, # Attention, seurat already use human gene name
    keepMetaCol=c("orig.ident","nCount_RNA","nFeature_RNA","Sample","brain_area","batch","specie","Class","SubClass","CellType"),
    commonMetaCol=c("orig.ident","nCount_RNA","nFeature_RNA","Sample","brain_area","batch","specie","Class","SubClass","CellType"),
    speciesPrefix="human_", clusterCol="SubClass"
)
mouseHomoSeurat
macaqueHomoSeurat
humanHomoSeurat
#merge
combineSeuratList = list(
	human=humanHomoSeurat,
    macaque=macaqueHomoSeurat,
    mouse=mouseHomoSeurat
)
saveRDS(combineSeuratList, paste0(outdir, "/",name,".combineSeuratList.rds"))
combineSeuratList = lapply(combineSeuratList, function(seurat) {
	DefaultAssay(seurat)="RNA"
	seurat = NormalizeData(seurat, verbose = FALSE) %>%
		FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = combineSeuratList,nfeatures = 3000)
length(features)
out_name1=paste0(name,".cca")
topN=2000
features <- SelectIntegrationFeatures(object.list = combineSeuratList,nfeatures = topN)
length(features)
out_name1=paste0(name,".topN",topN,".cca")
anchors <- FindIntegrationAnchors(object.list = combineSeuratList,
	anchor.features=features,
	normalization.method="LogNormalize",
	reduction = "cca", #"cca", "rpca", "jpca", "rlsi"
	)
saveRDS(anchors, paste0(outdir,"/",out_name1,".integrateAnchors.rds"))
integrateSeurat = IntegrateData(anchorset = anchors, normalization.method = "LogNormalize",dims=1:30)
DefaultAssay(integrateSeurat)
integrateSeurat = ScaleData(integrateSeurat, verbose = FALSE)
integrateSeurat = RunPCA(integrateSeurat, npcs = 100, verbose = FALSE)
pdf(paste0(outdir,"/",out_name1,".integrateSeurat.PCA.pdf"))
options(repr.plot.width=7, repr.plot.height=4)
print(ElbowPlot(integrateSeurat, ndims = 100))
dev.off()
pdf(paste0(outdir,"/",out_name1,".pca.pdf"),w=6,h=6)
print(DimPlot(integrateSeurat, reduction ="pca",group.by="species",shuffle = TRUE,label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="pca",group.by="SubClass",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="pca",group.by="batch",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
dev.off()
PC=20 #
resolution=0.8
k=20
out_name2=paste0(out_name1,".PC",PC)
integrateSeurat = FindNeighbors(integrateSeurat, reduction = "pca", dims = 1:PC,k.param = k) %>%
    FindClusters(resolution = resolution, n.start = 10, algorithm=1) %>%
    RunUMAP(dims = 1:PC)
table(integrateSeurat$seurat_clusters)
meta=integrateSeurat@meta.data
saveRDS(meta, paste0(outdir,"/",out_name2,".integrateSeurat.meta.rds"))
pdf(paste0(outdir,"/",out_name2,".umapsplit.pdf"),w=12,h=6)
print(DimPlot(integrateSeurat, reduction ="umap",group.by="SubClass",split.by="species",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="Class",split.by="species",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="batch",split.by="species",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
dev.off()
pdf(paste0(outdir,"/",out_name2,".umapsplit.celltype.pdf"),w=20,h=10)
print(DimPlot(integrateSeurat, reduction ="umap",group.by="CellType",split.by="species",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1)+NoLegend())
dev.off()
pdf(paste0(outdir,"/",out_name2,".umap1.pdf"),w=6,h=6)
print(DimPlot(integrateSeurat, reduction ="umap",group.by="species",shuffle = TRUE,label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="SubClass",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="Class",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="seurat_clusters",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="batch",shuffle = TRUE,label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
dev.off()
pdf(paste0(outdir,"/",out_name2,".umap2.pdf"),w=25,h=15)
print(DimPlot(integrateSeurat, reduction ="umap",group.by="clusterName",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="CellType",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="brain_area",shuffle = TRUE,label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="Sample",shuffle = TRUE,label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
dev.off()

meta=integrateSeurat@meta.data
umap=integrateSeurat@reductions$umap@cell.embeddings
meta$UMAP_X=umap[rownames(meta),"UMAP_1"]
meta$UMAP_Y=umap[rownames(meta),"UMAP_2"]
meta$Iterative=out_name2
saveRDS(meta,file=paste0(outdir,"/",out_name2,".integrateSeurat.meta.rds"))
saveRDS(integrateSeurat,file=paste0(outdir,"/",out_name2,".integrateSeurat.rds"))