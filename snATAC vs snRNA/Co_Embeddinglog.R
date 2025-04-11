Co_Embeddinglog<-function(snATAC,snRNA,name){
library(parallel)
library(Seurat)
library(dplyr)
library(ggplot2)
out_name1=paste0(name,".log")
#1、downsampling from each cluster or each Subclass
#2、Integrate snATAC & snRNA
print(snATAC)
print(snRNA)
#choose the intersection genes
gene=intersect(rownames(snATAC),rownames(snRNA));length(gene)
combineSeuratList = list(
    ATAC=snATAC,
    RNA=snRNA
)
combineSeuratList = lapply(combineSeuratList, function(seurat) {
	seurat=CreateSeuratObject(seurat@assays$RNA@counts[gene,],
		meta.data=seurat@meta.data)
	seurat = NormalizeData(seurat, verbose = FALSE) %>%
		FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
})
#
HVG=list(ATAC=VariableFeatures(combineSeuratList$ATAC),
             RNA=VariableFeatures(combineSeuratList$RNA))
print("Intersections HVG：")
print(length(intersect(HVG$ATAC,HVG$RNA)))
print(length(union(HVG$ATAC,HVG$RNA)))
useHVG = Reduce(union, HVG);length(useHVG)
length(useHVG)
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = combineSeuratList,nfeatures = 3000)
length(features)
out_name1=paste0(out_name1,".cca")
#
anchors <- FindIntegrationAnchors(object.list = combineSeuratList,
	anchor.features=features,
	normalization.method="LogNormalize",
	reduction = "cca", 
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
print(DimPlot(integrateSeurat, reduction ="pca",group.by="type",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="pca",group.by="SubClass",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="pca",group.by="Class",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="pca",group.by="batch",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
dev.off()
PC=25
resolution=0.3
k=20
out_name2=paste0(out_name1,".PC",PC)
integrateSeurat = FindNeighbors(integrateSeurat, reduction = "pca", dims = 1:PC,k.param = k) %>%
    FindClusters(resolution = resolution, n.start = 10, algorithm=1) %>%
    RunUMAP(dims = 1:PC)
table(integrateSeurat$seurat_clusters)
meta=integrateSeurat@meta.data
saveRDS(meta, paste0(outdir,"/",out_name2,".integrateSeurat.meta.rds"))
pdf(paste0(outdir,"/",out_name2,".umapsplit.pdf"),w=12,h=6)
print(DimPlot(integrateSeurat, reduction ="umap",group.by="SubClass",split.by="type",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="Class",split.by="type",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="batch",split.by="type",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
dev.off()
pdf(paste0(outdir,"/",out_name2,".umapsplit.celltype.pdf"),w=20,h=10)
print(DimPlot(integrateSeurat, reduction ="umap",group.by="CellType",split.by="type",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1)+NoLegend())
dev.off()
pdf(paste0(outdir,"/",out_name2,".umap1.pdf"),w=6,h=6)
print(DimPlot(integrateSeurat, reduction ="umap",group.by="type",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="SubClass",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="Class",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="seurat_clusters",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="batch",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
dev.off()
pdf(paste0(outdir,"/",out_name2,".umap2.pdf"),w=25,h=15)
print(DimPlot(integrateSeurat, reduction ="umap",group.by="CellType",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="brain_area",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
print(DimPlot(integrateSeurat, reduction ="umap",group.by="Sample",label = TRUE, label.size=2,pt.size =0.5)+coord_fixed(ratio = 1))
dev.off()
#Save data
meta=integrateSeurat@meta.data
umap=integrateSeurat@reductions$umap@cell.embeddings
meta$UMAP_X=umap[rownames(meta),"UMAP_1"]
meta$UMAP_Y=umap[rownames(meta),"UMAP_2"]
meta$Iterative=out_name2
saveRDS(meta,file=paste0(outdir,"/",out_name2,".integrateSeurat.meta.rds"))
saveRDS(integrateSeurat,file=paste0(outdir,"/",out_name2,".integrateSeurat.rds"))
return(meta)
}