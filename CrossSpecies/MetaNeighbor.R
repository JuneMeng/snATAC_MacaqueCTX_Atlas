seurat=paste0(outdir,"/",out_name2,".integrateSeurat.rds")
library(MetaNeighbor)
library(SummarizedExperiment)
options(future.globals.maxSize=20000*1024^3)
data=seurat@assays$integrated@data 
out_name3=paste0(name,".integrateData.SubClass")
meta=seurat@meta.data
table(meta$clusterName)
meta$group=meta$SubClass
table(meta$group)
var_genes=seurat@assays$integrated@var.features
length(var_genes)
sce=SummarizedExperiment(assays=list(gene_matrix=data[,rownames(meta)]),colData=meta)
table(sce$species,sce$group)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
	dat = sce,
	study_id = sce$species,
	cell_type = sce$group,
	symmetric_output=T,
	fast_version = TRUE #
	)
dim(celltype_NV)
saveRDS(celltype_NV,file=paste0(outdir,"/",out_name3,".celltype_NV.rds"))
diff_data=celltype_NV
diff_data[is.na(diff_data)]=0
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
pdf(paste0(outdir,"/",out_name3,".MetaNeighborUS.pdf"),w=ceiling(nrow(diff_data)/6+3),h=ceiling(nrow(diff_data)/6+3))
p=gplots::heatmap.2(diff_data,
	margins=c(8,8),keysize=1,key.xlab="AUROC",key.title="NULL",
	trace = "none",density.info = "none",
	col = cols,breaks = breaks,
	offsetRow=0.1,offsetCol=0.1,cexRow = 0.7,cexCol = 0.7)
dev.off()
#resort?
tmp=diff_data[p$rowInd,p$colInd]
write.table(tmp,file=paste0(outdir,"/",out_name3,".celltype_NV.txt"),quote=F,sep="\t")
#
row_sort=read.table(paste0(outdir,"/",out_name3,".celltype_NV.txt"),header=TRUE,sep="\t")
row_sort=row_sort$X
diff_data2=diff_data[row_sort,rev(row_sort)]
anno=data.frame(group=row_sort)
anno$Specie=sapply(anno$group,function(x) unlist(strsplit(x,"[|]"))[1] )
rownames(anno)=anno$group
anno$group=NULL
anno_color=list(Specie=SpecieCols)
library(pheatmap)
heatmap_corbar=rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
bk <- seq(0.01,1,by=0.01);length(bk)
p2=pheatmap(diff_data2, scale="none",cluster_rows = F, cluster_cols = F,
         color =  colorRampPalette(heatmap_corbar)(100), border_color = NA,
		 annotation_colors=anno_color,
		 annotation_row=anno,annotation_col=anno,
		 legend_breaks=c(0,0.5,1),breaks=bk,
         fontsize = 8)
pdf(paste0(outdir,"/",out_name3,".MetaNeighborUS.sort.pdf"),w=ceiling(nrow(diff_data)/9+2),h=ceiling(nrow(diff_data)/9+2))
print(p2)
dev.off()