source('/cluster/share/atac_group/mafas5/chen_ws/prepare.R')
addArchRThreads(threads = 30)

proj <- loadArchRProject('/cluster/share/atac_group/mafas5/chen_ws/SST/Save-SST-project')
ArrowFiles <- getArrowFiles(proj)

library(parallel)

calPM <- function(ArrowFile, allCells=proj$cellNames, features=proj@peakSet, ceiling =4){
  sampleName <- names(ArrowFile)
  cellNames <- ArchR:::.availableCells(ArrowFile)
  cellNames <- cellNames[cellNames %in% allCells]
  strand(features) <- "*"
  uniqueChr <- as.character(unique(seqnames(features)@values))
  mat.lst <- list()
  
  for(z in seq_along(uniqueChr)){
    message('runing chr: ', uniqueChr[z])
	chr <- uniqueChr[z]
	featurez <- features[BiocGenerics::which(seqnames(features)==chr)]
	fragments <- ArchR:::.getFragsFromArrow(ArrowFile, chr = chr, out = "IRanges", cellNames = cellNames)
	tabFrags <- table(mcols(fragments)$RG)
	
	temp <- IRanges(start = start(fragments), width = 1)
	stopifnot(length(temp) == length(fragments))
	oleft <- findOverlaps(ranges(featurez), temp)
	oleft <- DataFrame(queryHits=Rle(queryHits(oleft)), subjectHits = subjectHits(oleft))
	
	temp <- IRanges(start = end(fragments), width = 1)
	stopifnot(length(temp) == length(fragments))
	oright <- findOverlaps(ranges(featurez), temp)
	oright <- DataFrame(queryHits=Rle(queryHits(oright)), subjectHits = subjectHits(oright))
	remove(temp)
	
	oleft$queryHits@values <- mcols(featurez)$idx[oleft$queryHits@values]
	oright$queryHits@values <- mcols(featurez)$idx[oright$queryHits@values]
	oleft$subjectHits <- as.integer(BiocGenerics::match(mcols(fragments)$RG[oleft$subjectHits], cellNames))
	oright$subjectHits <- as.integer(BiocGenerics::match(mcols(fragments)$RG[oright$subjectHits], cellNames))
	remove(fragments)
	
	#Create Sparse Matrix
    mat <- Matrix::sparseMatrix(
      i = c( oleft$queryHits, oright$queryHits ),
      j = c( oleft$subjectHits, oright$subjectHits ),
      x = rep(1, nrow(oleft) + nrow(oright)),
      dims = c(max(mcols(featurez)$idx), length(cellNames))
      )
    colnames(mat) <- cellNames
	mat@x[mat@x > ceiling] <- ceiling
	rownames(mat) <- paste(seqnames(featurez),paste(start(featurez),end(featurez),sep='_'),sep='_')
	mat.lst[[z]] <- mat
  }
  pm <- do.call("rbind", mat.lst)
  return(pm)
}


pkmtx.lst <- mclapply(ArrowFiles, function(x) calPM(x, allCells=proj$cellNames, features=proj@peakSet), mc.cores=30)
mtx <- do.call('cbind',pkmtx.lst)
meta=as.data.frame(getCellColData(proj,select=c('Anno2st','CellType')))
saveRDS(mtx,paste0('/cluster/share/atac_group/mafas5/chen_ws/SST/','SST-peakMtx.Rds'))
saveRDS(meta,paste0('/cluster/share/atac_group/mafas5/chen_ws/SST/','SST-meta.Rds'))

