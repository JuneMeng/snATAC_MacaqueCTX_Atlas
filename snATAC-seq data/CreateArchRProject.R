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