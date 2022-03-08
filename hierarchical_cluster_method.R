suppressMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
len<-length(args)
if(len < 3 ){
  cat("\nRscript Seurat_Rds  outdir prefix")
  quit()
}
Seurat_Rdsfile<-args[1]
outdir<-args[2]
prefix<-args[3]
path<-getwd()
setwd(path)
dir.create(outdir,recursive=T)

sc.sub<-readRDS(Seurat_Rdsfile)
pheno.sub<-sc.sub@meta.data
exprs.sub<-sc.sub[['RNA']]@data
sc.meanexprs.by.sample <-sapply(split(rownames(pheno.sub),pheno.sub$sampleID),function(cellID){
  rowMeans(exprs.sub[,cellID])
})

hcc <- hclust(parallelDist::parDist(t(sc.meanexprs.by.sample),threads =4, method = "euclidean"),
              method = "ward.D2")
png(paste0(outdir,prefix,"_ward.D2_cluster.png"),width = 800, height = 500,pointsize = 20)
plot(hcc,xlab="",ylab="",cex = 1.2, hang = -1)
dev.off()
