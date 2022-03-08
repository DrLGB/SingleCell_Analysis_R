suppressMessages(library(CytoTRACE))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
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


sc.epi<-readRDS(Seurat_Rdsfile)
pheno.sc<-sc.epi@meta.data
count<-as.matrix(sc.epi[['RNA']]@counts)

results <- CytoTRACE(count,ncores = 6,subsamplesize=2000)

##plot
cluster.num<-pheno.sc$seurat_clusters
cluster.num<-as.character(cluster.num)
names(cluster.num)<-rownames(pheno.sc)
plotCytoTRACE(results, phenotype = cluster.num)

stemscore=data.frame(barcode=names(results$CytoTRACE),CytoTRACE=results$CytoTRACE,stringsAsFactors = F)

pheno.sc<-pheno.sc[match(stemscore$barcode,rownames(pheno.sc),0L),]

stemscore<-data.frame(stemscore,cluster=pheno.sc$seurat_clusters,SampleType=pheno.sc$Group,stringsAsFactors = F)
stemscore$cluster<-as.character(stemscore$cluster)


p1<-ggplot(stemscore,aes(x=cluster,y=CytoTRACE,colour=SampleType,Group=SampleType,shape=SampleType))+geom_point()+
  scale_shape_manual(values = c(1,3))+
  theme_bw()+scale_shape_manual(values = c(1,3))+geom_jitter()+
  scale_color_manual(values = c("#CC0000","#00DD00"))+labs(x="",y="Differentiation potential")
ggsave(paste0(out_dir,prefix,"_cytotrace.png"),plot = p1,width = 6,height = 5)
