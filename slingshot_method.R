suppressMessages(library(slingshot))
suppressMessages(library(destiny))
suppressMessages(library(scater))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
len<-length(args)

if(len < 3 ){
  cat("\nRscript Seurat_Rds  output_dir prefix")
  quit()
}
Seurat_Rdsfile<-args[1]
outdir<-args[2]
prefix<-args[3]
path<-getwd()
setwd(path)

dir.create(outdir,recursive=T)
sc.epi<-readRDS(Seurat_Rdsfile)
sc.scater.epi<-as.SingleCellExperiment(sc.epi)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

dm <- DiffusionMap(sc.scater.epi,k = 3)
dm.df<-as.data.frame(dm)
p1<-qplot(DC1, DC2, data = dm.df, colour = factor(seurat_clusters)) +theme_bw()+
  scale_color_manual(values = allcolour)+theme(legend.title = element_blank())

ggsave(paste0(outdir,prefix,"_destiny_DCplot.png"),plot = p1,width = 6,height = 6)

DC.meta<-data.frame(dm.df[,c(1:2)])
DC.meta1<-DC.meta[rownames(sc.epi@meta.data),]
dc.array<-array(DC.meta1,dim=c(1,2))
colnames(dc.array)<-c("DC1","DC2")
dc.array.mat<-as.matrix(dc.array) 

dc.dr <- CreateDimReducObject(
  embeddings = dc.array.mat,
  loadings =dc.array.mat,
  stdev = 0,
  key = "DC_",
  assay = "RNA"
)

sc.epi@reductions$DC<-dc.dr

#########slingshot analysis###########
#set start cell
sce <- slingshot(sc.scater.epi, clusterLabels = 'seurat_clusters', reducedDim = "DC",
                 allow.breaks = FALSE,start.clus=0)
# get the lineages:
lnes <- getLineages(reducedDim(sce,"DC"),
                    sce$seurat_clusters,start.clus=0,)

#------------------------------VISION--------------------#
suppressMessages(library(Polychrome))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(ggthemes))
suppressMessages(library(VISION))

slingshot_df <- data.frame(sce@colData)
slingshot_df<-cbind(slingshot_df,DC.meta1)
crv1 <- getCurves(lnes)
p2<-ggplot(slingshot_df, aes(x = DC1, y = DC2, 
                             colour = slingPseudotime_1)) +geom_point()+theme_bw()+
  scale_color_viridis()+labs(colour = "Pseudotime")
p3<-ggplot(slingshot_df, aes(x = DC1, y = DC2, 
                             colour = slingPseudotime_2)) +geom_point()+theme_bw()+
  scale_color_viridis()+labs(colour = "Pseudotime")

ggsave(paste0(outdir,prefix,"_destiny_DCplot_by_Pseudotime1.png"),plot = p2,width = 6,height = 6)
ggsave(paste0(outdir,prefix,"_destiny_DCplot_by_Pseudotime2.png"),plot = p3,width = 6,height = 6)
