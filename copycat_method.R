suppressMessages(library(copykat))
suppressMessages(library(Seurat))
suppressMessages(library(grDevices))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

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

SC.sub<-readRDS(Seurat_Rdsfile)

exp.rawdata<-as.matrix(SC.sub[['RNA']]@counts)

pheno.sc<-SC.sub@meta.data

rawcout<-exp.rawdata
barcode.select<-colnames(rawcout)
copykat.test <- copykat(rawmat=exp.rawdata, 
                        id.type="S",
                        ngene.chr=5, 
                        win.size=25,
                        KS.cut=0.1, 
                        sam.name=prefix,
                        distance="euclidean",
                        norm.cell.names="",
                        n.cores=6)
pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)

copycat<-readRDS(paste0(outdir,prefix,"_copykat_clustering_results.rds"))
CNA.test<-read.table(paste0(outdir,prefix,"_copykat_CNA_results.txt"),sep = "\t",header = T,stringsAsFactors = F)
CNA.test.raw<-read.table(paste0(outdir,prefix,"_copykat_CNA_raw_results_gene_by_cell.txt"),sep = "\t",header = T,stringsAsFactors = F)
predi.test<-read.table(paste0(outdir,prefix,"_copykat_prediction.txt"),sep = "\t",header = T,stringsAsFactors = F)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
CNA.test<-CNA.test
pred.test<-predi.test
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

pheno.cnv<-pheno.sc[pred.test$cell.names,]
pred.test.df<-data.frame(pred.test,Group=pheno.cnv$Group,Cluster=as.character(pheno.cnv$seurat_clusters))

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

##sampletype
mycolor<-c("#DDAA00","#EEEE00")
sampletype.col<-mycolor[factor(pred.test.df$Group)]

##cluster###
mycolor1<-allcolour[1:7]
cluster.col<-mycolor1[factor(pred.test.df$Cluster)]

cells <- rbind(pred,sampletype.col)
cells1<-rbind(cells,cluster.col)

rownames(cells1)<-c("Prediction","SampleType","Cluster")

col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

pdf(paste0(outdir,prefix,".CNV.heatmap.pdf"),width = 8,height = 6,onefile = F)
heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =6, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D"),
          ColSideColors=chr1,RowSideColors=cells1,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",labRow = NULL,labCol = NULL,
          cexRow=0.01,cexCol=0.01,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1, bty="n")
legend("top", paste(names(table(pred.test.df$Group)),sep=""), pch=15,col=mycolor, lwd= 2,cex=1, bty="n")
legend("bottomright", paste(names(table(factor(pred.test.df$Cluster))),sep=""), pch=15,lwd= 2,col=mycolor1, cex=1, bty="n")
dev.off()
