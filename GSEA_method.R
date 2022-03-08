suppressMessages(library(fgsea))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
len<-length(args)
if(len < 4 ){
  cat("\nRscript Seurat_Rds metabolism.pathway.db outdir prefix")
  quit()
}
Seurat_Rdsfile<-args[1]
metabolism.db.file<-args[2]
outdir<-args[3]
prefix<-args[4]
path<-getwd()
setwd(path)
dir.create(outdir,recursive=T)

sc.sub<-readRDS(Seurat_Rdsfile)
deg.all<-FindAllMarkers(sc.sub,logfc.threshold = 0,only.pos = T,test.use = "MAST")
pathway<-read.table(metabolism.db.file,sep = "\t",header = T,check.names = F)
pathway.name<-colnames(pathway)

pathway.list<-list()
for(i in c(1:length(pathway.name))){
  # i<-1
  mypath<-pathway[,i]
  mypath1<-mypath[mypath!=""]
  pathway.list<-list.append(pathway.list,mypath1)
}         
names(pathway.list)<-pathway.name

deg.all<-deg
fgseaRes.all<-as.numeric()
for(i in as.character(unique(deg.all$cluster))){
  deg.all1<-deg.all[deg.all$cluster==i,]
  
  S4table<-deg.all1
  
  gene_list = S4table$avg_logFC
  names(gene_list) = S4table$gene
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  
  fgseaRes <- fgsea(pathways = pathway.list,
                    stats    = gene_list,
                    minSize  = 0,
                    maxSize  = 5000)
  fgseaRes<-as.data.frame(fgseaRes)
  fgseaRes1<-fgseaRes[,c(1:7)]
  fgseaRes1<-fgseaRes1[which(fgseaRes1$NES>0),]
  fgseaRes1$cluster<-i
  fgseaRes2<-fgseaRes1[order(fgseaRes1$NES,decreasing = T),]
  fgseaRes.all<-rbind(fgseaRes.all,fgseaRes2)
}
fgseaRes.all$cluster<-as.character(fgseaRes.all$cluster)
dim(fgseaRes.all)
fgseaRes.all1<-fgseaRes.all[fgseaRes.all$pval<0.1&fgseaRes.all$NES>0,]
dim(fgseaRes.all1)
fgseaRes.all
p1<-ggplot(fgseaRes.all,aes(x=cluster,y=pathway,size=-log10(pval),colour= NES))+geom_point()+theme_bw()+
  scale_color_gradient( low = "#FFFFFF",high = "#CC0000") 

ggsave(paste0(outdir,prefix,"-GSEA.dotplot.png"),plot = p1,width = 10,height = 12)                                                  write.table(fgseaRes.all,paste0(outdir,prefix,"-GSEA.png"),sep = ",",col.names = T,row.names = F,quote = F)
