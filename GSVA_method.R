suppressMessages(library(GSVA))
suppressMessages(library(msigdb))
suppressMessages(library(Seurat))
suppressMessages(library(limma))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
len<-length(args)
if(len < 4 ){
  cat("\nRscript Seurat_Rds pathway.db outdir prefix")
  quit()
}
Seurat_Rdsfile<-args[1]
pathway.db.file<-args[2]
outdir<-args[3]
prefix<-args[4]
path<-getwd()
setwd(path)
dir.create(outdir,recursive=T)

sc<-readRDS(Seurat_Rdsfile)

sc.exprs<-as.matrix(sc[['RNA']]@data)
sc.pheno<-sc@meta.data

pathway<-read.table(pathway.db.file,sep = "\t",header = T,check.names = F)
pathway.name<-colnames(pathway)

pathway.list<-list()
for(i in c(1:length(pathway.name))){
  # i<-1
  mypath<-pathway[,i]
  mypath1<-mypath[mypath!=""]
  pathway.list<-list.append(pathway.list,mypath1)
}         
names(pathway.list)<-pathway.name

res.data=gsva(sc.exprs,pathway.list,method="ssgsea",parallel.sz=10)

group<-data.frame(Group=sc.pheno$seurat_clusters,row.names =rownames(sc.pheno))
group<-group[colnames(res.data),]
names(group)<-colnames(res.data)
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(res.data)

#-compare----#
contrast.matrix<-makeContrasts("Tumor-Epi.Normal",levels=design)

##step1
fit <- lmFit(res.data,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
descript<-rownames(nrDEG)
splitfu<-function(x){
  myvect<-unlist(strsplit(x,"_"))
  descript1<-paste(tolower(myvect)[2:length(myvect)],collapse  = " " )
  return(descript1)
}
descript1<-sapply(descript, splitfu)
rownames(nrDEG)<-descript1

nrDEG1<-nrDEG[abs(nrDEG$t)>5,]
nrDEG1$Group[nrDEG1$t>0]<-"Tumor"
nrDEG1$Group[nrDEG1$t<0]<-"Epi.Normal"

nrDEG.up<-nrDEG1[nrDEG1$t>0,]
nrDEG.down<-nrDEG1[nrDEG1$t<0,]

nrDEG.all<-rbind(nrDEG.up[1:20,],nrDEG.down[1:20,])
nrDEG.all$Regulated<-c(rep("Tumor",20),rep("Epi.Normal",20))

nrDEG.all$Descript<-rownames(nrDEG.all)
nrDEG.all$Descript<-factor(nrDEG.all$Descript,levels = rev(nrDEG.all$Descript))
nrDEG.all$Regulated<-factor(nrDEG.all$Regulated,levels = c("Tumor","Epi.Normal"))
p.pathway<-ggplot(nrDEG.all)+ geom_bar(aes(x = Descript, y = t,  fill = factor(Regulated)),stat = "identity" ,
                                       position = position_dodge(), binwidth = 25)+
  theme_classic()+coord_flip()+  scale_fill_manual(values =c("#FF0000","#00BBFF") )+theme(legend.title = element_blank())
print(p.pathway)

ggsave(paste0(outdir,prefix,"DE_GSVAscore_limma_barplot.png"),plot =p.pathway,width = 15,height = 10 )
write.table(nrDEG1,paste0(outdir,prefix,"DE_GSVAscore_limma.csv"),sep = ",",col.names = T,row.names = T,quote = F)
