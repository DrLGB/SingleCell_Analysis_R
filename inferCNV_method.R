suppressMessages(library(infercnv))
suppressMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
len<-length(args)
if(len < 4 ){
  cat("\nRscript Seurat_Rds grch38.bed outdir prefix")
  quit()
}
Seurat_Rdsfile<-args[1]
db.grch38<-args[2]
outdir<-args[3]
prefix<-args[4]
path<-getwd()
setwd(path)
dir.create(outdir,recursive=T)

sc<-readRDS(Seurat_Rdsfile)

barcode.fibro<-rownames(sc@meta.data)[sc$seurat_clusters=="Fibroblast"]
sc.fibro<-sc[,barcode.fibro]

barcode.epithe<-rownames(sc@meta.data)[sc$seurat_clusters=="Epithelial"]
sc.epithe<-sc[,barcode.epithe]


sc.merge<-merge(sc.epithe, y = sc.fibro, add.cell.ids = c("Epithelial", "Fibroblast"), project = "CNV")

reference.bedfile<-db.grch38
referencebed<-read.table(reference.bedfile,sep = "\t",header = T,stringsAsFactors = F,check.names = F)
referencebed<-referencebed[,c(1:5)]

detect.gene<-rownames(sc.merge[['RNA']]@counts)
detect.gene1<-referencebed[match(detect.gene,referencebed$`symbol/transcriptID`,0L),]

merge.sc1<-sc.merge[detect.gene1$`symbol/transcriptID`,]
sc.count<-merge.sc1[['RNA']]@counts
sc.cell.annote<-data.frame(row.names(merge.sc1@meta.data), merge.sc1@meta.data$seurat_clusters,stringsAsFactors = F)

sc.cell.annote$merge.sc1.meta.data.seurat_clusters[grep("Fibroblast",sc.cell.annote$row.names.merge.sc1.meta.data.)]<-"Fibroblast"

sc.gene.annotate<-data.frame(symbol=detect.gene1$`symbol/transcriptID`,chr=paste0("chr",detect.gene1$chr),
                           start=detect.gene1$start,end=detect.gene1$end,stringsAsFactors = F)

write.table(sc.count,paste0(outdire,"/",prefix,"-counts.txt"),sep = "\t",col.names = T,row.names = T,quote = F)
write.table(sc.cell.annote,paste0(outdire,"/",prefix,"-cellannotate.txt"),sep = "\t",col.names = F,row.names = F,quote = F)
write.table(sc.gene.annotate,paste0(outdire,"/",prefix,"-geneannotate.txt"),sep = "\t",col.names = F,row.names = F,quote = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=countfile,
                                    annotations_file= cellannotatefile,
                                    delim="\t",
                                    gene_order_file=genefile,
                                    ref_group_names=controlcelltype)
