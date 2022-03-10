suppressMessages(library(SCENIC))
suppressMessages(library(SCopeLoomR))
suppressMessages(library(Seurat))
suppressMessages(library(SCENIC))
suppressMessages(library(pheatmap))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
len<-length(args)

if(len < 4 ){
  cat("\nRscript Seurat_Rds  hg38_cistarget_db output_dir prefix")
  quit()
}
Seurat_Rdsfile<-args[1]
db_path<-args[2]
outdir<-args[3]
prefix<-args[4]

path<-getwd()
setwd(path)

dir.create(outdir,recursive=T)


sc<-readRDS(Seurat_Rdsfile)

exprMat <- as.matrix(sc[['RNA']]@counts)
cellInfo <- sc@meta.data

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

### Initialize settings 初始设置，导入评分数据库

cellInfo <- data.frame(cellInfo)
cellInfo$seurat_clusters<-as.character(cellInfo$seurat_clusters)
cluster.num <- "seurat_clusters"
colnames(cellInfo)[which(colnames(cellInfo)==cluster.num)] <- "CellType"

colVars <- list(CellType=allcolour[1:length( unique(cellInfo$CellType) )])
names(colVars$CellType)<-sort(unique(cellInfo$CellType))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]


hg38_dbs <- list('500bp'= 'hg38_cistarget_refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather','10kb' = 'hg38_cistarget_refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
db_mcVersion <- 'v9'
db_path <- db_path
scenicOptions <- initializeScenic(org = 'hgnc', dbDir = db_path, dbs = hg38_dbs, datasetTitle='SCENIC', nCores=6)

saveRDS(colVars, file="int/colVars.Rds")
saveRDS(cellInfo1, file="int/cellInfo.Rds")

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds")

### 共表达网络
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
scenicOptions <- initializeScenic(org = 'hgnc', dbDir = db_path, dbs = hg38_dbs, datasetTitle='SCENIC', nCores=1)
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


regulon<-readRDS(paste0("int/4.2_binaryRegulonActivity_nonDupl.Rds"))

binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
p1<-pheatmap::pheatmap(regulonActivity_byCellType_Binarized, 
                       breaks=seq(0, 1, length.out = 100),color = colorRampPalette(c("white","pink","red"))(100),
                       treeheight_row=10, treeheight_col=10, border_color=NA)
ggsave(paste0(outdir,prefix,"_binaryRegulonActivity_by_cluster.png"),plot = p1,width = 6,height = 6)
