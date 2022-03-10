# Step1 basied on analysis of SingleCells
R script for 10X SingleCell data analyzed by Seurat (include: QC、PCA、TSNE、UMAP、Cluster、FindAllMarkers、etc.)


Seurat_Single_Sample is used for analysing the single data from 10X cellranger.
Seurat_multi_Sample is used for integrating the results from Seurat_Single_Sample.

Usage:

      Rscript Seurat_Single_Sample.R 10X_matrix_dir output_dir  Sample_name
      
      Rscript Seurat_multi_Sample.R A.rds,B.rds,C.rds,D.rds  output_dir
      
      (All Parameters in the R scripts can be changed by user.)
      

# Step2 Advanced analysis of SingleCells
R scripts for 10x SingleCelll data include CytoTRACE, slingshot,SCENIC,copycat inferCNV,GSVA, GSEA analysis and hierarchical clustering etc.

CytoTRACE_method is used to predict differentiation potential of cells

The method of slingshot is used to  analysis  pseudotim trajectory

SCNICIC_method is used to analysis to transcription factor activity in single cells

copycat_method is used for CNV analysis to discrimination between malignant and nonmaligant cells 

inferCNV_method is another method to analysis CNV 

GSVA_method is used for pathway activity 

GSEA_method is used to analysis metabolism activity

hierarchical_cluster_method is used for hierarchical clustering


Usage: 
```
      Rscripts CytoTRACE_method.R 10x.seurat.rds outdir prefix
      Rscripts slingshot_method.R 10x.seurat.rds outdir prefix
      Rscripts SCENIC_method.R 10x.seurat.rds hg38_cistarget_db outdir prefix
      Rscripts copycat_method.R 10x.seurat.rds outdir prefix
      Rscripts inferCNV_method.R 10x.seurat.rds grch38.bed outdir prefix
      Rscripts GSVA_method.R 10x.seurat.rds pathway.db outdir prefix
      Rscripts GSEA_method.R 10x.seurat.rds metabolism.pathway.db outdir prefix
      Rscripts hierarchical_cluster_method.R 10x.seurat.rds  outdir prefix
      
```

Note:
    Before using this two R script, you need to install R(version >= 3.5) and some libraries(DoubletFinder、Seurat).
    
    R: https://www.r-project.org/
    BoubletFinder: https://github.com/chris-mcginnis-ucsf/DoubletFinder
    Seurat: https://github.com/satijalab/seurat
    
    
