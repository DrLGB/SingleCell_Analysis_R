# SingleCell_Analysis_R
R script for 10X SingleCell data analyzed by Seurat (include: QC、PCA、TSNE、UMAP、Cluster、FindAllMarkers、etc.)


Seurat_Single_Sample is used for analysing the single data from 10X cellranger.
Seurat_multi_Sample is used for integrating the results from Seurat_Single_Sample.

Usage:

      Rscript Seurat_Single_Sample.R 10X_matrix_dir output_dir  Sample_name
      
      Rscript Seurat_multi_Sample.R A.rds,B.rds,C.rds,D.rds  output_dir
      
      (All Parameters in the R scripts can be changed by user.)
      
      
Note:
    Before using this two R script, you need to install R(version >= 3.5) and some libraries(DoubletFinder、Seurat).
    
    R: https://www.r-project.org/
    BoubletFinder: https://github.com/chris-mcginnis-ucsf/DoubletFinder
    Seurat: https://github.com/satijalab/seurat
    
    
