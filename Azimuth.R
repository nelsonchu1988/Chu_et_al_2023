library(Seurat)
library(SeuratDisk)

Convert('C://Users/rodri/Downloads/nelson_main_online.h5ad', dest = 'h5Seurat')

nelson <- LoadH5Seurat('C://Users/rodri/Downloads/nelson_main_online.h5Seurat', assays = 'RNA', meta.data = F)

library(Azimuth)

nelson_azimuth <- RunAzimuth(nelson, reference = "bonemarrowref")

DimPlot(nelson_azimuth, raster = F, group.by = 'predicted.celltype.l3')

SaveH5Seurat(nelson_azimuth, 'C://Bioinf/nelson_azimuth.h5Seurat')

DefaultAssay(nelson_azimuth) <- 'refAssay'
nelson_azimuth <- FindNeighbors(nelson_azimuth)
nelson_azimuth <- FindClusters(nelson_azimuth)

nelson_azimuth <- RenameIdents(nelson_azimuth, '23' = 'Endothelial', '17' = 'Osteo-lineage', '28' = 'Osteo-lineage', '27' = 'Smooth-muscle', '7' = 'Smooth-muscle')

# Assuming nelson_azimuth is your Seurat object

# 1. Identify cell IDs for the renamed clusters
cells_endothelial <- WhichCells(nelson_azimuth, idents = "Endothelial")
cells_osteolineage <- WhichCells(nelson_azimuth, idents = "Osteo-lineage")
cells_smoothmuscle <- WhichCells(nelson_azimuth, idents = "Smooth-muscle")

# Combine all cell IDs
all_renamed_cells <- c(cells_endothelial, cells_osteolineage, cells_smoothmuscle)

# 2. Extract current annotations (optional step if you want to check the current annotations before replacing)
current_annotations_l1 <- FetchData(nelson_azimuth, vars = "predicted.celltype.l1", cells = all_renamed_cells)
current_annotations_l2 <- FetchData(nelson_azimuth, vars = "predicted.celltype.l2", cells = all_renamed_cells)

# 3. Replace annotations in predicted.celltype.l1 and predicted.celltype.l2 with those from active.ident
# Get annotations from active.ident for the identified cell IDs
new_annotations <- nelson_azimuth@active.ident[all_renamed_cells]

# Replace in predicted.celltype.l1 and predicted.celltype.l2
nelson_azimuth@meta.data[all_renamed_cells, "predicted.celltype.l1"] <- new_annotations
nelson_azimuth@meta.data[all_renamed_cells, "predicted.celltype.l2"] <- new_annotations

# Note: Depending on your specific Seurat version and object structure, you might need to adjust the code.
# For example, if meta.data is a DataFrame, you might need to use `[[` instead of `$` to access columns.
DimPlot(nelson_azimuth, label = T, repel = T, raster = F, group.by = 'predicted.celltype.l2', reduction = 'umap')

nelson_azimuth <- SetIdent(nelson_azimuth, value = 'predicted.celltype.l1')
nelson_azimuth <- RenameIdents(nelson_azimuth, '1' = 'Endothelial', '2' = 'Osteo-lineage', '3' = 'Smooth-muscle')
nelson_azimuth$predicted.celltype.l1 <- nelson_azimuth@active.ident

nelson_azimuth <- SetIdent(nelson_azimuth, value = 'predicted.celltype.l2')
nelson_azimuth <- RenameIdents(nelson_azimuth, '1' = 'Endothelial', '2' = 'Osteo-lineage', '3' = 'Smooth-muscle')
nelson_azimuth$predicted.celltype.l2 <- nelson_azimuth@active.ident
DefaultAssay(nelson_azimuth) <- 'refAssay'
SaveH5Seurat(nelson_azimuth, 'C://Bioinf/nelson_azimuth_upd_EK_anno.h5Seurat', overwrite = T)
Convert('C://Bioinf/nelson_azimuth_upd_EK_anno.h5Seurat', dest = 'h5ad')

celltype_data <- data.frame(
  CellID = rownames(nelson_azimuth@meta.data),
  CellTypeL2 = nelson_azimuth$predicted.celltype.l2
)
celltype_data
write.csv(celltype_data, 'C://Bioinf/nelson_celltypel2.csv')
