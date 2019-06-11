library(Seurat)

# -------------------------------------------------------------------------------------------------------------------------
# Using snRNA-seq mammalian brain
# Read data
gem <- scID::loadfast("~/ownCloud/DocSyncUoE/scID/TEST/snRNA_mammalian_brain/expression.txt")

cell_info <- read.delim("~/ownCloud/DocSyncUoE/scID/TEST/snRNA_mammalian_brain/metadata.txt", stringsAsFactors = F, row.names = 1)
labels <- cell_info$Cluster
names(labels) <- rownames(cell_info)
labels <- labels[-1]
rm(cell_info)

sobj <- CreateSeuratObject(counts = gem)
sobj <- FindVariableFeatures(sobj)
sobj <- ScaleData(sobj)
Idents(sobj) <- factor(labels)

markers <- FindAllMarkers(sobj, only.pos = FALSE, test.use = "MAST", logfc.threshold = 0.3)
write.table(markers, file="~/scID_manuscript_figures/data/Figure1/Hu_markers_26May_0.3logFC.txt", sep = "\t", quote = F)
