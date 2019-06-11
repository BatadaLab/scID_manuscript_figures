#devtools::install_version(package = 'Seurat', version = package_version('2.3.0'))
library(Seurat)
library(dplyr)
library(pheatmap)
library(scran)
library(scmap)
library(ggsci)
library(cowplot)
library(scID)
source("~/ownCloud/DocSyncUoE/scID/scripts/resubmission_scripts_May2019/run_other_methods.R")

# Read data
ctrl.data <- read.table(file = "~/scID_manuscript_figures/data/Kang_data/immune_control_expression_matrix.txt.gz", sep = "\t")
ctrl.data.cpm <- counts_to_cpm(ctrl.data)

stim.data <- read.table(file = "~/scID_manuscript_figures/data/Kang_data/immune_stimulated_expression_matrix.txt.gz", sep = "\t")
stim.data.cpm <- counts_to_cpm(stim.data)

# -------------------------------------------------------------------------
# Cluster with Seurat v2 CCA
# ctrl <- CreateSeuratObject(raw.data = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
# ctrl@meta.data$stim <- "CTRL"
# ctrl <- FilterCells(ctrl, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
# ctrl <- NormalizeData(ctrl)
# ctrl <- ScaleData(ctrl, display.progress = F)
# # Set up stimulated object
# stim <- CreateSeuratObject(raw.data = stim.data, project = "IMMUNE_STIM", min.cells = 5)
# stim@meta.data$stim <- "STIM"
# stim <- FilterCells(stim, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
# stim <- NormalizeData(stim)
# stim <- ScaleData(stim, display.progress = F)
# 
# # Gene selection for input to CCA
# ctrl <- FindVariableGenes(ctrl, do.plot = F)
# stim <- FindVariableGenes(stim, do.plot = F)
# g.1 <- head(rownames(ctrl@hvg.info), 1000)
# g.2 <- head(rownames(stim@hvg.info), 1000)
# genes.use <- unique(c(g.1, g.2))
# genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
# genes.use <- intersect(genes.use, rownames(stim@scale.data))
# 
# immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30)
# 
# p3 <- MetageneBicorPlot(immune.combined, grouping.var = "stim", dims.eval = 1:30, display.progress = FALSE)
# 
# immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim", dims.align = 1:20)
# 
# immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
# immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)
# 
# p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "stim")
# p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
# plot_grid(p1, p2)
# 
# FeaturePlot(object = immune.combined, features.plot = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
# 
# new.ident <- c("CD14 Mono", "CD4 Naive T", "CD4 Memory T", "B", "CD16 Mono", "T activated", "CD8 T", "NK", "DC", "B activated", "Mk", "pDC", "Eryth")
# for (i in 0:12) {
#   immune.combined <- RenameIdent(object = immune.combined, old.ident.name = i, new.ident.name = new.ident[i + 1])
# }
# TSNEPlot(immune.combined, do.label = T, pt.size = 0.5)
# 
# saveRDS(immune.combined, file="~/scID_manuscript_figures/data/Kang_data/immuned_combined.rds")
immune.combined <- readRDS("~/scID_manuscript_figures/data/Kang_data/immuned_combined.rds")

# -------------------------------------------------------------------------
# Map stimulated cells using control as reference with scID
scID_res <- scid_multiclass(target_gem = stim.data.cpm, reference_gem = ctrl.data, normalize_reference = TRUE,
                            reference_clusters = immune.combined@ident[colnames(ctrl.data)], 
                            logFC = 0.5, only_pos = T, estimate_weights_from_target = T)
scID_labels <- scID_res$labels
table(scID_labels)
table(immune.combined@ident[colnames(ctrl.data)])

scID_ARI <- mclust::adjustedRandIndex(scID_labels[names(which(scID_labels != "unassigned"))], immune.combined@ident[names(which(scID_labels != "unassigned"))])
scID_VI <- mcclust::vi.dist(scID_labels[names(which(scID_labels != "unassigned"))], immune.combined@ident[names(which(scID_labels != "unassigned"))])

markers <- scID_res$markers[which(scID_res$markers$avg_logFC > 0), ]
markers_filt <- markers[which(markers$gene %in% rownames(stim.data.cpm)), ]

gem_avg <- matrix(NA, length(unique(immune.combined@ident)), length(unique(immune.combined@ident)))
for (i in 1:length(unique(immune.combined@ident))) {
  cells <- names(scID_labels)[which(scID_labels == unique(immune.combined@ident)[i])]
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(stim.data.cpm[markers_filt$gene, cells, drop=FALSE])
  } else {
    next
  }
  for (j in 1:length(unique(immune.combined@ident))) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == unique(immune.combined@ident)[j])]]))
  }
}
rownames(gem_avg) <- unique(immune.combined@ident)
colnames(gem_avg) <- unique(immune.combined@ident)

pheatmap(gem_avg, border="white", color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), cluster_rows = F,
         cluster_cols = F, show_colnames = T, fontsize_row = 3, border_color = F, show_rownames = T,
         scale = "row", fontsize_col = 3, width = 2, height = 2, legend = F,
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/New plots/Kang_figure/scID_heatmap.pdf")

# -------------------------------------------------------------------------
# Show expression of markers in original clusters from CCA
celltypes <- unique(immune.combined@ident)
gem_avg_ds <- matrix(NA, nrow = length(celltypes), ncol = length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- names(immune.combined@ident)[which(immune.combined@ident == celltypes[i])]
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(immune.combined@scale.data[intersect(markers$gene, rownames(immune.combined@scale.data)), cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg_ds[j,i] <- mean(na.omit(avg_exp[toupper(markers$gene[which(markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg_ds) <- celltypes
colnames(gem_avg_ds) <- celltypes

# Figure 2D left
pheatmap(gem_avg_ds, border="white", color = colorspace::diverge_hsv(50), cluster_rows = F,
         cluster_cols = F, show_colnames = T, fontsize_row = 3, border_color = F,
         show_rownames = T, fontsize_col = 3, scale = "row", width = 2, height = 2, legend = F,
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/New plots/Kang_figure/CCA_heatmap.pdf")


# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scmap (Fig 4e right)
common_genes <- intersect(rownames(stim.data.cpm), rownames(ctrl.data.cpm))

ann <- data.frame(cell_type1=immune.combined@ident[colnames(ctrl.data)])

# Create index
ref_sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(ctrl.data.cpm[common_genes, ])), colData = ann)
# use gene names as feature symbols
rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
# remove features with duplicated names
ref_sce <- ref_sce[!duplicated(rownames(ref_sce)), ]
ref_sce <- selectFeatures(ref_sce, suppress_plot = TRUE)#, n_features = 150)
ref_sce <- indexCluster(ref_sce)

# Create sce for testing data (SS2)
sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(stim.data.cpm[common_genes, ])))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]
sce <- selectFeatures(sce, suppress_plot = TRUE)#, n_features = 150)

scmapCluster_results <- scmapCluster(
  projection = sce, 
  index_list = list(
    DropSeq = metadata(ref_sce)$scmap_cluster_index
  )
)

scmap_labels <- scmapCluster_results$combined_labs
names(scmap_labels) <- colnames(stim.data.cpm)
table(scmap_labels)

scmap_ARI <- mclust::adjustedRandIndex(immune.combined@ident[names(which(scmap_labels != "unassigned"))], scmap_labels[names(which(scmap_labels != "unassigned"))])
scmap_VI <- mcclust::vi.dist(immune.combined@ident[names(which(scmap_labels != "unassigned"))], scmap_labels[names(which(scmap_labels != "unassigned"))])

gem_avg <- matrix(NA, length(unique(immune.combined@ident)), length(unique(immune.combined@ident)))
for (i in 1:length(unique(immune.combined@ident))) {
  cells <- names(scmap_labels)[which(scmap_labels == unique(immune.combined@ident)[i])]
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(stim.data.cpm[markers_filt$gene, cells, drop=FALSE])
  } else {
    next
  }
  for (j in 1:length(unique(immune.combined@ident))) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == unique(immune.combined@ident)[j])]]))
  }
}
rownames(gem_avg) <- unique(immune.combined@ident)
colnames(gem_avg) <- unique(immune.combined@ident)

pheatmap(gem_avg, border="white", color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), cluster_rows = F,
         cluster_cols = F, show_colnames = T, fontsize_row = 3, border_color = F,show_rownames = T,
         scale = "row", fontsize_col = 3, width = 2, height = 2, legend = F,
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/New plots/Kang_figure/scmap_heatmap.pdf")

# ---------------------------------------------------------------
# use CaSTLe

ctrl.cells <- intersect(names(immune.combined@ident), colnames(ctrl.data))
castle_labels <- runCastle(stim.data, ctrl.data[, ctrl.cells], immune.combined@ident[ctrl.cells])

# Get ARI and VI (variation of information) with marker based for different thresholds
castle_ARI <- mclust::adjustedRandIndex(castle_labels, immune.combined@ident[names(castle_labels)])
castle_VI <- mcclust::vi.dist(castle_labels, immune.combined@ident[names(castle_labels)])

# average expression per celltype
gem_avg <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- na.omit(names(castle_labels)[which(castle_labels == celltypes[i])])
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(stim.data.cpm[toupper(markers_filt$gene), cells, drop=FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(markers$gene[which(markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg) <- celltypes
colnames(gem_avg) <- celltypes

# Figure 2F right
pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 3, legend = F,
         border_color = F, show_rownames = T, scale = "row", width = 2, height = 2,
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/New plots/Kang_figure/castle_heatmap.pdf")

# ------------------------------------------------------------------------------------------------------------------------
# Use MNN
# Keep only cells retained from CCA too!
ctrl.cells <- intersect(colnames(ctrl.data.cpm), names(immune.combined@ident))
stim.cells <- intersect(colnames(stim.data.cpm), names(immune.combined@ident))

ctrl.data.filt <- log(ctrl.data.cpm[common_genes, ctrl.cells]+1)
stim.data.filt <- log(stim.data.cpm[common_genes, stim.cells]+1)

mnn_corrected <- mnnCorrect(ctrl.data.filt, stim.data.filt)

# Cluster mnn_corrected gem with Seurat
gem <- as.data.frame(mnn_corrected$corrected)
colnames(gem)[grep("X", colnames(gem))] <- colnames(stim.data.filt)

seurat <- CreateSeuratObject(counts = gem)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat,  do.print = FALSE)
seurat <- FindNeighbors(seurat)
seurat <- FindClusters(seurat)

ncol <- unique(Idents(seurat))
gem_avg <- matrix(NA, nrow = length(celltypes), ncol = length(ncol))
for (i in 1:length(ncol)) {
  cells <- WhichCells(seurat, idents = ncol[i])
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(gem[toupper(markers$gene), cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(markers$gene[which(markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg) <- celltypes
colnames(gem_avg) <- paste("Cl", ncol, sep = "_")

# Figure 2D right
pheatmap(gem_avg, border="white", color = colorspace::diverge_hsv(50), cluster_rows = F,
         cluster_cols = F, show_colnames = T, fontsize_row = 5, border_color = F,
         show_rownames = T, fontsize = 3, scale = "row", width = 2, height = 2, legend = F,
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/New plots/Kang_figure/MNN_heatmap.pdf")

#MNN_ARI <- mclust::adjustedRandIndex(Idents(seurat)[colnames(stim.data.filt)], immune.combined@ident[colnames(stim.data.filt)])
#MNN_VI <- mcclust::vi.dist(Idents(seurat)[colnames(stim.data.filt)], immune.combined@ident[colnames(stim.data.filt)])
MNN_ARI <- mclust::adjustedRandIndex(Idents(seurat), immune.combined@ident[colnames(gem)])
MNN_VI <- mcclust::vi.dist(Idents(seurat), immune.combined@ident[colnames(gem)])

# ------------------------------------------------------------------------------------------------------------------------
# Plot ARI and  VI
df <- data.frame(ARI = c(scID_ARI, scmap_ARI, castle_ARI, MNN_ARI),
                 VI = c(scID_VI, scmap_VI, castle_VI, MNN_VI),
                 method=c("scID", "scmap", "CaSTLe", "MNN"))
df.m <- reshape2::melt(df, identities="method")

ggplot(df.m, aes(y=value, x=method, fill=method)) + facet_wrap(~variable, scales = "free") + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() + scale_y_continuous(limits = c(0,0.8)) + 
  theme(legend.position="top", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")
