library(Seurat)
library(dplyr)
library(pheatmap)
library(scran)
library(scmap)
library(ggsci)
library(scID)
library(cowplot)

# Read data
ctrl.data <- read.table(file = "~/scID_manuscript_figures/data/Kang_data/immune_control_expression_matrix.txt.gz", sep = "\t")
stim.data <- read.table(file = "~/scID_manuscript_figures/data/Kang_data/immune_stimulated_expression_matrix.txt.gz", sep = "\t")

# ------------------------------------------------------------------------------------------------------------------------
# Cluster with Seurat CCA and use labels as reference
# Set up control object
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
ctrl$stim <- "CTRL"
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim$stim <- "STIM"
stim <- subset(stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(immune.combined, reduction = "umap", split.by = "stim")

FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", 
                                          "CCL2", "PPBP"), min.cutoff = "q9")


immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
                                `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", 
                                `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")

DimPlot(immune.combined, label = TRUE)






# ------------------------------------------------------------------------------------------------------------------------
# Find markers of each cluster from DropSeq or read pre-computed markers
so_ds <- CreateSeuratObject(raw.data = Ds_gem)
so_ds <- NormalizeData(so_ds)
so_ds <- ScaleData(so_ds)
so_ds@ident <- factor(identities, levels = celltypes)
# markers <- FindAllMarkers(so_ds, only.pos = TRUE, test.use = "MAST", logfc.threshold = 0.7)

markers <- read.delim("~/scID_manuscript_figures/data/Figure2/markers.txt", stringsAsFactors = F)

gem_avg_ds <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- WhichCells(so_ds, celltypes[i])
  if (length(cells) > 1) {
    avg_exp <- rowMeans(Ds_gem[toupper(markers$gene), cells])
  } else if (length(cells) == 1) {
    avg_exp <- gem[toupper(markers$gene), cells]
    names(avg_exp) <- toupper(markers$gene)
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg_ds[j,i] <- mean(na.omit(avg_exp[toupper(markers$gene[which(markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg_ds) <- celltypes
colnames(gem_avg_ds) <- celltypes

# Figure 2B left
pheatmap(gem_avg_ds, border="white", color = colorspace::diverge_hsv(50), cluster_rows = F,
         cluster_cols = F, show_colnames = T, fontsize_row = 5, border_color = F,show_rownames = T,
         scale = "row", fontsize_col = 5)

gem_avg_norm <- apply(gem_avg_ds, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
Seurat_ref_counts <- table(colSums(gem_thr))/length(unique(so_ds@ident))
Seurat_ref_counts <- c(1-sum(Seurat_ref_counts), Seurat_ref_counts)
names(Seurat_ref_counts) <- c(0, names(Seurat_ref_counts[-1]))

# ------------------------------------------------------------------------------------------------------------------------
# Cluster SS2 data with Seurat
SS2_clustered <- CreateSeuratObject(raw.data = log2(ss2_gem+1))
SS2_clustered <- NormalizeData(SS2_clustered)
SS2_clustered <- ScaleData(SS2_clustered)
SS2_clustered <- FindVariableGenes(SS2_clustered, do.contour=F)
SS2_clustered <- RunPCA(SS2_clustered,  do.print = T)
SS2_clustered <- ProjectPCA(SS2_clustered)
PCElbowPlot(SS2_clustered)
dims <- 3
SS2_clustered <- FindClusters(SS2_clustered, dims.use = 1:dims)
SS2_clustered <- RunTSNE(SS2_clustered, dims.use = 1:dims, do.fast = F)

ncol <- length(unique(SS2_clustered@ident))

gem_avg <- matrix(NA, length(celltypes), ncol)
for (i in 1:ncol) {
  cells <- WhichCells(SS2_clustered, i-1)
  if (length(cells) > 1) {
    avg_exp <- rowMeans(ss2_gem[toupper(markers$gene), cells])
  } else if (length(cells) == 1) {
    avg_exp <- ss2_gem[toupper(markers$gene), cells]
    names(avg_exp) <- toupper(markers$gene)
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(markers$gene[which(markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg) <- celltypes
colnames(gem_avg) <- paste("Cl", unique(SS2_clustered@ident), sep = "_")

# Figure 2B right
pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, border_color = F,
         show_rownames = T, scale = "row")

gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
Seurat_target_counts <- table(colSums(gem_thr))/length(unique(SS2_clustered@ident))
Seurat_target_counts <- c(1-sum(Seurat_target_counts), Seurat_target_counts)
names(Seurat_target_counts) <- c(0, names(Seurat_target_counts[-1]))

rm(SS2_clustered)

# -----------------------------------------------------------------------------------------------------------------------
# Marker based identification
# remove duplicated markers
markers_filt <- markers[-which(duplicated(markers$gene)), ]
top_markers <- as.data.frame(markers_filt %>% group_by(cluster) %>% top_n(1, avg_logFC))
top_markers <- top_markers[which(top_markers$gene %in% rownames(ss2_gem)), ]
rownames(top_markers) <- top_markers$gene

# Normalize data by 99th percentile
ss2_gem_norm <- t(apply(ss2_gem[top_markers$gene, ], 1, function(x) scID::normalize_gem(x)))

thresholds <- c(0, 0.25, 0.5, 0.75)
data <- data.frame(matrix(NA, length(thresholds), 3))
biomarker_labels <- list()
for (i in 1:length(thresholds)) {
  biomarker_labels[[i]] <- list()
  for (cell in colnames(ss2_gem_norm)) {
    #expressed_genes <- which(ss2_gem_norm[, cell] > thresholds[i])
    expressed_genes <- c()
    for (gene in rownames(ss2_gem_norm)) {
      if (ss2_gem_norm[gene, cell] > quantile(unlist(ss2_gem_norm[gene, ]), thresholds[i])) {
        expressed_genes <- c(expressed_genes, gene)
      }
    }
    ct <- unlist(top_markers[expressed_genes, "cluster"])
    if (length(unique(ct)) == 1) {
      biomarker_labels[[i]] <- c(biomarker_labels[[i]], unique(ct))
    } else if (length(unique(ct)) == 0) {
      biomarker_labels[[i]] <- c(biomarker_labels[[i]], "orphan")
    } else {
      biomarker_labels[[i]] <- c(biomarker_labels[[i]], "ambiguous")
    }
  }
  biomarker_labels[[i]] <- unlist(biomarker_labels[[i]])
  data[i, 1] <- round(sum(biomarker_labels[[i]] == "orphan")*100/ncol(ss2_gem_norm), 2)
  data[i, 3] <- round(sum(biomarker_labels[[i]] == "ambiguous")*100/ncol(ss2_gem_norm), 2)
  data[i, 2] <- 100 - (data[i,1] + data[i,3])
}
colnames(data) <- c("Orphan", "Identified", "Ambiguous")
data$threshold <- thresholds

# Figure 2C
df <- reshape2::melt(data, id.vars = c("threshold"))
ggplot(df, aes(fill=variable, y=value, x=threshold)) + 
  geom_bar(stat="identity") + scale_fill_jco() + 
  theme(legend.position="top", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_y_continuous(limits=c(0, 100)) + labs(title = "", x="", y="")

# ------------------------------------------------------------------------------------------------------------------------
# Find cells with Seurat batch effect correction
# Set up DropSeq object
dropseq <- CreateSeuratObject(raw.data = Ds_gem, project = "Dropseq", min.cells = 5)
dropseq@meta.data$stim <- "Dropseq"
dropseq <- FilterCells(dropseq, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
dropseq <- NormalizeData(dropseq)
dropseq <- ScaleData(dropseq, display.progress = F)
# Set up SS2 object
ss2 <- CreateSeuratObject(raw.data = log2(ss2_gem+1), project = "SS2", min.cells = 5)
ss2@meta.data$stim <- "SS2"
ss2 <- FilterCells(ss2, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
ss2 <- NormalizeData(ss2)
ss2 <- ScaleData(ss2, display.progress = F)
# Gene selection for input to CCA
dropseq <- FindVariableGenes(dropseq, do.plot = F)
ss2 <- FindVariableGenes(ss2, do.plot = F)
g.1 <- head(rownames(dropseq@hvg.info), 1000)
g.2 <- head(rownames(ss2@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(dropseq@scale.data))
genes.use <- intersect(genes.use, rownames(ss2@scale.data))
# Run CCA to combine datasets
cells.combined <- RunCCA(dropseq, ss2, genes.use = genes.use, num.cc = 20)
rm(dropseq)
rm(ss2)
p3 <- MetageneBicorPlot(cells.combined, grouping.var = "stim", dims.eval = 1:20, display.progress = FALSE)
# Align the CCA subspaces
cells.combined <- AlignSubspace(cells.combined, reduction.type = "cca", grouping.var = "stim", dims.align = 1:20)
# Cluster combined dataset
#cells.combined <- RunTSNE(cells.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
cells.combined <- FindClusters(cells.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)

ncol <- length(unique(cells.combined@ident)) - 1
gem_avg_ds <- matrix(NA, nrow = length(celltypes), ncol = ncol)
gem <- cells.combined@scale.data
for (i in 1:ncol) {
  cells <- intersect(WhichCells(cells.combined, i), colnames(Ds_gem))
  if (length(cells) > 1) {
    avg_exp <- rowMeans(gem[toupper(markers$gene), cells])
  } else if (length(cells) == 1) {
    avg_exp <- gem[toupper(markers$gene), cells]
    names(avg_exp) <- toupper(markers$gene)
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg_ds[j,i] <- mean(na.omit(avg_exp[toupper(markers$gene[which(markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg_ds) <- celltypes
colnames(gem_avg_ds) <- paste("Cl", unique(cells.combined@ident), sep = "_")

# Figure 2D left
pheatmap(gem_avg_ds, border="white", color = colorspace::diverge_hsv(50), cluster_rows = F,
         cluster_cols = F, show_colnames = T, fontsize_row = 5, border_color = F,
         show_rownames = T, fontsize_col = 5, scale = "row")

gem_avg_norm <- apply(gem_avg_ds, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
CCA_ref_counts <- table(colSums(gem_thr))/length(unique(cells.combined@ident))
CCA_ref_counts <- c(1-sum(CCA_ref_counts), CCA_ref_counts)
names(CCA_ref_counts) <- c(0, names(CCA_ref_counts[-1]))

rm(cells.combined)

# ------------------------------------------------------------------------------------------------------------------------
# Find cells with MNN batch effect correction
common_genes <- intersect(rownames(Ds_gem), rownames(ss2_gem))
mnn_corrected <- mnnCorrect(as.matrix(Ds_gem[common_genes, ]), as.matrix(log2(ss2_gem[common_genes, ]+1)))

# Cluster mnn_corrected gem with Seurat
gem <- as.data.frame(mnn_corrected$corrected)
colnames(gem)[grep("X", colnames(gem))] <- colnames(ss2_gem)

seurat <- CreateSeuratObject(raw.data = gem)
seurat <- NormalizeData(seurat)
seurat <- ScaleData(seurat)
seurat <- FindVariableGenes(seurat, x.low.cutoff = 0.05, x.high.cutoff = 1, y.cutoff = 0.5, do.contour = F)
seurat <- RunPCA(seurat,  do.print = FALSE)
seurat <- ProjectPCA(seurat)
seurat <- FindClusters(seurat, dims.use = 1:5, save.SNN = F)

ncol <- unique(seurat@ident)
gem_avg <- matrix(NA, nrow = length(celltypes), ncol = length(ncol))
for (i in 1:length(ncol)) {
  cells <- intersect(WhichCells(seurat, ncol[i]), colnames(Ds_gem))
  if (length(cells) > 1) {
    avg_exp <- rowMeans(gem[toupper(markers$gene), cells])
  } else if (length(cells) == 1) {
    avg_exp <- gem[toupper(markers$gene), cells]
    names(avg_exp) <- toupper(markers$gene)
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
         cluster_cols = F, show_colnames = F, fontsize_row = 5, border_color = F,
         show_rownames = T, fontsize = 5, scale = "row")

gem_avg_norm <- apply(gem_avg[, -7], 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
MNN_ref_counts <- table(colSums(gem_thr))/length(unique(seurat@ident))
MNN_ref_counts <- c(1-sum(MNN_ref_counts), MNN_ref_counts)
names(MNN_ref_counts) <- c(0, names(MNN_ref_counts[-1]))

rm(seurat)

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scID
# remove bulk cells
scID_output <- scid_multiclass(target_gem = ss2_gem, markers = markers, likelihood_threshold = 0.99)
scID_labels <- scID_output$labels

# Get ARI and VI (variation of information) with marker based for different thresholds
scID_ARI <- c()
scID_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(biomarker_labels[[i]])
  names(marker_labels) <- colnames(ss2_gem)
  cells <- names(marker_labels)[which(!marker_labels %in% c("ambiguous", "orphan"))]
  scID_ARI[i] <- mclust::adjustedRandIndex(marker_labels[cells], scID_labels[cells])
  scID_VI[i] <- mcclust::vi.dist(marker_labels[cells], scID_labels[cells])
}

# average expression per celltype
gem_avg <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- na.omit(names(scID_labels)[which(scID_labels == celltypes[i])])
  if (length(cells) > 1) {
    avg_exp <- rowMeans(ss2_gem[toupper(markers$gene), cells])
  } else if (length(cells) == 1) {
    avg_exp <- ss2_gem[toupper(markers$gene), cells]
    names(avg_exp) <- toupper(markers$gene)
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
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, 
         border_color = F, show_rownames = T, scale = "row")

# Remove NA columns
NA_cols <- which(apply(gem_avg, 2, function(x) all(is.na(x))))
gem_avg_norm <- apply(gem_avg[, -NA_cols], 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
scID_counts <- table(colSums(gem_thr))/length(unique(celltypes))
scID_counts <- c(1-sum(scID_counts), scID_counts)
names(scID_counts) <- c(0, names(scID_counts[-1]))

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scmap
common_genes <- intersect(rownames(Ds_gem), rownames(ss2_gem))

ann <- data.frame(cell_type1=identities)

# Create index
ref_sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(Ds_gem[common_genes, ])), colData = ann)
# use gene names as feature symbols
rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
# remove features with duplicated names
ref_sce <- ref_sce[!duplicated(rownames(ref_sce)), ]
ref_sce <- selectFeatures(ref_sce, suppress_plot = FALSE, n_features = 150)
ref_sce <- indexCluster(ref_sce)

# Create sce for testing data (SS2)
sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(log2(ss2_gem[common_genes, ]+1))))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]
sce <- selectFeatures(sce, suppress_plot = FALSE, n_features = 150)

scmapCluster_results <- scmapCluster(
  projection = sce, 
  index_list = list(
    DropSeq = metadata(ref_sce)$scmap_cluster_index
  )
)

scmap_labels <- scmapCluster_results$combined_labs
names(scmap_labels) <- colnames(ss2_gem)

# Get ARI with marker based for different thresholds
scmap_ARI <- c()
scmap_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(marker_based_labels[[i]])
  names(marker_labels) <- colnames(ss2_gem)
  cells <- names(marker_labels)[which(!marker_labels %in% c("ambiguous", "orphan"))]
  scmap_ARI[i] <- mclust::adjustedRandIndex(marker_labels[cells], scmap_labels[cells])
  scmap_VI[i] <- mcclust::vi.dist(marker_labels[cells], scmap_labels[cells])
}

gem_avg <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- names(scmap_labels)[which(scmap_labels == celltypes[i])]
  if (length(cells) > 1) {
    avg_exp <- rowMeans(ss2_gem[toupper(markers$gene), cells])
  } else if (length(cells) == 1) {
    avg_exp <- ss2_gem[toupper(markers$gene), cells]
    names(avg_exp) <- toupper(markers$gene)
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(markers$gene[which(markers$cluster == celltypes[j])])]))
  }
}
colnames(gem_avg) <- celltypes
rownames(gem_avg) <- celltypes

# Figure 2F left
pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, 
         border_color = F,show_rownames = T, scale = "row")

NA_cols <- which(apply(gem_avg, 2, function(x) all(is.na(x))))
gem_avg_norm <- apply(gem_avg[ , -NA_cols], 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
scmap_counts <- table(colSums(gem_thr))/length(unique(celltypes))
scmap_counts <- c(1-sum(scmap_counts), scmap_counts)
names(scmap_counts) <- c(0, names(scmap_counts[-1]))

# ---------------------------------------------------------------
# plot ARI - Figure 2 H
df <- data.frame(value=c(scID_ARI, scmap_ARI), threshold=c(thresholds, thresholds), method=c(rep("scID", 4), rep("scmap", 4)))
df$threshold <- factor(df$threshold, levels = thresholds)

ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="none", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")

# ---------------------------------------------------------------
# plot VI - Figure 2 I
df <- data.frame(value=c(scID_VI, scmap_VI), threshold=c(thresholds, thresholds), method=c(rep("scID", 4), rep("scmap", 4)))
df$threshold <- factor(df$threshold, levels = thresholds)

ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="none", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")

# ---------------------------------------------------------------
# plot number of signatures per cluster (reference data) - Figure 2 E
df <- data.frame(counts = c(Seurat_ref_counts, MNN_ref_counts, CCA_ref_counts),
                 numbers = c(names(Seurat_ref_counts), names(MNN_ref_counts), names(CCA_ref_counts)), 
                 method = c(rep("Seurat", length(Seurat_ref_counts)), 
                            rep("MNN", length(MNN_ref_counts)), 
                            rep("CCA", length(CCA_ref_counts))))

ggplot(df, aes(y=counts, x=method, fill=numbers)) + 
  geom_bar(stat="identity") + scale_fill_grey(start = 1, end = 0) +
  theme(legend.position="top", text = element_text(size=5), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")

# ---------------------------------------------------------------
# plot number of signatures per cluster (target data) - Figure 2 G
df <- data.frame(counts = c(Seurat_target_counts, scID_counts, scmap_counts, Seurat_ref_counts),
                 numbers = c(names(Seurat_target_counts), names(scID_counts), names(scmap_counts), names(Seurat_ref_counts)), 
                 method = c(rep("Seurat", length(Seurat_target_counts)), 
                            rep("scID", length(scID_counts)), 
                            rep("scmap", length(scmap_counts)), 
                            rep("reference", length(Seurat_ref_counts))))

df$numbers <- factor(df$numbers, levels = c(0,1,2,3,4,5,11))


ggplot(df, aes(y=counts, x=method, fill=numbers)) + 
  geom_bar(stat="identity") + scale_fill_grey(start = 1, end = 0) +
  theme(legend.position="top", text = element_text(size=5), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")
