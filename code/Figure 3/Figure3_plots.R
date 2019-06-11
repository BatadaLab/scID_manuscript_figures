library(Seurat)
library(dplyr)
library(pheatmap)
library(scran)
library(scmap)
library(dplyr)
library(ggsci)
library(scID)
library(ggplot2)
library(ggpubr)
source("~/scID_manuscript_figures/code/run_other_methods.R")

# ------------------------------------------------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------------------------------------------------
E18_9k <- readRDS("~/scID_manuscript_figures/data/Figure3/reference_gem.rds")
E18_nuclei <- readRDS("~/scID_manuscript_figures/data/Figure3/target_gem.rds")

# CPM normalization
E18_9k_cpm <- counts_to_cpm(E18_9k)
E18_nuclei_cpm <- counts_to_cpm(E18_nuclei)

# Edit colnames of nuclei data to distinguish from E18_9k
colnames(E18_nuclei_cpm) <- paste(colnames(E18_nuclei_cpm), "nuc", sep = "_")
colnames(E18_nuclei) <- paste(colnames(E18_nuclei), "nuc", sep = "_")


# ------------------------------------------------------------------------------------------------------------------------
# Cluster E18_9k dataset with Seurat
# ------------------------------------------------------------------------------------------------------------------------
reference_clusters <- readRDS("~/scID_manuscript_figures/data/Figure3/reference_clusters.rds")
sobj_9k <- CreateSeuratObject(counts = E18_9k)
sobj_9k <- NormalizeData(sobj_9k)
sobj_9k <- FindVariableFeatures(sobj_9k)
sobj_9k <- ScaleData(sobj_9k)
sobj_9k <- RunPCA(sobj_9k,  do.print = FALSE)
Idents(sobj_9k) <- factor(reference_clusters)
sobj_9k <- RunTSNE(sobj_9k)

# Figure 3A left
TSNEPlot(sobj_9k, do.label = T, pt.size = 0.05, no.axes = T, no.legend = T)

# Find markers
markers <- FindAllMarkers(sobj_9k, only.pos = F, test.use = "MAST", logfc.threshold = 0.6, do.print = FALSE)
positive_markers <- markers[which(markers$avg_logFC > 0), ]
celltypes <- unique(Idents(sobj_9k))

gem_avg_9k <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- WhichCells(sobj_9k, idents = levels(celltypes)[i])
  if (length(celltypes) >= 1) {
    avg_exp <- rowMeans(E18_9k_cpm[positive_markers$gene, cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg_9k[j,i] <- mean(na.omit(avg_exp[positive_markers$gene[which(positive_markers$cluster == levels(celltypes)[j])]]))
  }
}
rownames(gem_avg_9k) <- paste("GS", levels(celltypes), sep = "_")
colnames(gem_avg_9k) <- paste("Cl", levels(celltypes), sep = "_")

# Figure 3B left
pheatmap(gem_avg_9k, border="white", color = colorspace::diverge_hsv(50), 
         cluster_rows = F, cluster_cols = F, show_colnames = F, fontsize = 5, 
         border_color = F, show_rownames = F, scale = "row")

gem_avg_norm <- apply(gem_avg_9k, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
Seurat_ref_counts <- table(colSums(gem_thr))/length(unique(Idents(sobj_9k)))
Seurat_ref_counts <- c(1-sum(Seurat_ref_counts), Seurat_ref_counts)
names(Seurat_ref_counts) <- c(0, names(Seurat_ref_counts[-1]))


# ------------------------------------------------------------------------------------------------------------------------
# Marker based identification
# ------------------------------------------------------------------------------------------------------------------------
markers_filt <- markers[-which(duplicated(markers$gene)), ]
top_markers <- as.data.frame(markers_filt %>% group_by(cluster) %>% top_n(2, avg_logFC))
top_markers <- top_markers[which(top_markers$gene %in% rownames(E18_nuclei_cpm)), ]
rownames(top_markers) <- top_markers$gene

# Normalize data by 99th percentile
nuclei_norm <- t(apply(E18_nuclei_cpm[intersect(top_markers$gene, rownames(E18_nuclei_cpm)), ], 1, function(x) normalize_gene(x)))

thresholds <- c(0.1, 0.25, 0.5, 0.75)
data <- data.frame(matrix(NA, length(thresholds), 3))
biomarker_labels <- list()
for (i in 1:length(thresholds)) {
  biomarker_labels[[i]] <- list()
  number_of_celltypes <- c()
  for (cell in colnames(nuclei_norm)) {
    ct <- c()
    for (gene in rownames(nuclei_norm)) {
      if (nuclei_norm[gene, cell] > thresholds[i]) {
        ct <- c(ct, top_markers$cluster[which(top_markers$gene == gene)])
      }
    }
    if (length(unique(ct )) == 0) {
      biomarker_labels[[i]] <- c(biomarker_labels[[i]], "orphan")
    } else if (length(unique(ct )) == 1) {
      biomarker_labels[[i]] <- c(biomarker_labels[[i]], unique(ct))
    } else {
      biomarker_labels[[i]] <- c(biomarker_labels[[i]], "ambiguous")
    }
    number_of_celltypes <- c(number_of_celltypes, length(unique(ct)))
  }
  data[i, 1] <- round(sum(number_of_celltypes == 0)*100/ncol(nuclei_norm), 2)
  data[i, 2] <- round(sum(number_of_celltypes == 1)*100/ncol(nuclei_norm), 2)
  data[i, 3] <- round(sum(number_of_celltypes > 1)*100/ncol(nuclei_norm), 2)
}

colnames(data) <- c("Orphan", "Identified", "Ambiguous")
data$threshold <- factor(thresholds, levels = thresholds)

# Figure 2C
df <- reshape2::melt(data, id.vars = c("threshold"))
ggplot(df, aes(fill=variable, y=value, x=threshold)) + 
  geom_bar(stat="identity") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits=c(0, 100)) + labs(title = "", x="", y="")


# ------------------------------------------------------------------------------------------------------------------------
# Cluster nuclei data with Seurat 
# ------------------------------------------------------------------------------------------------------------------------
nuclei_clustered <- CreateSeuratObject(counts = E18_nuclei)
nuclei_clustered <- NormalizeData(nuclei_clustered)
nuclei_clustered <- FindVariableFeatures(nuclei_clustered)
nuclei_clustered <- ScaleData(nuclei_clustered)
nuclei_clustered <- RunPCA(nuclei_clustered)
nuclei_clustered <- FindNeighbors(nuclei_clustered)
nuclei_clustered <- FindClusters(nuclei_clustered)
nuclei_clustered <- RunTSNE(nuclei_clustered)

# Figure 3A right
TSNEPlot(nuclei_clustered, pt.size = 0.05, label = T) + NoLegend() + NoAxes() 

markers_filt <- positive_markers[which(positive_markers$gene %in% rownames(E18_nuclei_cpm)), ]
nucClust <- unique(Idents(nuclei_clustered))
nuclei_avg <- matrix(NA, length(celltypes), length(nucClust))
for (i in 1:length(nucClust)) {
  cells <- WhichCells(nuclei_clustered, idents = levels(nucClust)[i])
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(E18_nuclei_cpm[markers_filt$gene, cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    nuclei_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == levels(celltypes)[j])]]))
  }
}
rownames(nuclei_avg) <- paste("GS", levels(celltypes), sep = "_")
colnames(nuclei_avg) <- paste("Cl", levels(nucClust), sep = "_")

# Figure 3B right
pheatmap(nuclei_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = F, fontsize = 5, border_color = F,
         show_rownames = F, scale = "row")

gem_avg_norm <- apply(nuclei_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
Seurat_target_counts <- table(colSums(gem_thr))/length(unique(Idents(nuclei_clustered)))
Seurat_target_counts <- c(1-sum(Seurat_target_counts), Seurat_target_counts)
names(Seurat_target_counts) <- c(0, names(Seurat_target_counts[-1]))


# ------------------------------------------------------------------------------------------------------------------------
# Find cells with Seurat batch effect correction 
# ------------------------------------------------------------------------------------------------------------------------
reference <- CreateSeuratObject(counts = E18_9k, project = "Reference", min.cells = 5)
reference$stim <- "Reference"
reference <- NormalizeData(reference, verbose = FALSE)
reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 2000)

target <- CreateSeuratObject(counts = E18_nuclei, project = "Target", min.cells = 5)
target$stim <- "Target"
target <- NormalizeData(target, verbose = FALSE)
target <- FindVariableFeatures(target, selection.method = "vst", nfeatures = 2000)

combined.anchors <- FindIntegrationAnchors(object.list = list(reference, target), dims = 1:20)
cells.combined <- IntegrateData(anchorset = combined.anchors, dims = 1:20)

DefaultAssay(cells.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
cells.combined <- ScaleData(cells.combined, verbose = FALSE)
cells.combined <- RunPCA(cells.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
cells.combined <- FindNeighbors(cells.combined, reduction = "pca", dims = 1:20)
cells.combined <- FindClusters(cells.combined, resolution = 0.5)

markers_filt <- positive_markers[which(positive_markers$gene %in% rownames(E18_nuclei_cpm)), ]

combClust <- unique(Idents(cells.combined)) 
gem_avg_ds <- matrix(NA, nrow = length(celltypes), ncol = length(combClust))
for (i in 1:length(combClust)) {
  cells <- intersect(WhichCells(cells.combined, idents = levels(combClust)[i]), colnames(E18_nuclei_cpm))
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(E18_nuclei_cpm[toupper(markers_filt$gene), cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg_ds[j,i] <- mean(na.omit(avg_exp[toupper(markers_filt$gene[which(markers_filt$cluster == levels(celltypes)[j])])]))
  }
}
rownames(gem_avg_ds) <- paste("GS", levels(celltypes), sep = "_")
colnames(gem_avg_ds) <- paste("Cl", levels(combClust), sep = "_")

gem_avg_ds <- gem_avg_ds[, complete.cases(t(gem_avg_ds))]
gem_avg_norm <- apply(gem_avg_ds, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
CCA_ref_counts <- table(colSums(gem_thr))/length(unique(Idents(cells.combined)))
CCA_ref_counts <- c(1-sum(CCA_ref_counts), CCA_ref_counts)
names(CCA_ref_counts) <- c(0, names(CCA_ref_counts[-1]))

CCA_labels <- Idents(cells.combined)[colnames(E18_nuclei)]
CCA_ARI <- c()
CCA_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(biomarker_labels[[i]])
  names(marker_labels) <- colnames(E18_nuclei)
  CCA_ARI[i] <- mclust::adjustedRandIndex(marker_labels, CCA_labels)
  CCA_VI[i] <- mcclust::vi.dist(marker_labels, CCA_labels)
}


# ------------------------------------------------------------------------------------------------------------------------
# Find cells with MNN batch effect correction
# ------------------------------------------------------------------------------------------------------------------------
mnn_labels <- runMNN(reference_gem = log(E18_9k_cpm+1), target_gem = log(E18_nuclei_cpm+1))

mnn_labels_reference <- mnn_labels$target_labels
ncol <- unique(mnn_labels_reference)
gem_avg <- matrix(NA, nrow = length(celltypes), ncol = length(levels(ncol)))
for (i in 1:length(levels(ncol))) {
  cells <- names(mnn_labels_reference)[which(mnn_labels_reference == levels(ncol)[i])]
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(E18_nuclei_cpm[markers_filt$gene, cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == levels(celltypes)[j])]]))
  }
}
rownames(gem_avg) <- paste("GS", levels(celltypes), sep = "_")
colnames(gem_avg) <- paste("Cl", levels(ncol), sep = "_")

gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
MNN_ref_counts <- table(colSums(gem_thr))/length(unique(Idents(cells.combined)))
MNN_ref_counts <- c(1-sum(MNN_ref_counts), MNN_ref_counts)
names(MNN_ref_counts) <- c(0, names(MNN_ref_counts[-1]))

MNN_ARI <- c()
MNN_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(biomarker_labels[[i]])
  names(marker_labels) <- colnames(E18_nuclei)
  MNN_ARI[i] <- mclust::adjustedRandIndex(marker_labels, mnn_labels_reference)
  MNN_VI[i] <- mcclust::vi.dist(marker_labels, mnn_labels_reference)
}


# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scID 
# ------------------------------------------------------------------------------------------------------------------------
scID_output <- scid_multiclass(target_gem = E18_nuclei_cpm, markers = markers, estimate_weights_from_target = T)
scID_labels <- scID_output$labels

# Get ARI and VI (variation of information) with marker based for different thresholds
scID_ARI <- c()
scID_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(biomarker_labels[[i]])
  names(marker_labels) <- colnames(E18_nuclei_cpm)
  cells <- names(marker_labels)[which(!marker_labels %in% c("ambiguous", "orphan"))]
  scID_ARI[i] <- mclust::adjustedRandIndex(marker_labels[cells], scID_labels[cells])
  scID_VI[i] <- mcclust::vi.dist(marker_labels[cells], scID_labels[cells])
}

gem_avg <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- names(scID_labels)[which(scID_labels == levels(celltypes)[i])]
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(E18_nuclei_cpm[markers_filt$gene, cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == levels(celltypes)[j])]]))
  }
}
rownames(gem_avg) <- paste("GS", levels(celltypes), sep = "_")
colnames(gem_avg) <- paste("Cl", levels(celltypes), sep = "_")

gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
scID_counts <- table(colSums(gem_thr))/length(unique(celltypes))
scID_counts <- c(1-sum(scID_counts), scID_counts)
names(scID_counts) <- c(0, names(scID_counts[-1]))


# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scmap
# ------------------------------------------------------------------------------------------------------------------------
scmap_labels <- run_scmap(reference_gem = log(E18_9k_cpm + 1), reference_labels = Idents(sobj_9k), 
                          target_gem = log(E18_nuclei_cpm + 1))
table(scmap_labels)

# Get ARI with marker based for different thresholds
scmap_ARI <- c()
scmap_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(biomarker_labels[[i]])
  names(marker_labels) <- colnames(E18_nuclei_cpm)
  cells <- names(which(marker_labels != "unclassified"))
  scmap_ARI[i] <- mclust::adjustedRandIndex(marker_labels[cells], scmap_labels[cells])
  scmap_VI[i] <- mcclust::vi.dist(marker_labels[cells], scmap_labels[cells])
}

gem_avg <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- names(scmap_labels)[which(scmap_labels == levels(celltypes)[i])]
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(E18_nuclei_cpm[markers_filt$gene, cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == levels(celltypes)[j])]]))
  }
}
rownames(gem_avg) <- paste("GS", levels(celltypes), sep = "_")
colnames(gem_avg) <- paste("Cl", levels(celltypes), sep = "_")

gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
scmap_counts <- table(colSums(gem_thr))/length(unique(celltypes))
scmap_counts <- c(1-sum(scmap_counts), scmap_counts)
names(scmap_counts) <- c(0, names(scmap_counts[-1]))


# ------------------------------------------------------------------------------------------------------------------------
# find equivalent cells with CaSTLe
# ------------------------------------------------------------------------------------------------------------------------
castle_labels <- runCastle(target_gem = log(E18_nuclei_cpm+1), source_gem = log(E18_9k_cpm+1), source_identities = Idents(sobj_9k))

# Get ARI and VI with marker based for different thresholds
castle_ARI <- c()
castle_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(biomarker_labels[[i]])
  names(marker_labels) <- colnames(E18_nuclei_cpm)
  castle_ARI[i] <- mclust::adjustedRandIndex(marker_labels, castle_labels)
  castle_VI[i] <- mcclust::vi.dist(marker_labels, castle_labels)
}

# average expression per celltype
gem_avg <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- names(castle_labels)[which(castle_labels == levels(celltypes)[i])]
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(E18_nuclei_cpm[markers_filt$gene, cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == levels(celltypes)[j])]]))
  }
}
rownames(gem_avg) <- paste("GS", levels(celltypes), sep = "_")
colnames(gem_avg) <- paste("Cl", levels(celltypes), sep = "_")

# Remove NA columns
gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
castle_counts <- table(colSums(gem_thr))/length(unique(celltypes))
castle_counts <- c(1-sum(castle_counts), castle_counts)
names(castle_counts) <- c(0, names(castle_counts[-1]))


# ------------------------------------------------------------------------------------------------------------------------
# plot ARI
# ------------------------------------------------------------------------------------------------------------------------
df <- data.frame(value=c(scID_ARI, scmap_ARI, castle_ARI, MNN_ARI, CCA_ARI), 
                 threshold=rep(thresholds, 5), 
                 method=c(rep("scID", length(thresholds)), rep("scmap", length(thresholds)), rep("CaSTLe", length(thresholds)), 
                          rep("MNN", length(thresholds)), rep("CCA", length(thresholds))))
df$threshold <- factor(df$threshold, levels = thresholds)

ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "", x="", y="")


# ------------------------------------------------------------------------------------------------------------------------
# plot VI 
# ------------------------------------------------------------------------------------------------------------------------
df <- data.frame(value=c(scID_VI, scmap_VI, castle_VI, MNN_VI, CCA_VI), 
                 threshold=rep(thresholds, 5), 
                 method=c(rep("scID", length(thresholds)), rep("scmap", length(thresholds)), rep("CaSTLe", length(thresholds)), 
                          rep("MNN", length(thresholds)), rep("CCA", length(thresholds))))
df$threshold <- factor(df$threshold, levels = thresholds)

ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "", x="", y="")
  

# ------------------------------------------------------------------------------------------------------------------------
# plot number of signatures per cluster 
# ------------------------------------------------------------------------------------------------------------------------
df <- data.frame(counts = c(Seurat_target_counts, scID_counts, scmap_counts, castle_counts, Seurat_ref_counts, MNN_ref_counts, CCA_ref_counts),
                 numbers = c(names(Seurat_target_counts), names(scID_counts), names(scmap_counts), names(castle_counts), 
                             names(Seurat_ref_counts), names(MNN_ref_counts), names(CCA_ref_counts)), 
                 method = c(rep("Seurat", length(Seurat_target_counts)), 
                            rep("scID", length(scID_counts)), 
                            rep("scmap", length(scmap_counts)), 
                            rep("CaSTLe", length(castle_counts)), 
                            rep("reference", length(Seurat_ref_counts)), 
                            rep("MNN", length(MNN_ref_counts)), 
                            rep("CCA", length(CCA_ref_counts))))

ggplot(df, aes(y=counts, x=method, fill=numbers)) + 
  geom_bar(stat="identity") + scale_fill_grey(start = 1, end = 0) +
  theme(legend.position="top", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "", x="", y="")
