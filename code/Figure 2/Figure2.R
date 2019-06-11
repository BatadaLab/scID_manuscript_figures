library(Seurat)
library(dplyr)
library(pheatmap)
library(scran)
library(scmap)
library(ggsci)
library(scID)
library(ggplot2)
source("~/scID_manuscript_figures/code/run_other_methods.R")

# ------------------------------------------------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------------------------------------------------
Ds_gem <- readRDS("~/scID_manuscript_figures/data/Figure2/Reference_gem.rds")
rownames(Ds_gem) <- toupper(rownames(Ds_gem))
identities <- readRDS("~/scID_manuscript_figures/data/Figure2/Reference_clusters.rds")

ss2_gem <- readRDS("~/scID_manuscript_figures/data/Figure2/Target_gem.rds")

# ------------------------------------------------------------------------------------------------------------------------
# Find markers of each cluster from DropSeq (Fig 4b left)
# ------------------------------------------------------------------------------------------------------------------------
so_ds <- CreateSeuratObject(counts = Ds_gem)
so_ds <- NormalizeData(so_ds)
so_ds <- ScaleData(so_ds)
Idents(so_ds) <- factor(identities, levels = c("RBC (Rod Bipolar cell)", "MG (Mueller Glia)", "BC5A (Cone Bipolar cell 5A)", 
                                             "BC7 (Cone Bipolar cell 7)", "BC6", "BC5C", "BC1A", "BC3B", "BC1B", 
                                             "BC2", "BC5D", "BC3A", "BC5B", "BC4", "BC8/9 (mixture of BC8 and BC9)", 
                                             "AC (Amacrine cell)", "Rod Photoreceptors", "Cone Photoreceptors"))

markers <- FindAllMarkers(so_ds, only.pos = FALSE, test.use = "MAST", logfc.threshold = 0.7)
# Keep balanced number of markers per cluster
positive_markers <- markers[which(markers$avg_logFC > 0), ]
top100_positive <- positive_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
negative_markers <- markers[which(markers$avg_logFC < 0), ]
top100_negative <- negative_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)

markers <- rbind(top100_negative, top100_positive)
write.table("~/scID_manuscript_figures/data/Figure2/markers_10June.txt", quote = F, sep = "\t")
positive_markers <- markers[which(markers$avg_logFC > 0), ]

celltypes <- c("RBC (Rod Bipolar cell)", "MG (Mueller Glia)", "BC5A (Cone Bipolar cell 5A)", 
               "BC7 (Cone Bipolar cell 7)", "BC6", "BC5C", "BC1A", "BC3B", "BC1B", 
               "BC2", "BC5D", "BC3A", "BC5B", "BC4", "BC8/9 (mixture of BC8 and BC9)", 
               "AC (Amacrine cell)", "Rod Photoreceptors", "Cone Photoreceptors")

gem_avg_ds <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- WhichCells(so_ds, idents = celltypes[i])
  if (length(cells) > 1) {
    avg_exp <- rowMeans(Ds_gem[toupper(positive_markers$gene), cells])
  } else if (length(cells) == 1) {
    avg_exp <- gem[toupper(positive_markers$gene), cells]
    names(avg_exp) <- toupper(positive_markers$gene)
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg_ds[j,i] <- mean(na.omit(avg_exp[toupper(positive_markers$gene[which(positive_markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg_ds) <- celltypes
colnames(gem_avg_ds) <- celltypes
pheatmap(gem_avg_ds, border="white", color = colorspace::diverge_hsv(50), cluster_rows = F,
         cluster_cols = F, show_colnames = F, fontsize_row = 3, border_color = F,show_rownames = F, scale = "row")

gem_avg_norm <- apply(gem_avg_ds, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
Seurat_ref_counts <- table(colSums(gem_thr))/length(unique(Idents(so_ds)))
Seurat_ref_counts <- c(1-sum(Seurat_ref_counts), Seurat_ref_counts)
names(Seurat_ref_counts) <- c(0, names(Seurat_ref_counts[-1]))


# ------------------------------------------------------------------------------------------------------------------------
# Cluster SS2 data with Seurat
# ------------------------------------------------------------------------------------------------------------------------
SS2_clustered <- CreateSeuratObject(counts = log2(ss2_gem+1))
SS2_clustered <- FindVariableFeatures(SS2_clustered, do.contour=F)
SS2_clustered <- ScaleData(SS2_clustered)
SS2_clustered <- RunPCA(SS2_clustered,  do.print = T)
SS2_clustered <- FindNeighbors(SS2_clustered)
SS2_clustered <- FindClusters(SS2_clustered)
SS2_clustered <- RunTSNE(SS2_clustered)

TSNEPlot(SS2_clustered, pt.size = 0.05, label = T) + NoLegend() + NoAxes() 

ncol <- unique(Idents(SS2_clustered))
gem_avg <- matrix(NA, length(celltypes), length(ncol))
for (i in 1:length(ncol)) {
  cells <- WhichCells(SS2_clustered, idents = levels(ncol)[i])
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(ss2_gem[toupper(positive_markers$gene), cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(positive_markers$gene[which(positive_markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg) <- celltypes
colnames(gem_avg) <- paste("Cl", levels(ncol), sep = "_")

# Figure 2B right
pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = F, fontsize = 5, border_color = F,
         show_rownames = F, scale = "row")

gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
Seurat_target_counts <- table(colSums(gem_thr))/length(unique(Idents(SS2_clustered)))
Seurat_target_counts <- c(1-sum(Seurat_target_counts), Seurat_target_counts)
names(Seurat_target_counts) <- c(0, names(Seurat_target_counts[-1]))


# -----------------------------------------------------------------------------------------------------------------------
# Marker based identification
# ------------------------------------------------------------------------------------------------------------------------
# Read signature genes
top_markers <- top100_positive %>% group_by(cluster) %>% top_n(1, avg_logFC)

# Normalize data by 99th percentile
ss2_gem_norm <- t(apply(ss2_gem[top_markers$gene, ], 1, function(x) scID::normalize_gene(x)))

thresholds <- c(0.1, 0.25, 0.5, 0.75)
data <- data.frame(matrix(NA, length(thresholds), 3))
marker_based_labels <- list()
for (i in 1:length(thresholds)) {
  marker_based_labels[[i]] <- list()
  number_of_celltypes <- c()
  for (cell in colnames(ss2_gem_norm)) {
    ct <- c()
    for (gene in rownames(ss2_gem_norm)) {
      if (ss2_gem_norm[gene, cell] > thresholds[i]) {
        ct <- c(ct, top_markers$cluster[which(top_markers$gene == gene)])
      }
    }
    if (length(unique(ct )) == 0) {
      marker_based_labels[[i]] <- c(marker_based_labels[[i]], "orphan")
    } else if (length(unique(ct )) == 1) {
      marker_based_labels[[i]] <- c(marker_based_labels[[i]], unique(ct))
    } else {
      marker_based_labels[[i]] <- c(marker_based_labels[[i]], "ambiguous")
    }
    number_of_celltypes <- c(number_of_celltypes, length(unique(ct)))
  }
  data[i, 1] <- round(sum(number_of_celltypes == 0)*100/ncol(ss2_gem_norm), 2)
  data[i, 2] <- round(sum(number_of_celltypes == 1)*100/ncol(ss2_gem_norm), 2)
  data[i, 3] <- round(sum(number_of_celltypes > 1)*100/ncol(ss2_gem_norm), 2)
}

colnames(data) <- c("Orphan", "Classified", "Ambiguous")
data$threshold <- factor(thresholds, levels = thresholds)

df <- reshape2::melt(data, id.vars = c("threshold"))
ggplot(df, aes(fill=variable, y=value, x=threshold)) + 
  geom_bar(stat="identity") + scale_fill_jco() + 
  theme(legend.position="top", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_y_continuous(limits=c(0, 101)) + labs(title = "", x="", y="")


# ------------------------------------------------------------------------------------------------------------------------
# Find cells with CCA
# ------------------------------------------------------------------------------------------------------------------------
# Set up DropSeq object
# Find cells with Seurat batch effect correction
# Set up DropSeq object
dropseq <- CreateSeuratObject(counts = Ds_gem, project = "Dropseq", min.cells = 5)
dropseq$stim <- "Dropseq"
dropseq <- FindVariableFeatures(dropseq, selection.method = "vst", nfeatures = 2000)

# Set up SS2 object
ss2 <- CreateSeuratObject(counts = log2(ss2_gem+1), project = "SS2", min.cells = 5)
ss2$stim <- "SS2"
ss2 <- FindVariableFeatures(ss2, selection.method = "vst", nfeatures = 2000)

combined.anchors <- FindIntegrationAnchors(object.list = list(dropseq, ss2), dims = 1:20)
cells.combined <- IntegrateData(anchorset = combined.anchors, dims = 1:20)

DefaultAssay(cells.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
cells.combined <- ScaleData(cells.combined, verbose = FALSE)
cells.combined <- RunPCA(cells.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
#cells.combined <- RunUMAP(cells.combined, reduction = "pca", dims = 1:20)
cells.combined <- FindNeighbors(cells.combined, reduction = "pca", dims = 1:20)
cells.combined <- FindClusters(cells.combined, resolution = 0.5)

ncol <- unique(Idents(cells.combined))
gem_avg_ds <- matrix(NA, nrow = length(celltypes), ncol = length(ncol))
for (i in 1:length(ncol)) {
  cells <- intersect(WhichCells(cells.combined, idents = levels(ncol)[i]), colnames(ss2_gem))
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(ss2_gem[toupper(positive_markers$gene), cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg_ds[j,i] <- mean(na.omit(avg_exp[toupper(positive_markers$gene[which(positive_markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg_ds) <- celltypes
colnames(gem_avg_ds) <- paste("Cl", levels(ncol), sep = "_")

gem_avg_ds <- gem_avg_ds[, complete.cases(t(gem_avg_ds))]
gem_avg_norm <- apply(gem_avg_ds, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
CCA_ref_counts <- table(colSums(gem_thr))/length(unique(Idents(cells.combined)))
CCA_ref_counts <- c(1-sum(CCA_ref_counts), CCA_ref_counts)
names(CCA_ref_counts) <- c(0, names(CCA_ref_counts[-1]))

CCA_labels <- Idents(cells.combined)[colnames(ss2_gem)]
CCA_ARI <- c()
CCA_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(marker_based_labels[[i]])
  names(marker_labels) <- colnames(ss2_gem)
  CCA_ARI[i] <- mclust::adjustedRandIndex(marker_labels, CCA_labels)
  CCA_VI[i] <- mcclust::vi.dist(marker_labels, CCA_labels)
}


# ------------------------------------------------------------------------------------------------------------------------
# Find cells with MNN batch effect correction (Fig 4c right)
# ------------------------------------------------------------------------------------------------------------------------
mnn_res <- runMNN(target_gem = log(ss2_gem+1), reference_gem = Ds_gem)

mnn_labels <- mnn_res$target_labels

ncol <- levels(mnn_labels)
gem_avg <- matrix(NA, nrow = length(celltypes), ncol = length(ncol))
for (i in 1:length(ncol)) {
  cells <- names(mnn_labels)[which(mnn_labels == ncol[i])]
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(ss2_gem[toupper(positive_markers$gene), cells, drop = F])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(positive_markers$gene[which(positive_markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg) <- paste("GS", celltypes, sep = "_")
colnames(gem_avg) <- paste("Cl", ncol, sep = "_")

gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
MNN_ref_counts <- table(colSums(gem_thr))/length(ncol)
MNN_ref_counts <- c(1-sum(MNN_ref_counts), MNN_ref_counts)
names(MNN_ref_counts) <- c(0, names(MNN_ref_counts[-1]))

MNN_ARI <- c()
MNN_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(marker_based_labels[[i]])
  names(marker_labels) <- colnames(ss2_gem)
  MNN_ARI[i] <- mclust::adjustedRandIndex(marker_labels, mnn_labels)
  MNN_VI[i] <- mcclust::vi.dist(marker_labels, mnn_labels)
}


# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scID
# ------------------------------------------------------------------------------------------------------------------------
scID_res <- scid_multiclass(target_gem = ss2_gem, markers = markers, estimate_weights_from_target = T)
scID_labels <- scID_res$labels

# Get ARI and VI with marker based for different thresholds
scID_ARI <- c()
scID_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(marker_based_labels[[i]])
  names(marker_labels) <- colnames(ss2_gem)
  cells <- names(marker_labels)[which(!marker_labels %in% c("ambiguous", "orphan"))]
  scID_ARI[i] <- mclust::adjustedRandIndex(marker_labels[cells], scID_labels[cells])
  scID_VI[i] <- mcclust::vi.dist(marker_labels[cells], scID_labels[cells])
}
scID_ARI
m <- scID_res$markers
positive_markers <- m[which(m$avg_logFC > 0), ]
# average expression per celltype
gem_avg <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- na.omit(names(scID_labels)[which(scID_labels == celltypes[i])])
  if (length(cells) > 1) {
    avg_exp <- rowMeans(ss2_gem[toupper(positive_markers$gene), cells])
  } else if (length(cells) == 1) {
    avg_exp <- ss2_gem[toupper(positive_markers$gene), cells]
    names(avg_exp) <- toupper(positive_markers$gene)
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(positive_markers$gene[which(positive_markers$cluster == celltypes[j])])]))
  }
}

gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
scID_counts <- table(colSums(gem_thr))/length(unique(celltypes))
scID_counts <- c(1-sum(scID_counts), scID_counts)
names(scID_counts) <- c(0, names(scID_counts[-1]))


# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scmap
# ------------------------------------------------------------------------------------------------------------------------
scmap_labels <- run_scmap(reference_gem = Ds_gem, reference_labels = identities, target_gem = log(ss2_gem+1), n_features = 150)
table(scmap_labels)

# Get ARI and VI with marker based for different thresholds
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
    avg_exp <- rowMeans(ss2_gem[toupper(positive_markers$gene), cells])
  } else if (length(cells) == 1) {
    avg_exp <- ss2_gem[toupper(positive_markers$gene), cells]
    names(avg_exp) <- toupper(positive_markers$gene)
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(positive_markers$gene[which(positive_markers$cluster == celltypes[j])])]))
  }
}

gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
scmap_counts <- table(colSums(gem_thr))/length(unique(celltypes))
scmap_counts <- c(1-sum(scmap_counts), scmap_counts)
names(scmap_counts) <- c(0, names(scmap_counts[-1]))


# ------------------------------------------------------------------------------------------------------------------------
# find equivalent cells with CaSTLe
# ------------------------------------------------------------------------------------------------------------------------
castle_labels <- runCastle(log(ss2_gem+1), Ds_gem, identities)

# Get ARI and VI with marker based for different thresholds
castle_ARI <- c()
castle_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(marker_based_labels[[i]])
  names(marker_labels) <- colnames(ss2_gem)
  castle_ARI[i] <- mclust::adjustedRandIndex(marker_labels, castle_labels)
  castle_VI[i] <- mcclust::vi.dist(marker_labels, castle_labels)
}

# average expression per celltype
gem_avg <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- na.omit(names(castle_labels)[which(castle_labels == celltypes[i])])
  if (length(cells) > 1) {
    avg_exp <- rowMeans(ss2_gem[toupper(positive_markers$gene), cells])
  } else if (length(cells) == 1) {
    avg_exp <- ss2_gem[toupper(positive_markers$gene), cells]
    names(avg_exp) <- toupper(positive_markers$gene)
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(positive_markers$gene[which(positive_markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg) <- celltypes
colnames(gem_avg) <- celltypes

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
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "", x="", y="")


# ------------------------------------------------------------------------------------------------------------------------
# plot number of signatures per cluster 
# ------------------------------------------------------------------------------------------------------------------------
df <- data.frame(counts = c(Seurat_target_counts, scID_counts, scmap_counts, castle_counts, Seurat_ref_counts, MNN_ref_counts, CCA_ref_counts),
                 numbers = c(names(Seurat_target_counts), names(scID_counts), names(scmap_counts), names(castle_counts), 
                             names(Seurat_ref_counts),  names(MNN_ref_counts), names(CCA_ref_counts)), 
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

