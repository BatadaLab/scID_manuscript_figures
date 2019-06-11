library(Seurat)
library(dplyr)
library(pheatmap)
library(scran)
library(scmap)
library(ggsci)
library(scID)
source("~/ownCloud/DocSyncUoE/scID/scripts/resubmission_scripts_May2019/run_other_methods.R")

# Read data
Ds_gem <- loadfast("~/Google Drive/Data/singlecell/published/mice/shekhar2016_retina/exp_matrix.txt")
rownames(Ds_gem) <- toupper(rownames(Ds_gem))

ss2_gem <- loadfast("~/Google Drive/Data/singlecell/published/mice/shekhar2016_retina/GSE80232_vsx2.RSEM.genes.tpm.matrix.filtered.txt")
ss2_gem <- ss2_gem[, -grep("bulk", colnames(ss2_gem))]

# Read labels of DropSeq data
dropseq_groups <- read.delim("~/Google Drive/Data/singlecell/published/mice/shekhar2016_retina/clust_retinal_bipolar.txt", stringsAsFactors = F, header = T, row.names = 1)
dropseq_groups <- dropseq_groups[-1, ]
doublets <- rownames(dropseq_groups)[which(dropseq_groups$CLUSTER == "Doublets/Contaminants")]
dropseq_groups <- dropseq_groups[setdiff(rownames(dropseq_groups), doublets), ]

Ds_gem <- Ds_gem[, setdiff(colnames(Ds_gem), doublets)]

# ------------------------------------------------------------------------------------------------------------------------
# Find markers of each cluster from DropSeq (Fig 4b left)
so_ds <- CreateSeuratObject(counts = Ds_gem)
so_ds <- NormalizeData(so_ds)
so_ds <- ScaleData(so_ds)

identities <- dropseq_groups[colnames(Ds_gem), "CLUSTER"]
names(identities) <- colnames(Ds_gem)
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
         cluster_cols = F, show_colnames = F, fontsize_row = 3, border_color = F,show_rownames = F,
         scale = "row", , width = 2, height = 2, legend = F,
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/DropSeq_heatmap.pdf")
dev.off()

gem_avg_norm <- apply(gem_avg_ds, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
Seurat_ref_counts <- table(colSums(gem_thr))/length(unique(Idents(so_ds)))
Seurat_ref_counts <- c(1-sum(Seurat_ref_counts), Seurat_ref_counts)
names(Seurat_ref_counts) <- c(0, names(Seurat_ref_counts[-1]))

# ------------------------------------------------------------------------------------------------------------------------
# Cluster SS2 data with Seurat
SS2_clustered <- CreateSeuratObject(counts = log2(ss2_gem+1))
SS2_clustered <- FindVariableFeatures(SS2_clustered, do.contour=F)
SS2_clustered <- ScaleData(SS2_clustered)
SS2_clustered <- RunPCA(SS2_clustered,  do.print = T)
SS2_clustered <- FindNeighbors(SS2_clustered)
SS2_clustered <- FindClusters(SS2_clustered)
SS2_clustered <- RunTSNE(SS2_clustered)

pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/SS2_tsne.pdf", width = 2, height = 2)
TSNEPlot(SS2_clustered, pt.size = 0.05, label = T) + NoLegend() + NoAxes() 
dev.off()

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
         show_rownames = F, scale = "row", width = 2, height = 2, legend = F,
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/SS2_Seurat_heatmap.pdf")

gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
Seurat_target_counts <- table(colSums(gem_thr))/length(unique(Idents(SS2_clustered)))
Seurat_target_counts <- c(1-sum(Seurat_target_counts), Seurat_target_counts)
names(Seurat_target_counts) <- c(0, names(Seurat_target_counts[-1]))

rm(SS2_clustered)

# -----------------------------------------------------------------------------------------------------------------------
# Marker based identification
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
pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/marker_based_identification_SS2.pdf", width = 2, height = 2.5)
ggplot(df, aes(fill=variable, y=value, x=threshold)) + 
  geom_bar(stat="identity") + scale_fill_jco() + 
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_y_continuous(limits=c(0, 101)) + labs(title = "", x="", y="")
dev.off()

# ------------------------------------------------------------------------------------------------------------------------
# Find cells with Seurat batch effect correction (Fig 4c left)
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
cells.combined <- RunUMAP(cells.combined, reduction = "pca", dims = 1:20)
cells.combined <- FindNeighbors(cells.combined, reduction = "pca", dims = 1:20)
cells.combined <- FindClusters(cells.combined, resolution = 0.5)

ncol <- unique(Idents(cells.combined))
gem_avg_ds <- matrix(NA, nrow = length(celltypes), ncol = length(ncol))
for (i in 1:length(ncol)) {
  cells <- intersect(WhichCells(cells.combined, idents = levels(ncol)[i]), colnames(Ds_gem))
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(Ds_gem[toupper(positive_markers$gene), cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg_ds[j,i] <- mean(na.omit(avg_exp[toupper(positive_markers$gene[which(positive_markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg_ds) <- celltypes
colnames(gem_avg_ds) <- paste("Cl", levels(ncol), sep = "_")

# Figure 2D left
pheatmap(gem_avg_ds, border="white", color = colorspace::diverge_hsv(50), cluster_rows = F,
         cluster_cols = F, show_colnames = F, fontsize_row = 5, border_color = F,
         show_rownames = F, fontsize_col = 5, scale = "row", width = 2, height = 2, legend = F,
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/CCA_heatmap.pdf")

gem_avg_ds <- gem_avg_ds[,colSums(is.na(gem_avg_ds))<nrow(gem_avg_ds)]

gem_avg_norm <- apply(gem_avg_ds, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
CCA_ref_counts <- table(colSums(gem_thr))/length(unique(Idents(cells.combined)))
CCA_ref_counts <- c(1-sum(CCA_ref_counts), CCA_ref_counts)
names(CCA_ref_counts) <- c(0, names(CCA_ref_counts[-1]))

rm(cells.combined)

# ------------------------------------------------------------------------------------------------------------------------
# Find cells with MNN batch effect correction (Fig 4c right)
mnn_res <- runMNN(target_gem = log(ss2_gem+1), reference_gem = Ds_gem)

mnn_labels <- mnn_res$reference_labels

ncol <- levels(mnn_labels)
gem_avg <- matrix(NA, nrow = length(celltypes), ncol = length(ncol))
for (i in 1:length(ncol)) {
  cells <- names(mnn_labels)[which(mnn_labels == ncol[i])]
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(Ds_gem[toupper(positive_markers$gene), cells, drop = F])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(positive_markers$gene[which(positive_markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg) <- paste("GS", celltypes, sep = "_")
colnames(gem_avg) <- paste("Cl", ncol, sep = "_")
pheatmap(gem_avg, border="white", color = colorspace::diverge_hsv(50), cluster_rows = F,
         cluster_cols = F, show_colnames = F, fontsize_row = 10, border_color = F,show_rownames = F,
         scale = "row", legend = F, width = 2, height = 2, 
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/MNN_heatmap.pdf")

gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
MNN_ref_counts <- table(colSums(gem_thr))/length(ncol)
MNN_ref_counts <- c(1-sum(MNN_ref_counts), MNN_ref_counts)
names(MNN_ref_counts) <- c(0, names(MNN_ref_counts[-1]))

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scID
# Below is working well
# scID_res <- scid_multiclass(target_gem = ss2_gem, reference_gem = Ds_gem, reference_clusters = identities, 
#                             normalize_reference = T, estimate_weights_from_target = T, logFC = 0.9, only_pos = F)
scID_res <- scid_multiclass(target_gem = ss2_gem, markers = markers, estimate_weights_from_target = T)
scID_labels <- scID_res$labels
table(scID_labels)

# Get ARI and VI (variation of information) with marker based for different thresholds
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

pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), cluster_rows = F,
         cluster_cols = F, show_colnames = F, fontsize_row = 10, border_color = F, show_rownames = F,
         scale = "row", height = 2, width = 2, legend = F,
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/scID_heatmap.pdf")
dev.off()
gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
scID_counts <- table(colSums(gem_thr))/length(unique(celltypes))
scID_counts <- c(1-sum(scID_counts), scID_counts)
names(scID_counts) <- c(0, names(scID_counts[-1]))

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scmap (Fig 4e right)
scmap_labels <- run_scmap(reference_gem = Ds_gem, reference_labels = identities, target_gem = log(ss2_gem+1), n_features = 150)
table(scmap_labels)

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

pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), cluster_rows = F,
         cluster_cols = F, show_colnames = F, fontsize_row = 3, border_color = F,show_rownames = F,
         scale = "row", height = 2, width = 2, legend = F,
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/scmap_heatmap.pdf")
dev.off()

gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
scmap_counts <- table(colSums(gem_thr))/length(unique(celltypes))
scmap_counts <- c(1-sum(scmap_counts), scmap_counts)
names(scmap_counts) <- c(0, names(scmap_counts[-1]))

# ---------------------------------------------------------------
# find equivalent cells with CaSTLe
castle_labels <- runCastle(log(ss2_gem+1), Ds_gem, identities)

# Get ARI and VI (variation of information) with marker based for different thresholds
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

# Figure 2F right
pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = F, fontsize = 5, 
         border_color = F, show_rownames = F, scale = "row", legend = F, width = 2, height = 2, 
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/castle_heatmap.pdf")

# Remove NA columns
gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
castle_counts <- table(colSums(gem_thr))/length(unique(celltypes))
castle_counts <- c(1-sum(castle_counts), castle_counts)
names(castle_counts) <- c(0, names(castle_counts[-1]))

# ---------------------------------------------------------------
# plot ARI
df <- data.frame(value=c(scID_ARI, scmap_ARI, castle_ARI), 
                 threshold=c(thresholds, thresholds, thresholds), 
                 method=c(rep("scID", length(thresholds)), rep("scmap", length(thresholds)), rep("CaSTLe", length(thresholds))))
df$threshold <- factor(df$threshold, levels = thresholds)

pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/ARI.pdf", width = 1.5, height = 1.5)
ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "", x="", y="")
dev.off()

# ---------------------------------------------------------------
# plot VI
df <- data.frame(value=c(scID_VI, scmap_VI, castle_VI), 
                 threshold=c(thresholds, thresholds, thresholds), 
                 method=c(rep("scID", length(thresholds)), rep("scmap", length(thresholds)), rep("CaSTLe", length(thresholds))))
df$threshold <- factor(df$threshold, levels = thresholds)

pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/VI.pdf", width = 1.5, height = 1.5)
ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "", x="", y="")
dev.off()

# ---------------------------------------------------------------
# plot number of signatures per cluster (reference data)
df <- data.frame(counts = c(Seurat_ref_counts, MNN_ref_counts, CCA_ref_counts),
                 numbers = c(names(Seurat_ref_counts), names(MNN_ref_counts), names(CCA_ref_counts)), 
                 method = c(rep("Seurat", length(Seurat_ref_counts)), 
                            rep("MNN", length(MNN_ref_counts)), 
                            rep("CCA", length(CCA_ref_counts))))

pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/Num_signatures_alignment.pdf", width = 1.5, height = 2)
ggplot(df, aes(y=counts, x=method, fill=numbers)) + 
  geom_bar(stat="identity") + scale_fill_grey(start = 1, end = 0) +
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "", x="", y="")
dev.off()

# ---------------------------------------------------------------
# plot number of signatures per cluster (target data)
df <- data.frame(counts = c(Seurat_target_counts, scID_counts, scmap_counts, castle_counts, Seurat_ref_counts),
                 numbers = c(names(Seurat_target_counts), names(scID_counts), names(scmap_counts), names(castle_counts), names(Seurat_ref_counts)), 
                 method = c(rep("Seurat", length(Seurat_target_counts)), 
                            rep("scID", length(scID_counts)), 
                            rep("scmap", length(scmap_counts)), 
                            rep("CaSTLe", length(castle_counts)), 
                            rep("reference", length(Seurat_ref_counts))))

pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/Num_signatures_mapping.pdf", width = 2, height = 2)
ggplot(df, aes(y=counts, x=method, fill=numbers)) + 
  geom_bar(stat="identity") + scale_fill_grey(start = 1, end = 0) +
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "", x="", y="")
dev.off()

