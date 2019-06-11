library(Seurat)
library(dplyr)
library(pheatmap)
library(scran)
library(scmap)
library(ggsci)
library(ggplot2)
library(scID)
source("~/ownCloud/DocSyncUoE/scID/scripts/resubmission_scripts_May2019/run_other_methods.R")

# Read data
Ds_gem <- readRDS("~/scID_manuscript_figures/data/Figure2/Reference_gem.rds")
ss2_gem <- readRDS("~/scID_manuscript_figures/data/Figure2/Target_gem.rds")
identities <- readRDS("~/scID_manuscript_figures/data/Figure2/Reference_clusters.rds")

celltypes <- c("RBC (Rod Bipolar cell)", "MG (Mueller Glia)", "BC5A (Cone Bipolar cell 5A)", 
               "BC7 (Cone Bipolar cell 7)", "BC6", "BC5C", "BC1A", "BC3B", "BC1B", 
               "BC2", "BC5D", "BC3A", "BC5B", "BC4", "BC8/9 (mixture of BC8 and BC9)", 
               "AC (Amacrine cell)", "Rod Photoreceptors", "Cone Photoreceptors")
# ------------------------------------------------------------------------------------------------------------------------
# Find markers of each cluster from DropSeq or read pre-computed markers
so_ds <- CreateSeuratObject(counts = Ds_gem)
so_ds <- ScaleData(so_ds)
Idents(so_ds) <- factor(identities, levels = celltypes)

markers <- FindAllMarkers(so_ds, only.pos = FALSE, test.use = "MAST", logfc.threshold = 0.5)
#markers <- read.delim("~/scID_manuscript_figures/data/Figure2/markers.txt", stringsAsFactors = F)

positive_markers <- markers[which(markers$avg_logFC > 0), ]

gem_avg_ds <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- WhichCells(so_ds, idents = celltypes[i])
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
colnames(gem_avg_ds) <- celltypes

# Figure 2B left
pheatmap(gem_avg_ds, border="white", color = colorspace::diverge_hsv(50), cluster_rows = F,
         cluster_cols = F, show_colnames = F, fontsize_row = 5, border_color = F,show_rownames = F,
         scale = "row", fontsize_col = 5, legend = F,
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/DropSeq_heatmap.pdf", 
         width = 2, height = 2)

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
gene_names <- c()
top_markers <- c()
for (celltype in celltypes) {
  positive_markers <- intersect(markers$gene[which(markers$cluster == celltype)], markers$gene[which(markers$avg_logFC > 0)])
  weights <- scID_weight(Ds_gem[positive_markers, ], true_cells = names(which(identities == celltype)), false_cells = names(which(identities != celltype)))
  gene_names <- c(top2_markers, names(sort(weights, decreasing = T)[1:2]))
  top_markers <- c(top_markers, rep(celltype, 2))
}
names(top_markers) <- gene_names

# Normalize data by 99th percentile
ss2_gem_norm <- t(apply(ss2_gem[names(top_markers), ], 1, function(x) scID::normalize_gene(x)))

thresholds <- c(0.1, 0.25, 0.5, 0.75)
data <- data.frame(matrix(NA, length(thresholds), 3))
biomarker_labels <- list()
for (i in 1:length(thresholds)) {
  biomarker_labels[[i]] <- list()
  number_of_celltypes <- c()
  for (cell in colnames(ss2_gem_norm)) {
    ct <- c()
    for (gene in rownames(ss2_gem_norm)) {
      if (ss2_gem_norm[gene, cell] > thresholds[i]) {
        ct <- c(ct, top_markers[gene])
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
  data[i, 1] <- round(sum(number_of_celltypes == 0)*100/ncol(ss2_gem_norm), 2)
  data[i, 2] <- round(sum(number_of_celltypes == 1)*100/ncol(ss2_gem_norm), 2)
  data[i, 3] <- round(sum(number_of_celltypes > 1)*100/ncol(ss2_gem_norm), 2)
}


colnames(data) <- c("Orphan", "Identified", "Ambiguous")
data$threshold <- factor(thresholds, levels = thresholds)

# Figure 2C
df <- reshape2::melt(data, id.vars = c("threshold"))
#pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/marker_based_identification_SS2.pdf", width = 2, height = 2.5)
ggplot(df, aes(fill=variable, y=value, x=threshold)) + 
  geom_bar(stat="identity") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits=c(0, 101)) + labs(title = "", x="", y="")
dev.off()

# ------------------------------------------------------------------------------------------------------------------------
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

ncol <- length(unique(Idents(cells.combined)))
gem_avg_ds <- matrix(NA, nrow = length(celltypes), ncol = ncol)
for (i in 1:ncol) {
  cells <- intersect(WhichCells(cells.combined, idents = i-1), colnames(Ds_gem))
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
colnames(gem_avg_ds) <- paste("Cl", unique(Idents(cells.combined)), sep = "_")

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
# Find cells with MNN batch effect correction
common_genes <- intersect(rownames(Ds_gem), rownames(ss2_gem))
mnn_corrected <- mnnCorrect(as.matrix(Ds_gem[common_genes, ]), as.matrix(log2(ss2_gem[common_genes, ]+1)))

# Cluster mnn_corrected gem with Seurat
gem <- as.data.frame(mnn_corrected$corrected)
colnames(gem)[grep("X", colnames(gem))] <- colnames(ss2_gem)

seurat <- CreateSeuratObject(counts = gem)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat,  do.print = FALSE)
seurat <- FindNeighbors(seurat)
seurat <- FindClusters(seurat)

ncol <- unique(Idents(seurat))
gem_avg <- matrix(NA, nrow = length(celltypes), ncol = length(ncol))
for (i in 1:length(ncol)) {
  cells <- intersect(WhichCells(seurat, idents = levels(ncol)[i]), colnames(Ds_gem))
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
colnames(gem_avg) <- paste("Cl", levels(ncol), sep = "_")

# Figure 2D right
pheatmap(gem_avg, border="white", color = colorspace::diverge_hsv(50), cluster_rows = F,
         cluster_cols = F, show_colnames = F, fontsize_row = 5, border_color = F,
         show_rownames = F, fontsize = 5, scale = "row", legend = F, width = 2, height = 2, 
         filename = "~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure2/MNN_heatmap.pdf")

gem_avg_ds <- gem_avg_ds[,colSums(is.na(gem_avg_ds))<nrow(gem_avg_ds)]
gem_avg_norm <- apply(gem_avg_ds, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
MNN_ref_counts <- table(colSums(gem_thr))/length(unique(Idents(seurat)))
MNN_ref_counts <- c(1-sum(MNN_ref_counts), MNN_ref_counts)
names(MNN_ref_counts) <- c(0, names(MNN_ref_counts[-1]))

rm(seurat)

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scID
# remove bulk cells
scID_output <- scid_multiclass(target_gem = ss2_gem, reference_gem = Ds_gem, reference_clusters = identities, 
                               estimate_weights_from_target = TRUE, normalize_reference = F, only_pos = T, logFC = 0.2)
scID_labels <- scID_output$labels
table(scID_labels)

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


#positive_markers <- markers[which(markers$avg_logFC > 0), ]
positive_markers <- scID_output$markers
# average expression per celltype
gem_avg <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- na.omit(names(scID_labels)[which(scID_labels == celltypes[i])])
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
colnames(gem_avg) <- celltypes

# Figure 2F right
pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, 
         border_color = F, show_rownames = T, scale = "row")

# Remove NA columns
gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
scID_counts <- table(colSums(gem_thr))/length(unique(celltypes))
scID_counts <- c(1-sum(scID_counts), scID_counts)
names(scID_counts) <- c(0, names(scID_counts[-1]))

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scmap
scmap_labels <- run_scmap(reference_gem = Ds_gem, reference_labels = identities, target_gem = log(ss2_gem+1), n_features = 150)
table(scmap_labels)

# Get ARI with marker based for different thresholds
scmap_ARI <- c()
scmap_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(biomarker_labels[[i]])
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
  marker_labels <- unlist(biomarker_labels[[i]])
  names(marker_labels) <- colnames(ss2_gem)
  castle_ARI[i] <- mclust::adjustedRandIndex(marker_labels, castle_labels)
  castle_VI[i] <- mcclust::vi.dist(marker_labels, castle_labels)
}

# average expression per celltype
gem_avg <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- na.omit(names(castle_labels)[which(castle_labels == celltypes[i])])
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
gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
castle_counts <- table(colSums(gem_thr))/length(unique(celltypes))
castle_counts <- c(1-sum(castle_counts), castle_counts)
names(castle_counts) <- c(0, names(castle_counts[-1]))

# ---------------------------------------------------------------
# plot ARI - Figure 2 H
df <- data.frame(value=c(scID_ARI, scmap_ARI, castle_ARI), 
                 threshold=c(thresholds, thresholds, thresholds), 
                 method=c(rep("scID", length(thresholds)), rep("scmap", length(thresholds)), rep("CaSTLe", length(thresholds))))
df$threshold <- factor(df$threshold, levels = thresholds)

ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")

# ---------------------------------------------------------------
# plot VI - Figure 2 I
df <- data.frame(value=c(scID_VI, scmap_VI, castle_VI), 
                 threshold=c(thresholds, thresholds, thresholds), 
                 method=c(rep("scID", length(thresholds)), rep("scmap", length(thresholds)), rep("CaSTLe", length(thresholds))))
df$threshold <- factor(df$threshold, levels = thresholds)

ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm")) + 
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
df <- data.frame(counts = c(Seurat_target_counts, scID_counts, scmap_counts, castle_counts, Seurat_ref_counts),
                 numbers = c(names(Seurat_target_counts), names(scID_counts), names(scmap_counts), names(castle_counts), names(Seurat_ref_counts)), 
                 method = c(rep("Seurat", length(Seurat_target_counts)), 
                            rep("scID", length(scID_counts)), 
                            rep("scmap", length(scmap_counts)), 
                            rep("CaSTLe", length(castle_counts)), 
                            rep("reference", length(Seurat_ref_counts))))

ggplot(df, aes(y=counts, x=method, fill=numbers)) + 
  geom_bar(stat="identity") + scale_fill_grey(start = 1, end = 0) +
  theme(legend.position="top", text = element_text(size=5), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")
