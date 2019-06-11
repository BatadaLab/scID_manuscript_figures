library(Seurat) # version 2.3.1
library(dplyr)
library(pheatmap)
library(scran)
library(scmap)
library(dplyr)
library(ggsci)
library(scID)
library(ggplot2)
library(ggpubr)

# --------------------------------------------------------------------------
counts_to_tpm <-function(GEM_umi_counts){
  library(scater)
  
  ## REMOVE ALL ROWS WITH ZERO VALUES
  IDX_ZEROS=apply(GEM_umi_counts, 1, function(row) all(row==0))
  GEM_umi_counts=GEM_umi_counts[!IDX_ZEROS,]
  message(sprintf("After discarding all zero rows, left with %s rows.", dim(GEM_umi_counts)[1]))
  
  counts_scater <- SingleCellExperiment(assays = list(counts = as.matrix(GEM_umi_counts)))
  n_genes=nrow(GEM_umi_counts)
  effective_length <- rep(1000, n_genes) # for paired end reads and tagmentation, the averge bioanalyzer fragment size is 1000
  GEM_tpm = calculateTPM(counts_scater, effective_length)
  
  return(GEM_tpm)
}

# ---------------------------------------------------------------------------
# Read data
E18_9k <- readRDS("~/scID_manuscript_figures/data/Figure3/reference_gem.rds")
E18_nuclei <- readRDS("~/scID_manuscript_figures/data/Figure3/target_gem.rds")

# TPM normalization
E18_9k_logtpm <- counts_to_tpm(E18_9k)
E18_nuclei_logtpm <- counts_to_tpm(E18_nuclei)

# Edit colnames of nuclei data to distinguish from E18_9k
colnames(E18_nuclei_logtpm) <- paste(colnames(E18_nuclei_logtpm), "nuc", sep = "_")
colnames(E18_nuclei) <- paste(colnames(E18_nuclei), "nuc", sep = "_")

# ------------------------------------------------------------------------------------------------------------------------
# Cluster E18_9k dataset with Seurat
sobj_9k <- CreateSeuratObject(raw.data = E18_9k)
sobj_9k <- NormalizeData(sobj_9k)
sobj_9k <- ScaleData(sobj_9k)
sobj_9k <- FindVariableGenes(sobj_9k)
sobj_9k <- RunPCA(sobj_9k,  do.print = FALSE)
sobj_9k <- ProjectPCA(sobj_9k)
PCElbowPlot(sobj_9k)
sobj_9k <- FindClusters(sobj_9k, dims.use = 1:5, save.SNN = F)
sobj_9k <- RunTSNE(sobj_9k, dims.use = 1:5, do.fast = F)

# Figure 3A left
TSNEPlot(sobj_9k, do.label = T, pt.size = 0.05, no.axes = T, no.legend = T)

# Find markers or use pre-computed list
# markers <- FindAllMarkers(sobj_9k, only.pos = TRUE, test.use = "MAST", logfc.threshold = 0.5, do.print = FALSE)
markers <- readRDS("~/scID_manuscript_figures/data/Figure3/markers.rds")

celltypes <- unique(sobj_9k@ident)

gem_avg_9k <- matrix(NA, length(unique(sobj_9k@ident)), length(unique(sobj_9k@ident)))
for (i in 1:length(unique(sobj_9k@ident))) {
  cells <- WhichCells(sobj_9k, i-1)
  if (length(unique(sobj_9k@ident)) > 1) {
    avg_exp <- rowMeans(E18_9k_logtpm[markers$gene, cells])
  } else if (length(cells) == 1) {
    avg_exp <- gem[markers$gene, cells]
    names(avg_exp) <- markers$gene
  } else {
    next
  }
  for (j in 1:length(unique(sobj_9k@ident))) {
    gem_avg_9k[j,i] <- mean(na.omit(avg_exp[markers$gene[which(markers$cluster == j-1)]]))
  }
}
rownames(gem_avg_9k) <- paste("GS", celltypes, sep = "_")
colnames(gem_avg_9k) <- paste("Cl", celltypes, sep = "_")

# Figure 3B left
pheatmap(gem_avg_9k, border="white", color = colorspace::diverge_hsv(50), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, 
         border_color = F, show_rownames = T, scale = "row")

pct_cells <- matrix(NA, length(unique(sobj_9k@ident)), length(unique(sobj_9k@ident)))
for (c in 1:length(celltypes)) {
  cells <- WhichCells(sobj_9k, celltypes[c])
  for (s in 1:length(celltypes)) {
    avg_exp <- colMeans(E18_9k_logtpm[markers[which(markers$cluster == celltypes[s]), "gene"], ])
    thr <- quantile(avg_exp, 0.95)
    pct_cells[s,c] <- round(length(intersect(names(which(avg_exp > thr)), cells))*100/length(cells), 2)
  }
}

gem_thr <- apply(pct_cells, 2, function(x) sum(x>=100/3))
Seurat_ref_counts <- table(gem_thr)/length(unique(sobj_9k@ident))

# -----------------------------------------------------------------------------------------------------------------------
# Marker based identification
markers_filt <- markers[-which(duplicated(markers$gene)), ]
top_markers <- markers_filt %>% group_by(cluster) %>% top_n(2, avg_logFC)
top_markers <- top_markers[which(top_markers$gene %in% rownames(E18_nuclei_logtpm)), ]
rownames(top_markers) <- top_markers$gene

# Normalize data by 99th percentile
nuclei_norm <- t(apply(E18_nuclei_logtpm[intersect(top_markers$gene, rownames(E18_nuclei_logtpm)), ], 1, function(x) normalize_gem(x)))

thresholds <- c(0, 0.25, 0.5, 0.75)
data <- data.frame(matrix(NA, length(thresholds), 3))
biomarker_labels <- list()
for (i in 1:length(thresholds)) {
  biomarker_labels[[i]] <- list()
  for (cell in colnames(nuclei_norm)) {
    expressed_genes <- c()
    for (gene in rownames(nuclei_norm)) {
      if (nuclei_norm[gene, cell] > quantile(unlist(nuclei_norm[gene, ]), thresholds[i])) {
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
  data[i, 1] <- round(sum(biomarker_labels[[i]] == "orphan")*100/ncol(nuclei_norm), 2)
  data[i, 3] <- round(sum(biomarker_labels[[i]] == "ambiguous")*100/ncol(nuclei_norm), 2)
  data[i, 2] <- 100 - (data[i,1] + data[i,3])
}
colnames(data) <- c("Orphan", "Identified", "Ambiguous")
data$threshold <- thresholds

# Figure 3C
df <- reshape2::melt(data, id.vars = c("threshold"))
ggplot(df, aes(fill=variable, y=value, x=threshold)) + 
  geom_bar(stat="identity") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=5), plot.margin = unit(c(0,0,0,0), "cm")) +labs(title = "", x="", y="")

# ------------------------------------------------------------------------------------------------------------------------
# Cluster nuclei data with Seurat 
nuclei_clustered <- CreateSeuratObject(raw.data = E18_nuclei)
nuclei_clustered <- NormalizeData(nuclei_clustered)
nuclei_clustered <- ScaleData(nuclei_clustered)
nuclei_clustered <- FindVariableGenes(nuclei_clustered, do.contour=F)
nuclei_clustered <- RunPCA(nuclei_clustered,  do.print = T)
nuclei_clustered <- ProjectPCA(nuclei_clustered)
PCElbowPlot(nuclei_clustered)
dims <- 10
nuclei_clustered <- FindClusters(nuclei_clustered, dims.use = 1:dims)
nuclei_clustered <- RunTSNE(nuclei_clustered, dims.use = 1:dims, do.fast = F)

# Figure 3A right
TSNEPlot(nuclei_clustered, do.label = T, pt.size = 0.05, no.axes = T, no.legend = T)

markers_filt <- markers[which(markers$gene %in% rownames(E18_nuclei_logtpm)), ]
ncol <- length(unique(nuclei_clustered@ident))
nuclei_avg <- matrix(NA, length(unique(sobj_9k@ident)), ncol)
for (i in 1:ncol) {
  cells <- WhichCells(nuclei_clustered, i-1)
  if (length(cells) > 1) {
    avg_exp <- rowMeans(E18_nuclei_logtpm[markers_filt$gene, cells])
  } else if (length(cells) == 1) {
    avg_exp <- E18_nuclei_logtpm[markers_filt$gene, cells]
    names(avg_exp) <- markers_filt$gene
  } else {
    next
  }
  for (j in 1:length(unique(sobj_9k@ident))) {
    nuclei_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == j-1)]]))
  }
}
rownames(nuclei_avg) <- paste("GS", celltypes, sep = "_")
colnames(nuclei_avg) <- paste("Cl", unique(sobj_9k@ident), sep = "_")

# Figure 3B right
pheatmap(nuclei_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, border_color = F,
         show_rownames = T, scale = "row")

pct_cells <- matrix(NA, length(unique(sobj_9k@ident)), length(unique(nuclei_clustered@ident)))
for (c in 1:length(unique(nuclei_clustered@ident))) {
  cells <- WhichCells(nuclei_clustered, unique(nuclei_clustered@ident)[c])
  for (s in 1:length(celltypes)) {
    avg_exp <- colMeans(E18_nuclei_logtpm[markers_filt[which(markers_filt$cluster == celltypes[s]), "gene"], ])
    thr <- quantile(avg_exp, 0.95)
    pct_cells[s,c] <- round(length(intersect(names(which(avg_exp > thr)), cells))*100/length(cells), 2)
  }
}

gem_thr <- apply(pct_cells, 2, function(x) sum(x>=100/3))
Seurat_target_counts <- table(gem_thr)/length(unique(nuclei_clustered@ident))

# ------------------------------------------------------------------------------------------------------------------------
# Find cells with Seurat batch effect correction 
# Set up DropSeq object
ref <- CreateSeuratObject(raw.data = E18_9k, project = "Reference")
ref@meta.data$stim <- "Reference"
ref <- FilterCells(ref, subset.names = "nGene")#, low.thresholds = 500, high.thresholds = Inf)
ref <- NormalizeData(ref)
ref <- ScaleData(ref, display.progress = F)
# Set up SS2 object
nuclei <- CreateSeuratObject(raw.data = E18_nuclei, project = "nuclei")#, min.cells = 5)
nuclei@meta.data$stim <- "nuclei"
nuclei <- FilterCells(nuclei, subset.names = "nGene")#, low.thresholds = 500, high.thresholds = Inf)
nuclei <- NormalizeData(nuclei)
nuclei <- ScaleData(nuclei, display.progress = F)
# Gene selection for input to CCA
ref <- FindVariableGenes(ref, do.plot = F)
nuclei <- FindVariableGenes(nuclei, do.plot = F)
g.1 <- head(rownames(ref@hvg.info), 1000)
g.2 <- head(rownames(nuclei@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ref@scale.data))
genes.use <- intersect(genes.use, rownames(nuclei@scale.data))
# Run CCA to combine datasets
cells.combined <- RunCCA(ref, nuclei, genes.use = genes.use, num.cc = 20)
p3 <- MetageneBicorPlot(cells.combined, grouping.var = "stim", dims.eval = 1:20, display.progress = T)
# Align the CCA subspaces
cells.combined <- AlignSubspace(cells.combined, reduction.type = "cca", grouping.var = "stim", dims.align = 1:10)
# Cluster combined dataset
cells.combined <- FindClusters(cells.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:10)

ncol <- length(unique(cells.combined@ident))
gem_avg_ds <- matrix(NA, nrow = length(unique(sobj_9k@ident)), ncol = ncol)
gem <- cells.combined@scale.data
for (i in 1:ncol) {
  cells <- intersect(WhichCells(cells.combined, i-1), colnames(E18_9k_logtpm))
  if (length(cells) > 1) {
    avg_exp <- rowMeans(gem[markers$gene, cells])
  } else if (length(cells) == 1) {
    avg_exp <- gem[markers$gene, cells]
    names(avg_exp) <- markers$gene
  } else {
    next
  }
  for (j in 1:length(unique(sobj_9k@ident))) {
    gem_avg_ds[j,i] <- mean(na.omit(avg_exp[markers$gene[which(markers$cluster == j-1)]]))
  }
}
rownames(gem_avg_ds) <- paste("GS", celltypes, sep="_")
colnames(gem_avg_ds) <- paste("Cl", unique(cells.combined@ident), sep="_")

# Figure 3D left
pheatmap(gem_avg_ds, border="white", color = colorspace::diverge_hsv(50), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, 
         border_color = F, show_rownames = T, scale = "row")

pct_cells <- matrix(NA, nrow = length(unique(sobj_9k@ident)), ncol = ncol)
for (c in 1:ncol) {
  cells <- intersect(WhichCells(cells.combined, i-1), colnames(E18_9k_logtpm))
  for (s in 1:length(celltypes)) {
    avg_exp <- colMeans(E18_9k_logtpm[markers[which(markers$cluster == celltypes[s]), "gene"], ])
    thr <- quantile(avg_exp, 0.95)
    pct_cells[s,c] <- round(length(intersect(names(which(avg_exp > thr)), cells))*100/length(cells), 2)
  }
}

gem_thr <- apply(pct_cells, 2, function(x) sum(x>=100/3))
CCA_ref_counts <- table(gem_thr)/length(unique(cells.combined@ident))

# ------------------------------------------------------------------------------------------------------------------------
# Find cells with MNN batch effect correction
common_genes <- intersect(rownames(E18_9k_logtpm), rownames(E18_nuclei_logtpm))
mnn_corrected <- mnnCorrect(as.matrix(E18_9k_logtpm[common_genes, ]), as.matrix(E18_nuclei_logtpm[common_genes, ]))

# Cluster mnn_corrected gem with Seurat
gem <- as.data.frame(mnn_corrected$corrected)
colnames(gem)[grep("X", colnames(gem))] <- colnames(E18_nuclei_logtpm)

seurat <- CreateSeuratObject(raw.data = gem)
seurat <- NormalizeData(seurat)
seurat <- ScaleData(seurat)
seurat <- FindVariableGenes(seurat, x.low.cutoff = 0.05, x.high.cutoff = 1, y.cutoff = 0.5, do.contour = F)
seurat <- RunPCA(seurat,  do.print = FALSE)
seurat <- ProjectPCA(seurat)
PCElbowPlot(seurat)
seurat <- FindClusters(seurat, dims.use = 1:7, save.SNN = F)

ncol <- unique(seurat@ident)
gem_avg <- matrix(NA, nrow = length(unique(sobj_9k@ident)), ncol = length(ncol))
for (i in 1:length(ncol)) {
  cells <- intersect(WhichCells(seurat, ncol[i]), colnames(E18_9k_logtpm))
  if (length(cells) > 1) {
    avg_exp <- rowMeans(gem[markers$gene, cells])
  } else if (length(cells) == 1) {
    avg_exp <- gem[markers$gene, cells]
    names(avg_exp) <- markers$gene
  } else {
    next
  }
  for (j in 1:length(unique(sobj_9k@ident))) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[markers$gene[which(markers$cluster == j-1)]]))
  }
}
rownames(gem_avg) <- paste("GS", unique(sobj_9k@ident), sep = "_")
colnames(gem_avg) <- paste("Cl", ncol, sep = "_")

# Figure 3D right
pheatmap(gem_avg, border="white", color = colorspace::diverge_hsv(50), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, 
         border_color = F, show_rownames = T, scale = "row")

pct_cells <- matrix(NA, nrow = length(unique(sobj_9k@ident)), ncol = length(ncol))
for (c in 1:length(ncol)) {
  cells <- intersect(WhichCells(seurat, ncol(c)), colnames(E18_9k_logtpm))
  for (s in 1:length(celltypes)) {
    avg_exp <- colMeans(E18_9k_logtpm[markers[which(markers$cluster == celltypes[s]), "gene"], ])
    thr <- quantile(avg_exp, 0.95)
    pct_cells[s,c] <- round(length(intersect(names(which(avg_exp > thr)), cells))*100/length(cells), 2)
  }
}

gem_thr <- apply(pct_cells, 2, function(x) sum(x>=100/3))
MNN_ref_counts <- table(gem_thr)/length(unique(seurat@ident))

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scID 
scID_output <- scid_multiclass(target_gem = E18_nuclei_logtpm, markers = markers, likelihood_threshold = 0.99)
scID_labels <- scID_output$labels

# Get ARI and VI (variation of information) with marker based for different thresholds
scID_ARI <- c()
scID_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(biomarker_labels[[i]])
  names(marker_labels) <- colnames(E18_nuclei_logtpm)
  cells <- names(marker_labels)[which(!marker_labels %in% c("ambiguous", "orphan"))]
  scID_ARI[i] <- mclust::adjustedRandIndex(marker_labels[cells], scID_labels[cells])
  scID_VI[i] <- mcclust::vi.dist(marker_labels[cells], scID_labels[cells])
}

# average expression per celltype
markers_filt <- markers[which(markers$gene %in% rownames(E18_nuclei_logtpm)), ]
gem_avg <- matrix(NA, length(unique(sobj_9k@ident)), length(unique(sobj_9k@ident)))
for (i in 1:length(unique(sobj_9k@ident))) {
  cells <- names(scID_labels)[which(scID_labels == i-1)]
  if (length(cells) > 1) {
    avg_exp <- rowMeans(E18_nuclei_logtpm[markers_filt$gene, cells])
  } else if (length(cells) == 1) {
    avg_exp <- E18_nuclei_logtpm[markers_filt$gene, cells]
    names(avg_exp) <- markers_filt$gene
  } else {
    next
  }
  for (j in 1:length(unique(sobj_9k@ident))) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == j-1)]]))
  }
}
rownames(gem_avg) <- paste("GS", unique(sobj_9k@ident), sep = "_")
colnames(gem_avg) <- paste("Cl", unique(sobj_9k@ident), sep = "_")

# Figure 3F right
pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, border_color = F,
         show_rownames = T, scale = "row")

pct_cells <- matrix(NA, length(unique(sobj_9k@ident)), length(unique(sobj_9k@ident)))
for (c in 1:length(unique(sobj_9k@ident))) {
  cells <- names(scID_labels)[which(scID_labels == unique(sobj_9k@ident)[c])]
  for (s in 1:length(celltypes)) {
    avg_exp <- colMeans(E18_nuclei_logtpm[markers_filt[which(markers_filt$cluster == celltypes[s]), "gene"], ])
    thr <- quantile(avg_exp, 0.95)
    pct_cells[s,c] <- round(length(intersect(names(which(avg_exp > thr)), cells))*100/length(cells), 2)
  }
}

NA_cols <- which(apply(pct_cells, 2, function(x) all(is.na(x))))
gem_thr <- apply(pct_cells[, -NA_cols], 2, function(x) sum(x>=100/3))
scID_counts <- table(gem_thr)/(length(unique(sobj_9k@ident))-2)

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scmap
common_genes <- intersect(rownames(E18_9k_logtpm), rownames(E18_nuclei_logtpm))

ann <- data.frame(cell_type1=sobj_9k@ident)#, row.names = rownames(dropseq_groups))

# Create index
ref_sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(E18_9k_logtpm[common_genes, ])), colData = ann)
# use gene names as feature symbols
rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
# remove features with duplicated names
ref_sce <- ref_sce[!duplicated(rownames(ref_sce)), ]
ref_sce <- selectFeatures(ref_sce, suppress_plot = TRUE)#, n_features = 150)
ref_sce <- indexCluster(ref_sce)

# Create sce for testing data (SS2)
sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(E18_nuclei_logtpm[common_genes, ])))
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
names(scmap_labels) <- colnames(E18_nuclei_logtpm)
table(scmap_labels)

# Get ARI with marker based for different thresholds
scmap_ARI <- c()
scmap_VI <- c()
for (i in 1:4) {
  marker_labels <- unlist(biomarker_labels[[i]])
  names(marker_labels) <- colnames(E18_nuclei_logtpm)
  cells <- names(which(marker_labels != "unclassified"))
  scmap_ARI[i] <- mclust::adjustedRandIndex(marker_labels[cells], scmap_labels[cells])
  scmap_VI[i] <- mcclust::vi.dist(marker_labels[cells], scmap_labels[cells])
}

gem_avg <- matrix(NA, length(unique(sobj_9k@ident)), length(unique(sobj_9k@ident)))
for (i in 1:length(unique(sobj_9k@ident))) {
  cells <- names(scmap_labels)[which(scmap_labels == i-1)]
  if (length(cells) > 1) {
    avg_exp <- rowMeans(E18_nuclei_logtpm[markers_filt$gene, cells])
  } else if (length(cells) == 1) {
    avg_exp <- E18_nuclei_logtpm[markers_filt$gene, cells]
    names(avg_exp) <- markers_filt$gene
  } else {
    next
  }
  for (j in 1:length(unique(sobj_9k@ident))) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == j-1)]]))
  }
}
rownames(gem_avg) <- paste("GS", unique(sobj_9k@ident), sep = "_")
colnames(gem_avg) <- paste("Cl", unique(sobj_9k@ident), sep = "_")

# Figure 3F left
pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, 
         border_color = F, show_rownames = T, scale = "row")

pct_cells <- matrix(NA, length(unique(sobj_9k@ident)), length(unique(sobj_9k@ident)))
for (c in 1:length(unique(sobj_9k@ident))) {
  cells <- names(scmap_labels)[which(scmap_labels == unique(sobj_9k@ident)[c])]
  for (s in 1:length(celltypes)) {
    avg_exp <- colMeans(E18_nuclei_logtpm[markers_filt[which(markers_filt$cluster == celltypes[s]), "gene"], ])
    thr <- quantile(avg_exp, 0.95)
    pct_cells[s,c] <- round(length(intersect(names(which(avg_exp > thr)), cells))*100/length(cells), 2)
  }
}

gem_thr <- apply(pct_cells, 2, function(x) sum(x>=100/3))
scmap_counts <- table(gem_thr)/length(unique(sobj_9k@ident))

# ---------------------------------------------------------------
# plot ARI -  Figure 3H
df <- data.frame(value=c(scID_ARI, scmap_ARI), threshold=c(thresholds, thresholds), method=c(rep("scID", 4), rep("scmap", 4)))
df$threshold <- factor(df$threshold, levels = thresholds)

ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=5), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")

# ---------------------------------------------------------------
# plot VI -  Figure 3I
df <- data.frame(value=c(scID_VI, scmap_VI), threshold=c(thresholds, thresholds), method=c(rep("scID", 4), rep("scmap", 4)))
df$threshold <- factor(df$threshold, levels = thresholds)

ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=5), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")

# ---------------------------------------------------------------
# plot number of signatures per cluster (reference data) -  Figure 3E
df <- data.frame(counts = c(Seurat_ref_counts, MNN_ref_counts, CCA_ref_counts),
                 numbers = c(names(Seurat_ref_counts), names(MNN_ref_counts), names(CCA_ref_counts)), 
                 method = c(rep("Seurat", length(Seurat_ref_counts)), 
                            rep("MNN", length(MNN_ref_counts)), 
                            rep("CCA", length(CCA_ref_counts))))

ggplot(df, aes(y=counts, x=method, fill=numbers)) + 
  geom_bar(stat="identity") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=5), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")

# ---------------------------------------------------------------
# plot number of signatures per cluster (target data) - Figure 3G
df <- data.frame(counts = c(Seurat_target_counts, scID_counts, scmap_counts),
                 numbers = c(names(Seurat_target_counts), names(scID_counts), names(scmap_counts)), 
                 method = c(rep("Seurat", length(Seurat_target_counts)), 
                            rep("scID", length(scID_counts)), 
                            rep("scmap", length(scmap_counts))))

ggplot(df, aes(y=counts, x=method, fill=numbers)) + 
  geom_bar(stat="identity") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=5), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")
