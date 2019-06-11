library(Seurat) # version 3
library(dplyr)
library(pheatmap)
library(scran)
library(scmap)
library(dplyr)
library(ggsci)
library(scID)
library(ggplot2)
library(ggpubr)

runCastle <- function(target_gem, source_gem, source_identities) {
  library(scater)  
  library(xgboost) 
  library(igraph)  
  
  BREAKS=c(-1, 0, 1, 6, Inf)
  nFeatures = 200
  
  source = SingleCellExperiment(assays = list(logcounts = as.matrix(source_gem)), 
                                colData = data.frame(cell_type1 = source_identities))
  
  target = SingleCellExperiment(assays = list(logcounts = as.matrix(target_gem)))
  
  ds1 = t(exprs(source)) 
  ds2 = t(exprs(target)) 
  sourceCellTypes = as.factor(colData(source)[,"cell_type1"])
  
  # 2. Unify sets, excluding low expressed genes
  source_n_cells_counts = apply(exprs(source), 1, function(x) { sum(x > 0) } )
  target_n_cells_counts = apply(exprs(target), 1, function(x) { sum(x > 0) } )
  common_genes = intersect(rownames(source)[source_n_cells_counts>10], rownames(target)[target_n_cells_counts>10])
  
  remove(source_n_cells_counts, target_n_cells_counts)
  ds1 = ds1[, colnames(ds1) %in% common_genes]
  ds2 = ds2[, colnames(ds2) %in% common_genes]
  ds = rbind(ds1[,common_genes], ds2[,common_genes])
  isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
  remove(ds1, ds2)
  
  # 3. Highest mean in both source and target
  topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
  
  # 4. Highest mutual information in source
  topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),sourceCellTypes,method = "nmi") }), decreasing = T))
  
  # 5. Top n genes that appear in both mi and avg
  selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )
  
  # 6. remove correlated features
  tmp = cor(ds[,selectedFeatures], method = "pearson")
  tmp[!lower.tri(tmp)] = 0
  selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
  remove(tmp)
  
  # 7,8. Convert data from continous to binned dummy vars
  # break datasets to bins
  dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
  # use only bins with more than one value
  nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
  # convert to dummy vars
  ds = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
  remove(dsBins, nUniq)
  
  # 9. Classify
  train = runif(nrow(ds[isSource,])) < 0.8
  # slightly different setup for multiclass and binary classification
  if (length(unique(sourceCellTypes)) > 2) {
    xg=xgboost(data=ds[isSource,][train, ] , 
               label=as.numeric(sourceCellTypes[train])-1,
               objective="multi:softmax", num_class=length(unique(sourceCellTypes)),
               eta=0.7 , nthread=5, nround=20, verbose=0,
               gamma=0.001, max_depth=5, min_child_weight=10)
  } else {
    xg=xgboost(data=ds[isSource,][train, ] , 
               label=as.numeric(sourceCellTypes[train])-1,
               eta=0.7 , nthread=5, nround=20, verbose=0,
               gamma=0.001, max_depth=5, min_child_weight=10)
  }
  
  # 10. Predict
  predictedClasses = predict(xg, ds[!isSource, ])
  names(predictedClasses) <- colnames(target_gem)
  
  celltype_dict <- levels(sourceCellTypes)
  predictedClasses <- unlist(lapply(predictedClasses, function(x) celltype_dict[x+1]))
  
  return(predictedClasses)
}

# ---------------------------------------------------------------------------
# Read data
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
reference_clusters <- readRDS("~/scID_manuscript_figures/data/Figure3/reference_clusters.rds")
sobj_9k <- CreateSeuratObject(counts = E18_9k)
sobj_9k <- NormalizeData(sobj_9k)
#sobj_9k <- FindVariableFeatures(sobj_9k)
sobj_9k <- ScaleData(sobj_9k)
#sobj_9k <- RunPCA(sobj_9k,  do.print = FALSE)
#sobj_9k <- FindNeighbors(sobj_9k)
#sobj_9k <- FindClusters(sobj_9k)
Idents(sobj_9k) <- factor(reference_clusters)
#sobj_9k <- RunTSNE(sobj_9k)

# Figure 3A left
#TSNEPlot(sobj_9k, do.label = T, pt.size = 0.05, no.axes = T, no.legend = T)

# Find markers or use pre-computed list
#markers <- FindAllMarkers(sobj_9k, only.pos = TRUE, test.use = "MAST", logfc.threshold = 0.5, do.print = FALSE)
# Changed the logFC threshold 28 May 2019
markers <- FindAllMarkers(sobj_9k, only.pos = TRUE, test.use = "MAST", logfc.threshold = 0.3, do.print = FALSE)
#markers <- readRDS("~/scID_manuscript_figures/data/Figure3/markers.rds")

celltypes <- unique(Idents(sobj_9k))

gem_avg_9k <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- WhichCells(sobj_9k, idents = celltypes[i])
  if (length(celltypes) >= 1) {
    avg_exp <- rowMeans(E18_9k_cpm[markers$gene, cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg_9k[j,i] <- mean(na.omit(avg_exp[markers$gene[which(markers$cluster == celltypes[j])]]))
  }
}
rownames(gem_avg_9k) <- paste("GS", celltypes, sep = "_")
colnames(gem_avg_9k) <- paste("Cl", celltypes, sep = "_")

# Figure 3B left
pheatmap(gem_avg_9k, border="white", color = colorspace::diverge_hsv(50), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, 
         border_color = F, show_rownames = T, scale = "row")

pct_cells <- matrix(NA, length(celltypes), length(celltypes))
for (c in 1:length(celltypes)) {
  cells <- WhichCells(sobj_9k, idents = celltypes[c])
  for (s in 1:length(celltypes)) {
    avg_exp <- colMeans(E18_9k_cpm[markers[which(markers$cluster == celltypes[s]), "gene"], ])
    thr <- quantile(avg_exp, 0.95)
    pct_cells[s,c] <- round(length(intersect(names(which(avg_exp > thr)), cells))*100/length(cells), 2)
  }
}

gem_thr <- apply(pct_cells, 2, function(x) sum(x>=100/3))
Seurat_ref_counts <- table(gem_thr)/length(celltypes)

# -----------------------------------------------------------------------------------------------------------------------
# Marker based identification
markers_filt <- markers[-which(duplicated(markers$gene)), ]
top_markers <- as.data.frame(markers_filt %>% group_by(cluster) %>% top_n(1, avg_logFC))
top_markers <- top_markers[which(top_markers$gene %in% rownames(E18_nuclei_cpm)), ]
rownames(top_markers) <- top_markers$gene

# Normalize data by 99th percentile
nuclei_norm <- t(apply(E18_nuclei_cpm[intersect(top_markers$gene, rownames(E18_nuclei_cpm)), ], 1, function(x) normalize_gene(x)))

thresholds <- c(0.1, 0.25, 0.5, 0.75)
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
nuclei_clustered <- CreateSeuratObject(counts = E18_nuclei)
nuclei_clustered <- NormalizeData(nuclei_clustered)
nuclei_clustered <- FindVariableFeatures(nuclei_clustered)
nuclei_clustered <- ScaleData(nuclei_clustered)
nuclei_clustered <- RunPCA(nuclei_clustered)
nuclei_clustered <- FindNeighbors(nuclei_clustered)
nuclei_clustered <- FindClusters(nuclei_clustered)
nuclei_clustered <- RunTSNE(nuclei_clustered)

# Figure 3A right
TSNEPlot(nuclei_clustered, do.label = T, pt.size = 0.05, no.axes = T, no.legend = T)

markers_filt <- markers[which(markers$gene %in% rownames(E18_nuclei_cpm)), ]
nucClust <- unique(Idents(nuclei_clustered))
nuclei_avg <- matrix(NA, length(celltypes), ncol)
for (i in 1:length(nucClust)) {
  cells <- WhichCells(nuclei_clustered, idents = nucClust[i])
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(E18_nuclei_cpm[markers_filt$gene, cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    nuclei_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == celltypes[j])]]))
  }
}
rownames(nuclei_avg) <- paste("GS", celltypes, sep = "_")
colnames(nuclei_avg) <- paste("Cl", nucClust, sep = "_")

# Figure 3B right
pheatmap(nuclei_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, border_color = F,
         show_rownames = T, scale = "row")

pct_cells <- matrix(NA, length(celltypes), length(nucClust))
for (c in 1:length(nucClust)) {
  cells <- WhichCells(nuclei_clustered, idents = nucClust[c])
  for (s in 1:length(celltypes)) {
    avg_exp <- colMeans(E18_nuclei_cpm[markers_filt[which(markers_filt$cluster == celltypes[s]), "gene"], ])
    thr <- quantile(avg_exp, 0.95)
    pct_cells[s,c] <- round(length(intersect(names(which(avg_exp > thr)), cells))*100/length(cells), 2)
  }
}

gem_thr <- apply(pct_cells, 2, function(x) sum(x>=100/3))
Seurat_target_counts <- table(gem_thr)/length(nucClust)

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
rownames(gem_avg_ds) <- paste("GS", unique(sobj_9k@ident), sep="_")
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
scID_output <- scid_multiclass(target_gem = E18_nuclei_cpm, markers = markers, estimate_weights_from_target = TRUE)
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

# average expression per celltype
markers_filt <- markers[which(markers$gene %in% rownames(E18_nuclei_cpm)), ]
gem_avg <- matrix(NA, length(celltypes), length(celltypes))
for (i in 1:length(celltypes)) {
  cells <- names(scID_labels)[which(scID_labels == celltypes[i])]
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(E18_nuclei_cpm[markers_filt$gene, cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == celltypes[j])]]))
  }
}
rownames(gem_avg) <- paste("GS", celltypes, sep = "_")
colnames(gem_avg) <- paste("Cl", celltypes, sep = "_")

# Figure 3F right
pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, border_color = F,
         show_rownames = T, scale = "row")

gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
scID_counts <- table(colSums(gem_thr))/length(unique(celltypes))
scID_counts <- c(1-sum(scID_counts), scID_counts)
names(scID_counts) <- c(0, names(scID_counts[-1]))

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scmap
common_genes <- intersect(rownames(E18_9k_cpm), rownames(E18_nuclei_cpm))

ann <- data.frame(cell_type1=Idents(sobj_9k))

# Create index
ref_sce <- SingleCellExperiment(assays = list(logcounts = log(E18_9k_cpm[common_genes, ] + 1)), colData = ann)
# use gene names as feature symbols
rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
# remove features with duplicated names
ref_sce <- ref_sce[!duplicated(rownames(ref_sce)), ]
ref_sce <- selectFeatures(ref_sce, suppress_plot = TRUE)#, n_features = 150)
ref_sce <- indexCluster(ref_sce)

# Create sce for testing data (SS2)
sce <- SingleCellExperiment(assays = list(logcounts = log(E18_nuclei_cpm[common_genes, ] + 1)))
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
names(scmap_labels) <- colnames(E18_nuclei_cpm)
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
  cells <- names(scmap_labels)[which(scmap_labels == celltypes[i])]
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(E18_nuclei_cpm[markers_filt$gene, cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == celltypes[j])]]))
  }
}
rownames(gem_avg) <- paste("GS", celltypes, sep = "_")
colnames(gem_avg) <- paste("Cl", celltypes, sep = "_")

# Figure 3F left
pheatmap(gem_avg, border="white",color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
         cluster_rows = F, cluster_cols = F, show_colnames = T, fontsize = 5, 
         border_color = F, show_rownames = T, scale = "row")

gem_avg <- gem_avg[, complete.cases(t(gem_avg))]
gem_avg_norm <- apply(gem_avg, 1, function(x) (x-min(x)) /(max(x) - min(x)))
gem_thr <- apply(gem_avg_norm, 1, function(x) {ifelse(x >= quantile(x, 0.95), 1, 0)})
scmap_counts <- table(colSums(gem_thr))/length(unique(celltypes))
scmap_counts <- c(1-sum(scmap_counts), scmap_counts)
names(scmap_counts) <- c(0, names(scmap_counts[-1]))

# ---------------------------------------------------------------
# find equivalent cells with CaSTLe
castle_labels <- runCastle(E18_nuclei_cpm, E18_9k_cpm, Idents(sobj_9k))

# Get ARI and VI (variation of information) with marker based for different thresholds
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
  cells <- na.omit(names(castle_labels)[which(castle_labels == celltypes[i])])
  if (length(cells) >= 1) {
    avg_exp <- rowMeans(E18_nuclei_cpm[toupper(markers_filt$gene), cells, drop = FALSE])
  } else {
    next
  }
  for (j in 1:length(celltypes)) {
    gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(markers_filt$gene[which(markers$cluster == celltypes[j])])]))
  }
}
rownames(gem_avg) <- celltypes
colnames(gem_avg) <- celltypes

# Figure 2F right
pheatmap(gem_avg, border="white", color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
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
# plot ARI -  Figure 3H
df <- data.frame(value=c(scID_ARI, scmap_ARI, castle_ARI), 
                 threshold=c(thresholds, thresholds, thresholds), 
                 method=c(rep("scID", length(thresholds)), rep("scmap", length(thresholds)), rep("CaSTLe", length(thresholds))))
df$threshold <- factor(df$threshold, levels = thresholds)

ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")

# ---------------------------------------------------------------
# plot VI -  Figure 3I
df <- data.frame(value=c(scID_VI, scmap_VI, castle_VI), 
                 threshold=c(thresholds, thresholds, thresholds), 
                 method=c(rep("scID", length(thresholds)), rep("scmap", length(thresholds)), rep("CaSTLe", length(thresholds))))
df$threshold <- factor(df$threshold, levels = thresholds)

ggplot(df, aes(y=value, x=threshold, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() +
  theme(legend.position="top", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")

# ---------------------------------------------------------------
# plot number of signatures per cluster (reference data) -  Figure 3E
df <- data.frame(counts = c(Seurat_ref_counts, MNN_ref_counts, CCA_ref_counts),
                 numbers = c(names(Seurat_ref_counts), names(MNN_ref_counts), names(CCA_ref_counts)), 
                 method = c(rep("Seurat", length(Seurat_ref_counts)), 
                            rep("MNN", length(MNN_ref_counts)), 
                            rep("CCA", length(CCA_ref_counts))))

ggplot(df, aes(y=counts, x=method, fill=numbers)) + 
  geom_bar(stat="identity") + scale_fill_grey()+
  theme(legend.position="top", text = element_text(size=5), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")

# ---------------------------------------------------------------
# plot number of signatures per cluster (target data) - Figure 3G
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
