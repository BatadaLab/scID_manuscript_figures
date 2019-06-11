# Using pancreatic data from CCA
library(Seurat)
library(cowplot)
library(ggplot2)
library(scID)
library(ggsci)
source("~/scID_manuscript_figures/code/run_other_methods.R")

# ------------------------------------------------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------------------------------------------------
pancreas.data <- readRDS(file = "~/Desktop/scID_benchmarking_datasets/data/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata <- readRDS(file = "~/Desktop/scID_benchmarking_datasets/data/pancreas_v3_files/pancreas_metadata.rds")
true_labels <- metadata$celltype
names(true_labels) <- rownames(metadata)  

pancreas <- CreateSeuratObject(pancreas.data, meta.data = metadata)
pancreas.list <- SplitObject(pancreas, split.by = "tech")

celseq_gem <- as.data.frame(GetAssayData(pancreas.list[["celseq"]], slot = "counts"))
celseq2_gem <- as.data.frame(GetAssayData(pancreas.list[["celseq2"]], slot = "counts"))
ss2_gem <- as.data.frame(GetAssayData(pancreas.list[["smartseq2"]], slot = "counts"))
ss2_gem_cpm <- counts_to_cpm(ss2_gem)

reference_clusters <- metadata[colnames(ss2_gem), "celltype"]
names(reference_clusters) <- colnames(ss2_gem)

# ------------------------------------------------------------------------------------------------------------------------
# scID
# ------------------------------------------------------------------------------------------------------------------------
celseq_res <- scid_multiclass(target_gem = celseq_gem, reference_gem = ss2_gem, reference_clusters = reference_clusters, 
                              logFC = 0.5, estimate_weights_from_target = TRUE, only_pos = TRUE, normalize_reference = TRUE)
scID_celseq_labels <- celseq_res$labels

celseq2_res <- scid_multiclass(target_gem = celseq2_gem, markers = celseq_res$markers, estimate_weights_from_target = TRUE)
scID_celseq2_labels <- celseq2_res$labels

scID_celseq_ARI <- mclust::adjustedRandIndex(scID_celseq_labels[names(which(scID_celseq_labels != "unassigned"))], true_labels[names(which(scID_celseq_labels != "unassigned"))])
scID_celseq2_ARI <- mclust::adjustedRandIndex(scID_celseq2_labels[names(which(scID_celseq2_labels != "unassigned"))], true_labels[names(which(scID_celseq2_labels != "unassigned"))])

# ------------------------------------------------------------------------------------------------------------------------
# scmap
# ------------------------------------------------------------------------------------------------------------------------
scmap_celseq_labels <- run_scmap(reference_gem = log2(ss2_gem_cpm+1), reference_labels = true_labels[colnames(ss2_gem)], target_gem = celseq_gem)
scmap_celseq_ARI <- mclust::adjustedRandIndex(scmap_celseq_labels[names(which(scmap_celseq_labels != "unassigned"))], true_labels[names(which(scmap_celseq_labels != "unassigned"))])

scmap_celseq2_labels <- run_scmap(reference_gem = log2(ss2_gem_cpm+1), reference_labels = true_labels[colnames(ss2_gem)], target_gem = celseq2_gem)
scmap_celseq2_ARI <- mclust::adjustedRandIndex(scmap_celseq2_labels[names(which(scmap_celseq2_labels != "unassigned"))], true_labels[names(which(scmap_celseq2_labels != "unassigned"))])

# ------------------------------------------------------------------------------------------------------------------------
# CaSTLE
# ------------------------------------------------------------------------------------------------------------------------
castle_celseq_labels <- runCastle(target_gem = celseq_gem, source_gem = log2(ss2_gem_cpm+1), source_identities = true_labels[colnames(ss2_gem)])
castle_celseq_ARI <- mclust::adjustedRandIndex(castle_celseq_labels, true_labels[names(castle_celseq_labels)])

castle_celseq2_labels <- runCastle(target_gem = celseq2_gem, source_gem = log2(ss2_gem_cpm+1), source_identities = true_labels[colnames(ss2_gem)])
castle_celseq2_ARI <- mclust::adjustedRandIndex(castle_celseq2_labels, true_labels[names(castle_celseq2_labels)])

# ------------------------------------------------------------------------------------------------------------------------
# CCA
# ------------------------------------------------------------------------------------------------------------------------
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
                                             verbose = FALSE)
}

reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

# switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- FindNeighbors(pancreas.integrated)
pancreas.integrated <- FindClusters(pancreas.integrated)

Seurat_labels <- Idents(pancreas.integrated)
table(Seurat_labels)

Seurat_celseq_ARI <- mclust::adjustedRandIndex(Seurat_labels[colnames(celseq_gem)], true_labels[colnames(celseq_gem)])
Seurat_celseq2_ARI <- mclust::adjustedRandIndex(Seurat_labels[colnames(celseq2_gem)], true_labels[colnames(celseq2_gem)])

# ------------------------------------------------------------------------------------------------------------------------
# MNN
# ------------------------------------------------------------------------------------------------------------------------
MNN_celseq_labels <- runMNN(celseq_gem, log2(ss2_gem_cpm+1))
labels <- MNN_celseq_labels$target
MNN_celseq_ARI <- mclust::adjustedRandIndex(labels, true_labels[names(labels)])

MNN_celseq2_labels <- runMNN(celseq2_gem, log2(ss2_gem_cpm+1))
labels <- MNN_celseq2_labels$target
MNN_celseq2_ARI <- mclust::adjustedRandIndex(labels, true_labels[names(labels)])

# ------------------------------------------------------------------------------------------------------------------------
# Plot results
# ------------------------------------------------------------------------------------------------------------------------
df <- data.frame(ARI = c(scID_celseq_ARI, scID_celseq2_ARI, Seurat_celseq_ARI, Seurat_celseq2_ARI, scmap_celseq_ARI, scmap_celseq2_ARI, 
                         castle_celseq_ARI, castle_celseq2_ARI, MNN_celseq_ARI, MNN_celseq2_ARI), 
                 method = c("scID", "scID", "CCA", "CCA", "scmap", "scmap", "CaSTLe", "CaSTLe", "MNN", "MNN"), 
                 dataset = c("CEL-Seq", "CEL-Seq2", "CEL-Seq", "CEL-Seq2", "CEL-Seq", "CEL-Seq2", "CEL-Seq", "CEL-Seq2", "CEL-Seq", "CEL-Seq2"))

ggplot(df, aes(y=ARI, x=method, fill=method)) + facet_wrap(~dataset) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() + 
  theme(legend.position="top", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "", x="", y="Adjusted Rand Index")
