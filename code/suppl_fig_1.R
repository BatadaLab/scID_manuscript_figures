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
# Compare weights and accuracy of scID using the Smart-Seq2 dataset as reference and the celseq datasets as targets
common_celltypes <- intersect(names(celseq_res$estimated_weights), names(celseq2_res$estimated_weights))
correlations <- c()
p_values <- c() 
pdf("~/ownCloud/DocSyncUoE/Thesis/thesis_full/chapter2/chapter/figs/pancreatic_correlations.pdf", width = 8, height = 6)
par(mfrow=c(3,4))
for (celltype in common_celltypes) {
  
  w1 <- celseq_res$estimated_weights[[celltype]]
  w2 <- celseq2_res$estimated_weights[[celltype]]
  
  common_genes <- intersect(names(w1), names(w2))
  
  c <- cor.test(w1[common_genes], w2[common_genes], method = "spearman")
  correlations <- c(correlations, c$estimate)
  p_values <- c(p_values, c$p.value)
  plot(w1[common_genes], w2[common_genes], pch=16, xlab = "Weights from CEL-Seq", ylab = "Weights from CEL-Seq2", 
       main = sprintf("%s - cor=%s \np-val=%s", celltype, round(c$estimate,2), c$p.value))
}
dev.off()
names(correlations) <- common_celltypes


# -----------------------------------------
# Show accuracy for each dataset
df <- data.frame(ARI = c(scID_celseq_ARI, scID_celseq2_ARI), dataset = c("CEL-Seq", "CEL-Seq2"))

ggplot(df, aes(y=ARI, x=dataset)) + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco() + 
  theme(legend.position="top", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "", x="", y="Adjusted Rand Index")
