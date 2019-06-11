# Comparison of accuracy of step 2 (precision-recall) versus step 3(scID score)

library(ggplot2)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(rfUtilities)
library(scID)

# -------------------------------------------------------------------------------------------------------------------------
# Using Shekhar data
gem <- scID::loadfast("~/Google Drive/Data/singlecell/published/mice/shekhar2016_retina/exp_matrix.txt")
rownames(gem) <- toupper(rownames(gem))

# Read labels of DropSeq data
dropseq_groups <- read.delim("~/Google Drive/Data/singlecell/published/mice/shekhar2016_retina/clust_retinal_bipolar.txt", stringsAsFactors = F, header = T, row.names = 1)
dropseq_groups <- dropseq_groups[-1, ]
doublets <- rownames(dropseq_groups)[which(dropseq_groups$CLUSTER == "Doublets/Contaminants")]
dropseq_groups <- dropseq_groups[setdiff(rownames(dropseq_groups), doublets), ]
gem <- gem[, setdiff(colnames(gem), doublets)]
labels <- dropseq_groups[colnames(gem), "CLUSTER"]
names(labels) <- rownames(dropseq_groups)

rm("dropseq_groups", "doublets")

rownames(gem) <- make.names(toupper(rownames(gem)), unique=TRUE)
gem <- gem[which(rowSums(gem) != 0), ]
markers <- find_markers(gem, labels, 0.5, only.pos=FALSE, normalize_reference=FALSE)
celltypes <- unique(markers$cluster)

# binarize gem
sink("aux");
binned_gem <- apply(gem, 1, function(x) biomod2::BinaryTransformation(x, threshold = quantile(x[which(x>0)], 0.25, na.rm = TRUE)))
sink(NULL);

TPR_step2 <- c()
FPR_step2 <- c()
TPR_step3 <- c()
FPR_step3 <- c()

gem_norm <- t(apply(gem[unique(markers$gene), ], 1, function(x) normalize_gene(x)))
gem_norm <- gem_norm[complete.cases(gem_norm), ]

for (i in 1:length(celltypes)) {
  celltype_markers <- markers[which(markers$cluster == celltypes[i]), ]
  positive_markers <- intersect(celltype_markers$gene[which(celltype_markers$avg_logFC > 0)], colnames(binned_gem))
  negative_markers <- intersect(celltype_markers$gene[which(celltype_markers$avg_logFC < 0)], colnames(binned_gem))
  
  # Step 2
  training_groups <- choose_training_set(gem, positive_markers, negative_markers)
  
  TPR_step2 <- c(TPR_step2, length(intersect(training_groups$in_pop, names(which(labels == celltypes[i]))))/sum(labels == celltypes[i]))
  FPR_step2 <- c(FPR_step2, length(setdiff(training_groups$in_pop, names(which(labels == celltypes[i]))))/sum(labels != celltypes[i]))
  
  # Step 3
  signature_genes <- c(positive_markers, negative_markers)
  gene.weights <- scID_weight(gem_norm[signature_genes, , drop=FALSE], training_groups$in_pop, training_groups$out_pop)
  gene.weights[is.infinite(gene.weights)] <- 0
  weighted_gem <- gene.weights[signature_genes] * gem_norm[signature_genes, ,drop=FALSE]
  if (!all(weighted_gem == 0)) {
    score <- colSums(weighted_gem)/sqrt(sum(gene.weights^2))
    step3_cells <- final_populations(score) 
    
    TPR_step3 <- c(TPR_step3, length(intersect(step3_cells, names(which(labels == celltypes[i]))))/sum(labels == celltypes[i]))
    FPR_step3 <- c(FPR_step3, length(setdiff(step3_cells, names(which(labels == celltypes[i]))))/sum(labels != celltypes[i]))
  } else {
    TPR_step3 <- c(TPR_step3, 0)
    FPR_step3 <- c(FPR_step3, 0)
  }

  
  svMisc::progress(i*100/length(celltypes))
  Sys.sleep(0.01)
  if (i==length(celltypes)) cat("Done!")
}

df <- data.frame(TPR = c(TPR_step2, TPR_step3),
                 FPR = c(FPR_step2, FPR_step3), 
                 step = c(rep("Step2", length(TPR_step2)), rep("Step3", length(TPR_step3))))

df.m <- reshape2::melt(df)

# True Positive Rate (sensitivity)
p_TPR <- ggplot(df, aes(x=step, y=TPR, fill=step)) + geom_boxplot(notch=TRUE, width = 0.3, outlier.size = 0.05, lwd = 0.5) +
  labs(title="", x="", y = "TPR") + 
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_compare_means(paired = TRUE) +
  scale_fill_jco() + theme(legend.position="none")

# specificity
p_FPR <- ggplot(df, aes(x=step, y=FPR, fill=step)) + geom_boxplot(notch=TRUE, width = 0.3, outlier.size = 0.05, lwd = 0.5) +
  labs(title="", x="", y = "FPR") + stat_compare_means(paired = TRUE) +
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_jco() + theme(legend.position="none")

#pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/New plots/step2_step3_comparison.pdf", width = 8, height = 4)
pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure1/step2_vs_step3.pdf", width = 3, height = 2)
gridExtra::grid.arrange(p_TPR, p_FPR, ncol=2)
dev.off()
