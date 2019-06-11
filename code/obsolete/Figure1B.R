# Comparison of accuracy of step 2 (precision-recall) versus step 3(scID score)
# This takes into accound multiclass resolution as well
library(ggplot2)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(rfUtilities)
library(scID)

get_stats <- function(gem, labels, logFC, markers = NULL) {
  
  #library(mclust)
  
  # Remove genes that are zero across all cells
  gem <- gem[which(rowSums(gem) != 0), ]
  
  if (is.null(markers)) {
    markers <- find_markers(gem, labels, logFC, only.pos=FALSE, normalize_reference = FALSE)
  }
  celltypes <- unique(markers$cluster)
  
  # Bin values to 0 and 1 for present (expressed) and absent genes
  sink("aux");
  binned_gem <- apply(gem, 1, function(x) biomod2::BinaryTransformation(x, threshold = quantile(x[which(x>0)], 0.25, na.rm = TRUE)))
  sink(NULL);
  rownames(binned_gem) <- colnames(gem)
  binned_gem <- binned_gem[complete.cases(binned_gem), ]
  
  # Step 2
  sensitivity_S2 <- c()
  specificity_S2 <- c()
  for (i in 1:length(celltypes)) {
    celltype_markers <- markers[which(markers$cluster == celltypes[i]), ]
    positive_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC > 0)]
    negative_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC < 0)]
    
    step2_labels <- choose_training_set(gem, positive_markers, negative_markers)

    predicted_labels <- c(rep("Other", length(step2_labels$out_pop)), rep(as.character(celltypes[i]), length(step2_labels$in_pop)))
    names(predicted_labels) <- c(step2_labels$out_pop, step2_labels$in_pop)
    predicted_labels <- factor(predicted_labels, levels = c(as.character(celltypes[i]), "Other"))
    
    true_labels <- labels
    true_labels[which(true_labels != celltypes[i])] <- "Other"
    true_labels <- factor(true_labels, levels = c(as.character(celltypes[i]), "Other"))
    
    sensitivity_S2 <- c(sensitivity_S2, table(Predictions = predicted_labels, TrueLabels = true_labels[names(predicted_labels)])[1,1]/sum(true_labels==as.character(celltypes[i])))
    specificity_S2 <- c(specificity_S2, table(Predictions = predicted_labels, TrueLabels = true_labels[names(predicted_labels)])[2,2]/sum(true_labels!=as.character(celltypes[i])))
    svMisc::progress(i*100/length(celltypes))
    Sys.sleep(0.01)
  }
  names(sensitivity_S2) <- celltypes
  names(specificity_S2) <- celltypes
  
  rm(binned_gem)
  
  
  stats <- data.frame(sensitivity = sensitivity_S2, 
                      specificity = specificity_S2, 
                      celltype = names(sensitivity_S2)) 
  
  return(stats)
}

# -------------------------------------------------------------------------------------------------------------------------
# Using snRNA-seq mammalian brain
# Read data
gem <- loadfast("~/Desktop/scID_benchmarking_datasets/data/snRNA_mammalian_brain/expression.txt")

cell_info <- read.delim("~/Desktop/scID_benchmarking_datasets/data/snRNA_mammalian_brain/metadata.txt", stringsAsFactors = F, row.names = 1)
labels <- cell_info$Cluster
names(labels) <- rownames(cell_info)
labels <- labels[-1]
rm(cell_info)

markers <- read.delim("~/ownCloud/DocSyncUoE/scID/TEST/Hu_markers_26May_0.3logFC.txt", stringsAsFactors = F)
res_mamBrain <- get_stats(gem, labels, 0.3, markers)
saveRDS(res_mamBrain, file="~/ownCloud/DocSyncUoE/scID/TEST/step2_accuracy_mamBrain.rds")

# -------------------------------------------------------------------------------------------------------------------
# For Shekhar data
gem <- scID::loadfast("~/Google Drive/Data/singlecell/published/mice/shekhar2016_retina/exp_matrix.txt")
rownames(gem) <- toupper(rownames(gem))

# Read labels of DropSeq data
dropseq_groups <- read.delim("~/Google Drive/Data/singlecell/mice/shekhar2016_retina/clust_retinal_bipolar.txt", stringsAsFactors = F, header = T, row.names = 1)
dropseq_groups <- dropseq_groups[-1, ]
doublets <- rownames(dropseq_groups)[which(dropseq_groups$CLUSTER == "Doublets/Contaminants")]
dropseq_groups <- dropseq_groups[setdiff(rownames(dropseq_groups), doublets), ]
gem <- gem[, setdiff(colnames(gem), doublets)]
labels <- dropseq_groups[colnames(gem), "CLUSTER"]
names(labels) <- rownames(dropseq_groups)

res_Shekhar <- get_stats(gem, labels, 0.7)
saveRDS(res_Shekhar, file="~/ownCloud/DocSyncUoE/scID/TEST/step2_accuracy_Shekhar.rds")

# -------------------------------------------------------------------------------------------------------------------------
# Using Tirosh non-malignant cells
gem <- scID::loadfast("~/Desktop/scID_benchmarking_datasets/data/melanoma_immunotherapy_res/tumors_tpm.txt")
cell.info <- read.delim("~/Desktop/scID_benchmarking_datasets/data/melanoma_immunotherapy_res/tumors.nonmal_tsne_anno.txt", stringsAsFactors = F)
cell.info <- cell.info[-1, ]
labels <- cell.info$cell.type
names(labels) <- cell.info$NAME
# Remove generic T.cells
labels <- labels[-which(labels == "T.cell")]

common_cells <- intersect(names(labels), colnames(gem))
labels <- labels[common_cells]
gem <- gem[, common_cells]

rm("common_cells", "cell.info")

res_Tirosh <- get_stats(gem, labels, 0.5)
saveRDS(res_Tirosh, file="~/ownCloud/DocSyncUoE/scID/TEST/step2_accuracy_Tirosh.rds")

# # -------------------------------------------------------------------------------------------------------------------------
# # Using intestinal stem cells (Biton)
# # Read data
# gem <- read.delim("~/Desktop/scID_benchmarking_datasets/data/intestinal_stem_cells/FigS4D_Log2TPM.txt", stringsAsFactors = F, row.names = 1)
# 
# labels <- unlist(lapply(colnames(gem), function(x) tail(strsplit(x, "_")[[1]], 1)))
# names(labels) <- colnames(gem)
# rnames <- toupper(rownames(gem))
# gem <- gem[-which(duplicated(rnames)), ]
# 
# res_Biton <- get_stats(gem, labels, 0.05)
# saveRDS(res_Biton, file="~/ownCloud/DocSyncUoE/scID/TEST/step2_vs_step3_Biton.rds")

# -------------------------------------------------------------------------------------------------------------------------
# Using airway epithelial (Montoro)
# Read data
gem <- read.delim("~/Desktop/scID_benchmarking_datasets/data/airway_epithelium/trachea_10x_log2TPM.txt", stringsAsFactors = F, row.names = 1)
rnames <- toupper(rownames(gem))
gem <- gem[-which(duplicated(rnames)), ]
# Remove genes that are 0 across all cells
gem <- gem[-which(rowSums(gem) == 0), ]

cell_info <- read.delim("~/Desktop/scID_benchmarking_datasets/data/airway_epithelium/trachea_10x_metadata.txt", stringsAsFactors = F, row.names = 1)
cell_info <- cell_info[-1, ]
labels <- cell_info$cluster
names(labels) <- rownames(cell_info)

rm(rnames, cell_info)

res_Montoro <- get_stats(gem, labels, 0.5)
saveRDS(res_Montoro, file="~/ownCloud/DocSyncUoE/scID/TEST/step2_accuracy_Montoro.rds")

# # -----------------------------------------------------------
# # Plot results
res_Shekhar <- readRDS("~/ownCloud/DocSyncUoE/scID/TEST/step2_accuracy_Shekhar.rds")
res_Tirosh <- readRDS("~/ownCloud/DocSyncUoE/scID/TEST/step2_accuracy_Tirosh.rds")
res_Montoro <- readRDS("~/ownCloud/DocSyncUoE/scID/TEST/step2_accuracy_Montoro.rds")
res_mamBrain <- readRDS("~/ownCloud/DocSyncUoE/scID/TEST/step2_accuracy_mamBrain.rds")
df <- data.frame(sensitivity = c(res_Shekhar$sensitivity, res_Tirosh$sensitivity, res_Montoro$sensitivity, res_mamBrain$sensitivity),
                 specificity = c(res_Shekhar$specificity, res_Tirosh$specificity, res_Montoro$specificity, res_mamBrain$specificity),
                 dataset = c(rep("Shekhar", nrow(res_Shekhar)), rep("Tirosh", nrow(res_Tirosh)),
                             rep("Montoro", nrow(res_Montoro)), rep("mammalian Brain", nrow(res_mamBrain))))

df$FPR <- 1 - df$specificity
df.m <- reshape2::melt(df[, -2])
# True Positive Rate (sensitivity)
p_sensitivity <- ggplot(df, aes(x=dataset, y=sensitivity)) + geom_boxplot(notch=TRUE, width = 0.5, outlier.size = 0.05, lwd = 0.5) +
  labs(title="", x="", y = "TPR") +
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_jco() + theme(legend.position="none")

# False Positive Rate
p_fpr <- ggplot(df, aes(x=dataset, y=FPR)) + geom_boxplot(notch=TRUE, width = 0.5, outlier.size = 0.05, lwd = 0.5) +
  labs(title="", x="", y = "FPR") +
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_jco() + theme(legend.position="none")

#pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/New plots/precision_recall_performance.pdf", width = 8, height = 4)
pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure1/step2_accuracy.pdf", width = 3, height = 2)
gridExtra::grid.arrange(p_sensitivity, p_fpr, ncol=2)
dev.off()




