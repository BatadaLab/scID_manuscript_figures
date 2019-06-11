# Accuracy of Stage 2 (precision-recall-like)

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(rfUtilities)
library(scID)

get_stats <- function(gem, labels, logFC, markers = NULL) {
  
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
  TPR <- c()
  FPR <- c()
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
    
    TPR <- c(TPR, table(Predictions = predicted_labels, TrueLabels = true_labels[names(predicted_labels)])[1,1]/sum(true_labels==as.character(celltypes[i])))
    FPR <- c(FPR, 1-table(Predictions = predicted_labels, TrueLabels = true_labels[names(predicted_labels)])[2,2]/sum(true_labels!=as.character(celltypes[i])))
    svMisc::progress(i*100/length(celltypes))
    Sys.sleep(0.01)
  }
  names(TPR) <- celltypes
  names(FPR) <- celltypes
  
  rm(binned_gem)
  
  stats <- data.frame(TPR = TPR, specificity = FPR, celltype = names(TPR)) 
  
  return(stats)
}

# ------------------------------------------------------------------------------------------------------------------------
# Using Hu 2017
# ------------------------------------------------------------------------------------------------------------------------
# Read data
gem <- readRDS("~/scID_manuscript_figures/data/Figure1/Hu2017_gem.rds")
labels <- readRDS("~/scID_manuscript_figures/data/Figure1/Hu2017_labels.rds")

res_mamBrain <- get_stats(gem, labels, 0.3)

# ------------------------------------------------------------------------------------------------------------------------
# Using Shekhar 2016 data
# ------------------------------------------------------------------------------------------------------------------------
gem <- readRDS("~/scID_manuscript_figures/data/Figure2/Reference_gem.rds")
labels <- readRDS("~/scID_manuscript_figures/data/Figure2/Reference_clusters.rds")

res_Shekhar <- get_stats(gem, labels, 0.7)

# -------------------------------------------------------------------------------------------------------------------------
# Using Tirosh 2016
# ------------------------------------------------------------------------------------------------------------------------
gem <- readRDS("~/scID_manuscript_figures/data/Figure1/Tirosh2016_gem.rds")
labels <- readRDS("~/scID_manuscript_figures/data/Figure1/Tirosh2016_labels.rds")

res_Tirosh <- get_stats(gem, labels, 0.5)

# -------------------------------------------------------------------------------------------------------------------------
# Using Montoro 2018
# ------------------------------------------------------------------------------------------------------------------------
# Read data
gem <- readRDS("~/scID_manuscript_figures/data/Figure1/Montoro2018_gem.rds")
labels <- readRDS("~/scID_manuscript_figures/data/Figure1/Montoro2018_labels.rds")

res_Montoro <- get_stats(gem, labels, 0.5)

# ------------------------------------------------------------------------------------------------------------------------
# Plot results
df <- data.frame(TPR = c(res_Shekhar$TPR, res_Tirosh$TPR, res_Montoro$TPR, res_mamBrain$TPR),
                 FPR = c(res_Shekhar$FPR, res_Tirosh$FPR, res_Montoro$FPR, res_mamBrain$FPR),
                 dataset = c(rep("Shekhar", nrow(res_Shekhar)), rep("Tirosh", nrow(res_Tirosh)),
                             rep("Montoro", nrow(res_Montoro)), rep("mammalian Brain", nrow(res_mamBrain))))

df.m <- reshape2::melt(df)

ggplot(df.m, aes(x=dataset, y=value, fill=variable)) + geom_boxplot(notch=TRUE, width = 0.5, outlier.size = 0.05, lwd = 0.5) +
  labs(title="", x="", y = "TPR") +
  theme(legend.position="top", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
