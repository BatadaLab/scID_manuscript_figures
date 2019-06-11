library(Seurat)
library(dplyr)
library(pheatmap)
library(scran)
library(scmap)
library(ggsci)
library(cowplot)
library(scID)
source("~/scID_manuscript_figures/code/run_other_methods.R")

# ------------------------------------------------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------------------------------------------------
ctrl.data <- read.table(file = "~/scID_manuscript_figures/data/Kang_data/immune_control_expression_matrix.txt.gz", sep = "\t")
ctrl.data.cpm <- counts_to_cpm(ctrl.data)

stim.data <- read.table(file = "~/scID_manuscript_figures/data/Kang_data/immune_stimulated_expression_matrix.txt.gz", sep = "\t")
stim.data.cpm <- counts_to_cpm(stim.data)

# ------------------------------------------------------------------------------------------------------------------------
# Reference labels from CCA
# ------------------------------------------------------------------------------------------------------------------------
immune.combined <- readRDS("~/scID_manuscript_figures/data/Kang_data/immuned_combined.rds")

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent cells with scID
# ------------------------------------------------------------------------------------------------------------------------
scID_res <- scid_multiclass(target_gem = stim.data.cpm, reference_gem = ctrl.data, normalize_reference = TRUE,
                            reference_clusters = immune.combined@ident[colnames(ctrl.data)], 
                            logFC = 0.5, only_pos = T, estimate_weights_from_target = T)
scID_labels <- scID_res$labels
table(scID_labels)
table(immune.combined@ident[colnames(ctrl.data)])

scID_ARI <- mclust::adjustedRandIndex(scID_labels[names(which(scID_labels != "unassigned"))], immune.combined@ident[names(which(scID_labels != "unassigned"))])
scID_VI <- mcclust::vi.dist(scID_labels[names(which(scID_labels != "unassigned"))], immune.combined@ident[names(which(scID_labels != "unassigned"))])

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent SS2 cells with scmap
# ------------------------------------------------------------------------------------------------------------------------
scmap_labels <- run_scmap(reference_gem = log(ctrl.data.cpm + 1), reference_labels = immune.combined@ident[colnames(ctrl.data)], 
                          target_gem = log(stim.data.cpm + 1), n_features = 50)
table(scmap_labels)

scmap_ARI <- mclust::adjustedRandIndex(immune.combined@ident[names(which(scmap_labels != "unassigned"))], scmap_labels[names(which(scmap_labels != "unassigned"))])
scmap_VI <- mcclust::vi.dist(immune.combined@ident[names(which(scmap_labels != "unassigned"))], scmap_labels[names(which(scmap_labels != "unassigned"))])

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent cells with CaSTLe
# ------------------------------------------------------------------------------------------------------------------------
castle_labels <- runCastle(target_gem = log(stim.data.cpm + 1), source_gem = log(ctrl.data.cpm + 1), 
                           source_identities = immune.combined@ident[colnames(ctrl.data)])

castle_ARI <- mclust::adjustedRandIndex(castle_labels, immune.combined@ident[names(castle_labels)])
castle_VI <- mcclust::vi.dist(castle_labels, immune.combined@ident[names(castle_labels)])

# ------------------------------------------------------------------------------------------------------------------------
# Find equivalent cells with MNN
# ------------------------------------------------------------------------------------------------------------------------
mnn_res <- runMNN(target_gem = log(stim.data.cpm + 1), reference_gem = log(ctrl.data.cpm + 1))
mnn_labels <- mnn_res$target

MNN_ARI <- mclust::adjustedRandIndex(mnn_labels, immune.combined@ident[names(mnn_labels)])
MNN_VI <- mcclust::vi.dist(mnn_labels, immune.combined@ident[names(mnn_labels)])

# ------------------------------------------------------------------------------------------------------------------------
# Plot ARI
# ------------------------------------------------------------------------------------------------------------------------
df <- data.frame(ARI = c(scID_ARI, scmap_ARI, castle_ARI, MNN_ARI),
                 VI = c(scID_VI, scmap_VI, castle_VI, MNN_VI),
                 method=c("scID", "scmap", "CaSTLe", "MNN"))

ggplot(df, aes(y=ARI, x=method, fill=method))  + 
  geom_bar(stat="identity", position = "dodge") + scale_fill_jco()  + scale_y_continuous(limits = c(0,1)) + 
  theme(legend.position="top", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(title = "", x="", y="ARI")



