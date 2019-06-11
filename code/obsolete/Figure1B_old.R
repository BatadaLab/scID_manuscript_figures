library(scID)

# ----------------------------------------------------------------------------------
# Read data
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

# ----------------------------------------------------------------------------------
# Find markers
markers <- find_markers(gem, labels, 0.7, only.pos=FALSE, normalize_reference = TRUE)
celltypes <- unique(markers$cluster)

# ----------------------------------------------------------------------------------
# Bin values to 0 and 1 for present (expressed) and absent genes
sink("aux");
binned_gem <- apply(gem, 1, function(x) biomod2::BinaryTransformation(x, threshold = quantile(x[which(x>0)], 0.25, na.rm = TRUE)))
sink(NULL);
rownames(binned_gem) <- colnames(gem)
binned_gem <- binned_gem[complete.cases(binned_gem), ]

# ----------------------------------------------------------------------------------
# Normalize gem
gem_norm <-  t(apply(gem[unique(markers$gene), ], 1, function(x) normalize_gene(x)))
gem_norm <- gem_norm[complete.cases(gem_norm), ]

# ----------------------------------------------------------------------------------
# Step 2: find training sets and calculate weights
celltype <- "MG (Mueller Glia)"

celltype_markers <- markers[which(markers$cluster == celltype), ]
positive_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC > 0)]
negative_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC < 0)]

# Find total number of expressed genes per cell (n_e)
n_e <- rowSums(binned_gem)
# Find total number of expressed positive marker genes per cell (n_pme)
if (length(positive_markers) >= 1) {
  n_pme <- rowSums(binned_gem[, positive_markers, drop = FALSE])
} else {
  n_pme <- rep(0, nrow(binned_gem))
  names(n_pme) <- rownames(binned_gem)
}
# Find total number of expressed negative marker genes per cell (n_nme)
if (length(negative_markers) >= 1) {
  n_nme <- rowSums(binned_gem[, negative_markers, drop = FALSE])
} else {
  n_nme <- rep(0, nrow(binned_gem))
  names(n_nme) <- rownames(binned_gem)
}
# Find total number of positive marker genes (n_pm)
n_pm <- max(length(positive_markers), 1)
# Find total number of negative marker genes (n_nm)
n_nm <- max(length(negative_markers), 1)

data <- data.frame(recall = (n_pme/n_pm) - (n_nme/n_nm), precision = (n_pme-n_nme)/n_e)
rownames(data) <- colnames(gem)
data <- data[complete.cases(data), ]

library(mclust)
fit <- Mclust(data)

plot(fit)


# Get centroids of each cluster
centroids <- data.frame(matrix(NA, length(unique(fit$classification)), 2), row.names = unique(fit$classification))
colnames(centroids) <- c("precision", "recall")
sds <- data.frame(matrix(NA, length(unique(fit$classification)), 2), row.names = unique(fit$classification))
colnames(sds) <- c("precision", "recall")
for (ID in rownames(centroids)) {
  centroids[ID, "precision"] <- mean(data[which(fit$classification == ID), "precision"]) 
  sds[ID, "precision"] <- sd(data[which(fit$classification == ID), "precision"]) 
  centroids[ID, "recall"] <- mean(data[which(fit$classification == ID), "recall"]) 
  sds[ID, "recall"] <- sd(data[which(fit$classification == ID), "recall"]) 
}

IN_candidates <- unique(c(rownames(centroids)[which(centroids$recall == max(centroids$recall))], rownames(centroids)[which(centroids$precision == max(centroids$precision))]))
E_dist <- apply(centroids, 1, function(x) sqrt((1-x[1])^2 + (1-x[2])^2))
IN_id <- names(E_dist)[which(E_dist == min(E_dist))]

IN_cells <- colnames(gem)[which(fit$classification %in% IN_id)]
OUT_cells <- colnames(gem)[which(fit$classification %in% setdiff(unique(fit$classification), IN_candidates))]

df <- data[c(IN_cells, OUT_cells), ]
df$label <- rep("negative", nrow(df))
df[IN_cells, "label"] <- "positive"
df$true_labels <- labels[rownames(df)]
df$true_labels[which(df$true_labels != celltype)] <- "Other"
df$new <- paste(df$true_labels, df$label, sep = "_")
df$new <- factor(df$new, levels = c(paste(celltype, "positive", sep = "_"), 
                                    "Other_positive", "Other_negative",
                                    paste(celltype, "negative", sep = "_"))) 

PR_plot <- ggplot(df, aes(x=precision, y=recall, color=new)) + scale_color_jco() + 
  geom_point(size=0.5) + 
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(title = "", x="precision", y="recall")


gene.weights <- scID_weight(gem_norm[c(positive_markers, negative_markers), , drop=FALSE], IN_cells, OUT_cells)

# ----------------------------------------------------------------------------------
signature <- intersect(names(gene.weights), rownames(gem_norm))
weighted_gem <- gene.weights[signature] * gem_norm[signature, ,drop=FALSE]

score <- colSums(weighted_gem)/sqrt(sum(gene.weights^2))
fit <- mclust::densityMclust(score)
plot(fit)

avgScore <- rep(NA, length(unique(fit$classification)))
names(avgScore) <- unique(fit$classification)
for (ID in names(avgScore)) {
  avgScore[ID] <- mean(score[names(which(fit$classification == ID))])
}

matches <- names(fit$classification)[which(fit$classification == names(which(avgScore == max(avgScore))))]
scID_labels <- rep("Other", length(score))
names(scID_labels) <- names(score)
scID_labels[matches] <- celltype
true_labels <- labels
true_labels[which(true_labels != celltype)] <- "Other"

df <- data.frame(score=score, true_label=true_labels[names(score)], predicted_label=scID_labels[names(score)])
df$new <- paste(df$true_label, df$predicted_label, sep = "_")
df$new <- factor(df$new, levels = c(paste(celltype, celltype, sep = "_"), 
                                    paste(celltype, "Other", sep = "_"),
                                    "Other_Other"))

score_plot <- ggplot(df, aes(x=score, fill=new)) + geom_histogram(bins=100) + scale_fill_jco() + 
  theme(legend.position="", text = element_text(size=3), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  #scale_y_log10() +
  labs(title = "", x="score", y="Number of cells")
score_plot

pdf("~/Dropbox/nizar_katerina2017/scID/cellreports/manuscript/Figures/Figure1/Fig1B.pdf", width = 2, height = 4)
gridExtra::grid.arrange(PR_plot, score_plot, ncol=1)
dev.off()
  



