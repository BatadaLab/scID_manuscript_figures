library(ggplot2)
library(scID)

ARI <- c()

# ------------------------------------------------------------------------------------------------------------------------
# Using Tirosh 2016
# ------------------------------------------------------------------------------------------------------------------------
gem <- readRDS("~/scID_manuscript_figures/data/Figure1/Tirosh2016_gem.rds")
labels <- readRDS("~/scID_manuscript_figures/data/Figure1/Tirosh2016_labels.rds")

scID_output <- scid_multiclass(target_gem = gem, reference_gem = gem, reference_clusters = labels, 
                               logFC = 0.5, only_pos = T, estimate_weights_from_target = T)
scID_labels <- scID_output$labels

ARI[1] <- mclust::adjustedRandIndex(labels[names(which(scID_labels != "unassigned"))], scID_labels[names(which(scID_labels != "unassigned"))])

# ------------------------------------------------------------------------------------------------------------------------
# Using Montoro 2018
# ------------------------------------------------------------------------------------------------------------------------
gem <- readRDS("~/scID_manuscript_figures/data/Figure1/Montoro2018_gem.rds")
labels <- readRDS("~/scID_manuscript_figures/data/Figure1/Montoro2018_labels.rds")

scID_output <- scid_multiclass(target_gem = gem, reference_gem = gem, reference_clusters = labels, 
                               logFC = 0.5, estimate_weights_from_target = T, only_pos = T)
scID_labels <- scID_output$labels

ARI[2] <- mclust::adjustedRandIndex(labels[names(which(scID_labels != "unassigned"))], scID_labels[names(which(scID_labels != "unassigned"))])

# ------------------------------------------------------------------------------------------------------------------------
# Using Hu 2017
# ------------------------------------------------------------------------------------------------------------------------
gem <- readRDS("~/scID_manuscript_figures/data/Figure1/Hu2017_gem.rds")
labels <- readRDS("~/scID_manuscript_figures/data/Figure1/Hu2017_labels.rds")

scID_output <- scid_multiclass(target_gem = gem, reference_gem = gem, reference_clusters = labels, 
                               logFC = 0.3, estimate_weights_from_target = T, only_pos = T)
scID_labels <- scID_output$labels

ARI[3] <- mclust::adjustedRandIndex(labels[names(which(scID_labels != "unassigned"))], scID_labels[names(which(scID_labels != "unassigned"))])

# ------------------------------------------------------------------------------------------------------------------------
# Using Shekar 2016
# ------------------------------------------------------------------------------------------------------------------------
gem <- readRDS("~/scID_manuscript_figures/data/Figure2/Reference_gem.rds")
labels <- readRDS("~/scID_manuscript_figures/data/Figure2/Reference_clusters.rds")

scID_output <- scid_multiclass(target_gem = gem, reference_gem = gem, reference_clusters = labels, 
                               logFC = 0.7, estimate_weights_from_target = T, only_pos = T)
scID_labels <- scID_output$labels

ARI[4] <- mclust::adjustedRandIndex(labels[names(which(scID_labels != "unassigned"))], scID_labels[names(which(scID_labels != "unassigned"))])

# ------------------------------------------------------------------------------------------------------------------------
# Plot results
# ------------------------------------------------------------------------------------------------------------------------
df <- data.frame(value=ARI, dataset = c("Tirosh", "Montoro", "Hu", "Shekhar"))
ggplot(df, aes(y=value, x=dataset)) + scale_y_continuous(limits = c(0, 1)) +
  geom_bar(stat="identity", position = "dodge") + 
  theme(legend.position="none", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm")) + 
  labs(title = "", x="", y="")


