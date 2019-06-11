library(kBET)
library(ggplot2)

# ------------------------------------------------------------------------------------------------------------------------
# Pair of Figure 2: Shekhar retinal DropSeq vs SS2
# ------------------------------------------------------------------------------------------------------------------------

Ds_gem <- readRDS("~/scID_manuscript_figures/data/Figure2/Reference_gem.rds")
ss2_gem <- readRDS("~/scID_manuscript_figures/data/Figure2/Target_gem.rds")

common_genes <- intersect(rownames(Ds_gem), rownames(ss2_gem))
combined_gem <- cbind(Ds_gem[common_genes, ], ss2_gem[common_genes, ])
batch_labels <- c(rep(1, ncol(Ds_gem)), rep(2, ncol(ss2_gem)))
names(batch_labels) <- c(colnames(Ds_gem), colnames(ss2_gem))
kBET_Shekhar <- kBET(t(combined_gem), batch = batch_labels)
plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(kBET_Shekhar$stats$kBET.observed)), 
                        data =  c(kBET_Shekhar$stats$kBET.observed,
                                  kBET_Shekhar$stats$kBET.expected))
g_Shekhar <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
  labs(x='', y='Rejection rate',title='Shekhar') +
  theme_bw() +  
  scale_y_continuous(limits=c(0,1))

# ------------------------------------------------------------------------------------------------------------------------
# Pair of Figure 3
# ------------------------------------------------------------------------------------------------------------------------
E18_9k <- sc_utils.load10x("~/Google Drive/Data/singlecell/normal/10x/E18_braincells_9k/filtered_gene_bc_matrices/mm10/")
E18_nuclei <- sc_utils.load10x("~/Google Drive/Data/singlecell/normal/10x/E18_brainnuclei_1k/filtered_gene_bc_matrices/mm10/")
# Edit colnames of nuclei data to distinguish from E18_9k
colnames(E18_nuclei) <- paste(colnames(E18_nuclei), "nuc", sep = "_")

common_genes <- intersect(rownames(E18_9k), rownames(E18_nuclei))
combined_gem <- cbind(E18_9k[common_genes, ], E18_nuclei[common_genes, ])
batch_labels <- c(rep(1, ncol(E18_9k)), rep(2, ncol(E18_nuclei)))
names(batch_labels) <- c(colnames(E18_9k), colnames(E18_nuclei))

kBET_10Xbrain <- kBET(t(combined_gem), batch = batch_labels)
plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(kBET_10Xbrain$stats$kBET.observed)), 
                        data =  c(kBET_10Xbrain$stats$kBET.observed,
                                  kBET_10Xbrain$stats$kBET.expected))
g_10Xbrain <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
  labs(x='', y='Rejection rate',title='10X mouse brain') +
  theme_bw() +  
  scale_y_continuous(limits=c(0,1))

# ------------------------------------------------------------------------------------------------------------------------
# Pair of Figure 1E
# ------------------------------------------------------------------------------------------------------------------------
ctrl.data <- read.table(file = "~/scID_manuscript_figures/data/Kang_data/immune_control_expression_matrix.txt.gz", sep = "\t")
stim.data <- read.table(file = "~/scID_manuscript_figures/data/Kang_data/immune_stimulated_expression_matrix.txt.gz", sep = "\t")

common_genes <- intersect(rownames(ctrl.data), rownames(stim.data))
# Subsample cells as it is not possible to run for whole datasets due to memory limits
ctrl_sbst <- sample(colnames(ctrl.data), 2000)
stim_sbst <- sample(colnames(stim.data), 2000)
combined_gem <- cbind(ctrl.data[common_genes, ctrl_sbst], stim.data[common_genes, stim_sbst])
batch_labels <- c(rep(1, length(ctrl_sbst)), rep(2, length(stim_sbst)))
names(batch_labels) <- c(ctrl_sbst, stim_sbst)

kBET_Kang <- kBET(t(combined_gem), batch = batch_labels)

plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(kBET_Kang$stats$kBET.observed)), 
                        data =  c(kBET_Kang$stats$kBET.observed,
                                  kBET_Kang$stats$kBET.expected))
g_Kang <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
  labs(x='', y='Rejection rate',title='Kang') +
  theme_bw() +  
  scale_y_continuous(limits=c(0,1))


# ------------------------------------------------------------------------------------------------------------------------
# Pair of Figure 1F,G
# ------------------------------------------------------------------------------------------------------------------------
pancreas.data <- readRDS(file = "~/Desktop/scID_benchmarking_datasets/data/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata <- readRDS(file = "~/Desktop/scID_benchmarking_datasets/data/pancreas_v3_files/pancreas_metadata.rds")
true_labels <- metadata$celltype
names(true_labels) <- rownames(metadata)  

library(Seurat)
pancreas <- CreateSeuratObject(pancreas.data, meta.data = metadata)
pancreas.list <- SplitObject(pancreas, split.by = "tech")

celseq_gem <- as.data.frame(GetAssayData(pancreas.list[["celseq"]], slot = "counts"))
celseq2_gem <- as.data.frame(GetAssayData(pancreas.list[["celseq2"]], slot = "counts"))
ss2_gem <- as.data.frame(GetAssayData(pancreas.list[["smartseq2"]], slot = "counts"))
ss2_gem_cpm <- log2(scID::counts_to_cpm(ss2_gem)+1)

# SS2 - CEL-Seq pair
common_genes <- intersect(rownames(ss2_gem_cpm), rownames(celseq_gem))
combined_gem <- cbind(ss2_gem_cpm[common_genes, ], celseq_gem[common_genes, ])
batch_labels <- c(rep(1, ncol(ss2_gem_cpm)), rep(2, ncol(celseq_gem)))
names(batch_labels) <- c(colnames(ss2_gem_cpm), colnames(celseq_gem))

kBET_celseq <- kBET(t(combined_gem), batch = batch_labels)

plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(kBET_celseq$stats$kBET.observed)), 
                        data =  c(kBET_celseq$stats$kBET.observed,
                                  kBET_celseq$stats$kBET.expected))
g_celseq <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
  labs(x='', y='Rejection rate',title='CEL-Seq') +
  theme_bw() +  
  scale_y_continuous(limits=c(0,1))

# SS2 - CEL-Seq2 pair
common_genes <- intersect(rownames(ss2_gem_cpm), rownames(celseq2_gem))
combined_gem <- cbind(ss2_gem_cpm[common_genes, ], celseq2_gem[common_genes, ])
batch_labels <- c(rep(1, ncol(ss2_gem_cpm)), rep(2, ncol(celseq2_gem)))
names(batch_labels) <- c(colnames(ss2_gem_cpm), colnames(celseq2_gem))

kBET_celseq2 <- kBET(t(combined_gem), batch = batch_labels)

plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(kBET_celseq2$stats$kBET.observed)), 
                        data =  c(kBET_celseq2$stats$kBET.observed,
                                  kBET_celseq2$stats$kBET.expected))
g_celseq2 <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
  labs(x='', y='Rejection rate',title='CEL-Seq2') +
  theme_bw() +  
  scale_y_continuous(limits=c(0,1))

# ------------------------------------------------------------------------------------------------------------------------
# Plot results
# ------------------------------------------------------------------------------------------------------------------------
gridExtra::grid.arrange(g_Kang, g_10Xbrain, g_Shekhar, g_celseq, g_celseq2, ncol=5)
