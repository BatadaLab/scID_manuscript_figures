run_scmap <- function(reference_gem, reference_labels, target_gem, n_features = 500) {
  
  library(SingleCellExperiment)
  library(scmap)
  # Keep only genes that are non-zero across cells
  reference_gem <- reference_gem[which(rowSums(reference_gem) != 0), ]
  target_gem <- target_gem[which(rowSums(target_gem) != 0), ]
  
  common_genes <- intersect(rownames(reference_gem), rownames(target_gem))
  
  ann <- data.frame(cell_type1=as.factor(reference_labels))
  
  # Create index
  ref_sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(reference_gem[common_genes, ])), colData = ann)
  # use gene names as feature symbols
  rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
  # remove features with duplicated names
  ref_sce <- ref_sce[!duplicated(rownames(ref_sce)), ]
  ref_sce <- selectFeatures(ref_sce, suppress_plot = TRUE, n_features = n_features)
  ref_sce <- indexCluster(ref_sce)
  
  # Create sce for testing data
  sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(target_gem[common_genes, ])))
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
  names(scmap_labels) <- colnames(target_gem)
  
  return(scmap_labels)
}

runCastle <- function(target_gem, source_gem, source_identities) {
  library(scater)  
  library(xgboost) 
  library(igraph)  
  
  BREAKS=c(-1, 0, 1, 6, Inf)
  nFeatures = 200
  
  # Keep only reference cells with available labels
  common_cells <- intersect(colnames(source_gem), names(source_identities))
  
  source = SingleCellExperiment(assays = list(logcounts = as.matrix(source_gem[, common_cells])), 
                                colData = data.frame(cell_type1 = source_identities[common_cells]))
  
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

runMNN <- function(target_gem, reference_gem) {
  
  library(SingleCellExperiment)
  library(scran)
  common_genes <- intersect(rownames(reference_gem), rownames(target_gem))
  mnn_corrected <- mnnCorrect(as.matrix(reference_gem[common_genes, ]), as.matrix(target_gem[common_genes, ]))
  
  # Cluster mnn_corrected gem with Seurat
  gem <- as.data.frame(mnn_corrected$corrected)
  colnames(gem)[grep("X", colnames(gem))] <- colnames(target_gem)
  
  seurat <- CreateSeuratObject(counts = gem)
  seurat <- FindVariableFeatures(seurat)
  seurat <- ScaleData(seurat)
  seurat <- RunPCA(seurat,  do.print = FALSE)
  seurat <- FindNeighbors(seurat)
  seurat <- FindClusters(seurat)
  
  mnn_labels_ref <- Idents(seurat)[colnames(reference_gem)]
  mnn_labels_target <- Idents(seurat)[colnames(target_gem)]
  
  return(list(reference_labels=mnn_labels_ref, target_labels=mnn_labels_target))
}
