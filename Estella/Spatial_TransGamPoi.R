library(scry)
library(BayesSpace)
library(spatialLIBD)

num_features = 3000

# Get data - DLPFC
# Connect to ExperimentHub
ehub <- ExperimentHub::ExperimentHub()
# Download the full real data (about 2.1 GB in RAM) use:
spe_all <- spatialLIBD::fetch_data(type = "spe", eh = ehub)

# Create three spe objects, one per sample:
spe1 <- spe_all[, colData(spe_all)$sample_id == '151673']
spe2 <- spe_all[, colData(spe_all)$sample_id == '151674']
rm(spe_all)

# For covenience, we use “gene names” instead of “gene ids”:
rownames(spe1) <- rowData(spe1)$gene_name
rownames(spe2) <- rowData(spe2)$gene_name

# Specify column names of spatial coordinates in colData(spe)
coordinates <- c("array_row", "array_col")
# Specify column names of spatial clusters in colData(spe)
spatial_cluster <- 'layer_guess_reordered'


# QC --------------------------------------------------
spe1 <- scuttle::addPerCellQC(spe1,)
# Remove combined set of low-quality spots
spe1 <- spe1[, !(colData(spe1)$sum < 25 |
                   colData(spe1)$detected < 25 |
                   colData(spe1)$expr_chrM_ratio > 0.3|
                   colData(spe1)$cell_count > 15)]

spe2 <- scuttle::addPerCellQC(spe2,)
# Remove combined set of low-quality spots
spe2 <- spe2[, !(colData(spe2)$sum < 25 |
                   colData(spe2)$detected < 25 |
                   colData(spe2)$expr_chrM_ratio > 0.3|
                   colData(spe2)$cell_count > 15)]

# For each sample i:
for(i in seq_len(2)){
  spe_i <- eval(parse(text = paste0("spe", i)))
  # Select QC threshold for lowly expressed genes: at least 20 non-zero spots:
  qc_low_gene <- rowSums(assays(spe_i)$counts > 0) >= 20
  # Remove lowly abundant genes
  spe_i <- spe_i[qc_low_gene,]
  assign(paste0("spe", i), spe_i)
  message("Dimension of spe", i, ": ", dim(spe_i)[1], ", ", dim(spe_i)[2])
}

# Dimension of spe1: 15116, 3611
# Dimension of spe2: 15879, 3655

UMI1 <- assay(spe1, "counts")
UMI2 <- assay(spe1, "counts")

# Transform (Normalize) ----------------------------------------





# Get variable features ----------------------------------------
get_high_dev_features <- function(UMI){
  dev <- scry::devianceFeatureSelection(UMI)
  var.genes <- rownames(UMI)[order(dev,decreasing=TRUE)[1:num_features]]

  return(var.genes)
}

get_high_exp_features <- function(UMI){
  dev <- scry::devianceFeatureSelection(UMI)
  var.genes <- rownames(UMI)[order(dev,decreasing=TRUE)[1:num_features]]

  return(var.genes)
}

get_high_var_features <- function(UMI){
  dev <- scry::devianceFeatureSelection(UMI)
  var.genes <- rownames(UMI)[order(dev,decreasing=TRUE)[1:num_features]]

  return(var.genes)
}

get_spa_var_features <- function(UMI){
  dev <- scry::devianceFeatureSelection(UMI)
  var.genes <- rownames(UMI)[order(dev,decreasing=TRUE)[1:num_features]]

  return(var.genes)
}


# kNN calculation -----------------------------------------------
str(out <- BiocSingular::runPCA(a, rank=10, get.rotation = FALSE)$x) # 5000 rows, 10 cols
KNN <- BiocNeighbors::findAnnoy(out, k = 50, warn.ties = FALSE)$index # 5000 rows, 50 cols

str(out2 <- BiocSingular::runPCA(b, rank=10, get.rotation = FALSE)$x) # 5000 rows, 10 cols
KNN2 <- BiocNeighbors::findAnnoy(out2, k = 50, warn.ties = FALSE)$index # 5000 rows, 50 cols

KNNs <- list(KNN, KNN2)

# row here are the cells
n_cells <- nrow(KNNs[[1]])

cons <- mean(
  sapply(seq_len(n_cells), function(cell_idx){
    length(intersect(KNNs[[1]][cell_idx,], KNNs[[2]][cell_idx,]))
  })
)


# Canoncial correlation with first 50 PCs

