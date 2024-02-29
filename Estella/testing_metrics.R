# so input dat should be gene by cell
# BiocSingular::runPCA() takes a cell by gene matrix, returns a cell by PCA reddim
# irlba::prcomp_irlba()  takes a cell by gene matrix, returns a cell by PCA reddim
# stats::prcomp()        takes a cell by gene matrix, returns a cell by PCA reddim
# scater::runPCA()       takes a gene by cell matrix, returns a cell by PCA reddim
BiocSingular::runPCA(t(dat), rank = pca_dim, get.rotation = FALSE, BSPARAM = BiocSingular::FastAutoParam())$x
BiocNeighbors::findAnnoy(red_dat, k = k_nearest_neighbors, warn.ties = FALSE)$index

a <- matrix(rnorm(100000), ncol=20) # 5000 rows, 20 cols
str(out <- BiocSingular::runPCA(a, rank=10, get.rotation = FALSE)$x) # 5000 rows, 10 cols
KNN <- BiocNeighbors::findAnnoy(out, k = 50, warn.ties = FALSE)$index # 5000 rows, 50 cols

b <- matrix(rnorm(100000), ncol=20) # 5000 rows, 20 cols
str(out2 <- BiocSingular::runPCA(b, rank=10, get.rotation = FALSE)$x) # 5000 rows, 10 cols
KNN2 <- BiocNeighbors::findAnnoy(out2, k = 50, warn.ties = FALSE)$index # 5000 rows, 50 cols

KNNs <- list(KNN, KNN2)

# row here are the cells
n_cells <- nrow(KNNs[[1]])

cons <-
  # mean(
  sapply(seq_len(n_cells), function(cell_idx){
  length(intersect(KNNs[[1]][cell_idx,], KNNs[[2]][cell_idx,]))
})
# )


# irlba PCA calculation ---------------------------------------------------
library(irlba)
library(matrixStats)
set.seed(1)
x  <- matrix(rnorm(200), nrow=20) # 20 rows for cells, 10 cols for genes
p1 <- prcomp_irlba(x, n=3)$x # takes a cell by gene matrix, returns cell by PCs reddim # 20, 3
summary(p1)

# Compare with
p2 <- prcomp(x, tol=0.7)$x # takes a cell by gene matrix, returns cell by PCs reddim # 20, 3
summary(p2)


# In paper ----------------------------------------------------------------
set.seed(1)
x  <- matrix(rnorm(20000), nrow=200) # 200 rows for genes, 100 cols for cells
# p <- prcomp_irlba(x, n=10)$x;  # takes a cell by gene matrix, so need to transpose
# p <- prcomp_irlba(t(x), n=10)$x # 100 cells by 10 PCs

# PCA calculation
pca <- irlba::prcomp_irlba(t(x), n = 10)$x # 100 cells 10 PCs

size_factors <- colSums2(x) # 100 cells 1 col
size_factors <- size_factors / mean(size_factors)
summary(size_factors)

cancor(pca, size_factors)$cor # just one value, so no matter the order of x y
cancor(size_factors, pca)$cor # just one value, so no matter the order of x y




# Canonical correlation ---------------------------------------------------
# Compute the canonical correlations between two data matrices.
# The canonical correlation analysis seeks linear combinations of the y variables
# which are well explained by linear combinations of the x variables.
## signs of results are random
pop <- LifeCycleSavings[, 2:3] # 50 2
oec <- LifeCycleSavings[, -(2:3)] # 50 3
cancor(pop, oec) # $cor has two values, I guess each value for a col in min ncol of x = pop or y = oec

x <- matrix(rnorm(150), 50, 3) # 50 3
y <- matrix(rnorm(250), 50, 5) # 50 5
(cxy <- cancor(x, y)) # $cor has three values, I guess each value for a col in min ncol of x or y
all(abs(cor(x %*% cxy$xcoef,
            y %*% cxy$ycoef)[,1:3] - diag(cxy $ cor)) < 1e-15)
all(abs(cor(x %*% cxy$xcoef) - diag(3)) < 1e-15)
all(abs(cor(y %*% cxy$ycoef) - diag(5)) < 1e-15)



# GLM-PCA attempt ---------------------------------------------------------
library(glmpca)

#create a simple dataset with two clusters
#mu<-rep(c(.5,3),each=10)
mu<-matrix(exp(rnorm(100*20)),nrow=100)
mu[,1:10]<-mu[,1:10]*exp(rnorm(100))
clust<-rep(c("red","black"),each=10)
Y<-matrix(rpois(prod(dim(mu)),mu),nrow=nrow(mu))

#visualize the latent structure
res<-glmpca(Y, 2)
factors<-res$factors
plot(factors[,1],factors[,2],col=clust,pch=19)


library(scry)

# Get variable features
dev <- scry::devianceFeatureSelection(data)
var.genes <- rownames(data)[order(dev,decreasing=TRUE)[1:num_features]]








