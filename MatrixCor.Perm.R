# Permutation tests for estimating correlation between two symmetric matrices
# Rows and columns are permuted instead of elements (according to Manly, 1985, Section 6.8)

# by Igor Dzeverin, 2025-03-16




# MatrixCor.AllPerm

# Description
# This function computes a two-tailed permutation test of a Pearson 
# correlation coefficient between two symmetric matrices.
# Diagonal elements are neglected. 
# Rows and columns are permuted instead of elements (according to Manly, 1985, Section 6.8).
# All possible permutations are used.

# Usage
# MatrixCor.AllPerm (X, Y)

# Arguments
# X and Y are symmetric matrices

# Value
# MatrixCor.AllPerm returns a vector containing 
# Pearson correlation, total number of permutations, and p-value for matrices X and Y

MatrixCor.AllPerm <- function(X,Y)
{
  require(combinat)
  cor.main <- cor(as.vector(as.dist(X)),as.vector(as.dist(Y)))
  
  s <- nrow(X)
  sf <- factorial(s)
  a <- 0
  P <- matrix(unlist(permn(s)),sf,s,byrow=TRUE)
  
  for (i in 1:sf)
  {
    y <- P[i,]
    Y.perm <- Y[y,y]
    cor.perm <- cor(as.vector(as.dist(X)),as.vector(as.dist(Y.perm)))
    if (abs(cor.perm) >= abs(cor.main)) a <- a+1
  }
  
  p.value <- a/sf
  result <- c(cor.main, sf, p.value)
  names(result) <- c("Pearson r", "Permut N", "p-value")
  return(result)
}




# MatrixCor.RandPerm

# Description
# This function computes a two-tailed permutation test of a Pearson 
# correlation coefficient between two symmetric matrices.
# Diagonal elements are neglected. 
# Rows and columns are permuted instead of elements (according to Manly, 1985, Section 6.8).
# 1000 permutations are used by default,
# non-permuted data and 999 random permutations.

# Usage
# MatrixCor.RandPerm (X, Y, nperm=999)

# Arguments
# X and Y are symmetric matrices
# nperm is a number of random permutations (total number of permutations is nperm+1)

# Value
# MatrixCor.RandPerm returns a vector containing
# Pearson correlation, total number of permutations, and p-value for matrices X and Y

MatrixCor.RandPerm <- function(X,Y,nperm=999)
{
  cor.main <- cor(as.vector(as.dist(X)),as.vector(as.dist(Y)))
  a <- 1
  s <- nrow(X)
  
  for (i in 1:nperm)
  {
    y <- sample(s)
    Y.perm <- Y[y,y]
    cor.perm <- cor(as.vector(as.dist(X)),as.vector(as.dist(Y.perm)))
    if (abs(cor.perm) >= abs(cor.main)) a <- a+1
  }
  
  p.value <- a/(nperm+1)
  result <- c(cor.main, nperm+1, p.value)
  names(result) <- c("Pearson r", "Permut N", "p-value")
  return(result)  
}




# MatrixCor.PhyloPerm

# Description
# This function computes a two-tailed permutation test of a Pearson 
# correlation coefficient between two symmetric matrices using the phylogenetic
# information.
# Diagonal elements are neglected. 
# Rows and columns are permuted instead of elements (according to Manly, 1985, Section 6.8).
# 1000 permutations are used by default:
# non-permuted data and 999 random permutations.
# Probabilities for phylogenetic permutations are calculated using 
# "EVO_973_sm_PhyloMantel.R" program (Harmon and Glor 2010) based on the 
# algorithm developed by Lapointe and Garland 2001.
# A parameter k determining the weighting given to the permutation probabilities
# is set by default to 2.

# Usage
# MatrixCor.PhyloPerm (X, Y, phy, k=2, nperm=999)

# Arguments
# X and Y are symmetric matrices
# phy is phylogenetic tree
# k is a parameter determining the weighting given to the permutation probabilities
# nperm is a number of random permutations (total number of permutations is nperm+1)

# Value
# MatrixCor.PhyloPerm returns a vector containing
# Pearson correlation, total number of permutations, k parameter, and p-value for matrices X and Y

MatrixCor.PhyloPerm <- function(X,Y,phy,k=2,nperm=999)
{
  require(ape)
  source("EVO_973_sm_PhyloMantel.R")
  
  cor.main <- cor(as.vector(as.dist(X)),as.vector(as.dist(Y)))
  a <- 1
  
  for (i in 1:nperm)
  {
    y <- phyloPermute(phy,k)
    Y.perm <- Y[y,y]
    cor.perm <- cor(as.vector(as.dist(X)),as.vector(as.dist(Y.perm)))
    if (abs(cor.perm) >= abs(cor.main)) a <- a+1
  }
  
  p.value <- a/(nperm+1)
  result <- c(cor.main, nperm+1, k, p.value)
  names(result) <- c("Pearson r", "Permut N", "k", "p-value")
  return(result)
}




# References
# Harmon L. J., Glor R. E. Poor statistical performance of the Mantel test in phylogenetic comparative analyses // Evolution. 2010. V. 64, No. 7. P. 2173-2178. https://doi.org/10.1111/j.1558-5646.2010.00973.x
# Lapointe F.-J., Garland T. A generalized permutation model for the analysis of cross-specific data // Journal of Classification. 2001. V. 18. P. 109-127. https://doi.org/10.1007/s00357-001-0007-0
# Manly B. F. J. The statistics of natural selection on animal populations. London – New York: Chapman and Hall. 1985. xvi+484 p.
