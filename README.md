# MatrixCor.Perm
Permutation tests for estimating correlation between two symmetric matrices. Rows and columns are permuted instead of elements.

This script contains the functions performing the permutation tests for correlations between two symmetric matrices. Unlike the standard permutation test, the columns and rows of the matrix are permuted instead of individual elements (according to the algorithm described in Manly, 1985, Section 6.8).

MatrixCor.Perm contains three functions.

1. MatrixCor.AllPerm
This function computes a two-tailed permutation test of a Pearson correlation coefficient between two symmetric matrices. Diagonal elements are neglected. Rows and columns are permuted instead of elements (according to Manly, 1985, Section 6.8). All possible permutations are used. Therefore, running the test for large matrix sizes may take a long time.
The function requires combinat R package.

2. MatrixCor.RandPerm
This function computes a two-tailed permutation test of a Pearson correlation coefficient between two symmetric matrices. Diagonal elements are neglected. 
Rows and columns are permuted instead of elements (according to Manly, 1985, Section 6.8). 1000 permutations are used by default: non-permuted data and 999 random permutations.

3. MatrixCor.PhyloPerm
This function computes a two-tailed permutation test of a Pearson correlation coefficient between two symmetric matrices using the phylogenetic information. Diagonal elements are neglected. Rows and columns are permuted instead of elements (according to Manly, 1985, Section 6.8). 1000 permutations are used by default: non-permuted data and 999 random permutations.
Probabilities for phylogenetic permutations are calculated using the "EVO_973_sm_PhyloMantel.R" function (Harmon and Glor, 2010) based on the algorithm developed by Lapointe and Garland, 2001. A parameter k determining the weighting given to the permutation probabilities is set by default to 2.
The function requires ape R package and the "EVO_973_sm_PhyloMantel.R" function (Harmon and Glor, 2010).

References
Harmon L. J., Glor R. E. Poor statistical performance of the Mantel test in phylogenetic comparative analyses // Evolution. 2010. V. 64, No. 7. P. 2173-2178. https://doi.org/10.1111/j.1558-5646.2010.00973.x

Lapointe F.-J., Garland T. A generalized permutation model for the analysis of cross-specific data // Journal of Classification. 2001. V. 18. P. 109-127. https://doi.org/10.1007/s00357-001-0007-0

Manly B. F. J. The statistics of natural selection on animal populations. London – New York: Chapman and Hall. 1985. xvi+484 p.
