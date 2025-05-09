#ifndef LIBFUSE_H_
#define LIBFUSE_H_

#include <stddef.h> 
#include <math.h>  // For NAN

#define FUSE_NA_DOUBLE NAN

static inline int is_na_double(double x) {
    return isnan(x);
}

typedef struct {
	double k0, k1;
} data_elem_t;

typedef struct {
	double nllr, d;
} distance_t;

/* NB. size in doubles */

/**
 * @brief Computes x / y safely.
 * @param x Numerator.
 * @param y Denominator.
 * @return Result of x / y or 0.5 if y is zero.
 */
double bino_div(double x, double y);

/**
 * @brief Computes x / y with integer numerator.
 * @param x Integer numerator.
 * @param y Double denominator.
 * @return Result of x / y or 0.5 if y is zero.
 */
double bino_div_id(int x, double y);

/**
 * @brief Computes x * log(y) safely.
 * @param x Multiplicand.
 * @param y Value inside log function.
 * @return x * log(y) or 0 if y <= 0.
 */
double bino_xlogy(double x, double y);

/**
 * @brief Computes x * log(y) with integer x.
 * @param x Integer multiplicand.
 * @param y Double value inside log function.
 * @return x * log(y) or 0 if y is zero.
 */
double bino_xlogy_id(int x, double y);

/**
 * @brief Computes cost function based on binomial likelihood.
 * @param k0 Unmethylated counts.
 * @param k1 Methylated counts.
 * @return Negative binomial-likelihood.
 */
double cost(double k0, double k1);

/**
 * @brief Computes negative log-likelihood ratio for data elements.
 * @param K1 Lefthand cluster to merge.
 * @param K2 Righthand cluster to merge.
 * @param LOG_L1 Log-likelihood of lefthand cluster.
 * @param LOG_L2 Log-likelihood of righthand cluster.
 * @param n Number of samples in the data = length of vectors.
 * @return Negative log-likelihood ratio.
 */
double nLLR(const data_elem_t *K1, const data_elem_t *K2, double LOG_L1, double LOG_L2, size_t n);

/**
 * @brief Computes merged cost between two elements.
 * @param K1 Lefthand cluster to merge.
 * @param K2 Righthand cluster to merge.
 * @param LOG_L1 Log-likelihood of lefthand cluster.
 * @param LOG_L2 Log-likelihood of righthand cluster.
 * @param chr1 Chromosome of lefthand cluster.
 * @param chr2 Chromosome of righthand cluster.
 * @param pos1 Position of lefthand cluster.
 * @param pos2 Position of righthand cluster.
 * @param n Number of samples in the data = length of vectors.
 * @return Distance structure with negative log-likelihood ratio and genomic distance.
 */
distance_t merged_cost(const data_elem_t *K1, const data_elem_t *K2, double LOG_L1, double LOG_L2, const int chr1, const int chr2, const int pos1, const int pos2, size_t n);

/**
 * @brief Initializes the state for binomial clustering.
 * 
 * This function computes the binomial log-likelihoods for each site, initializes 
 * distance structures, and sets up data for hierarchical clustering.
 * 
 * @param D Array of distance structures to store pairwise distances and log-likelihood ratios.
 * @param L Array representing cluster labels, initialized with negative indices.
 * @param PL Array storing left neighbor indices for merging.
 * @param PR Array storing right neighbor indices for merging.
 * @param K Input data elements containing binomial counts for clustering.
 * @param LOG_L Array to store computed binomial log-likelihoods for each site.
 * @param tracked_cost Array to track the cumulative cost associated with each cluster.
 * @param chr Array of chromosome identifiers for genomic positions.
 * @param pos Array of genomic positions corresponding to data elements.
 * @param m Number of elements (sites) in the dataset.
 * @param n Number of samples in the data = length of vectors.
 */
void init_state(distance_t *D, double *L, size_t *PL, size_t *PR, const data_elem_t *K, double *LOG_L, double *tracked_cost, const int *chr, const int *pos, size_t m, size_t n);

/**
 * @brief Finds the best neighboring cluster to merge.
 * 
 * This function scans left or right from a given position `current` to find the 
 * neighboring cluster with the smallest combined genomic distance and 
 * negative log-likelihood ratio (nllr). It determines the best merge candidate 
 * based on the minimum cost.
 * 
 * @param p Index of the current cluster.
 * @param D Array of distance structures containing genomic distances and nllr values.
 * @param PL Array storing left neighbor indices for merging.
 * @param PR Array storing right neighbor indices for merging.
 * @return Index of the best neighboring cluster to merge with.
 */
size_t new_pair(size_t p, distance_t *D, size_t *PL, size_t *PR);

/**
 * @brief Returns the maximum of two double values.
 * @param A First double value.
 * @param B Second double value.
 * @return The larger of A and B, or A if they are equal.
 */
double max_d(double A, double B);

/**
 * @brief Updates cluster data by merging two data clusters.
 * 
 * This function merges the counts of clusters (`p1` and `p2`) 
 * in the data array `K`. 
 * It also updates the corresponding log-likelihood values in `LOG_L`.
 * 
 * @param K Array of counts.
 * @param LOG_L Array of log-likelihood values.
 * @param p1 Index of the lefthand cluster (destination of merge).
 * @param p2 Index of the righthand cluster.
 * @param n Number of samples in the data = length of vectors.
 */
void update_counts(data_elem_t *K, double *LOG_L, size_t p1, size_t p2, size_t n);

/**
 * @brief Cuts a hierarchical clustering tree into `k` clusters.
 *
 * This function assigns cluster labels to `n` data points based on a hierarchical 
 * clustering tree (`Z`). It processes the tree by assigning labels to each leaf, 
 * propagating labels from merged nodes, and compressing labels into a range of `[0, k)`.
 *
 * @param[out] L Pointer to an array of size `n` that will store cluster labels.
 * @param[in] Z Pointer to a `double` array representing the hierarchical clustering tree.
 * @param[in] k The desired number of clusters (must be between `1` and `n`).
 * @param[in] n The number of original data points.
 */
void cutree(size_t *L, const double *Z, size_t k, size_t n);

/**
 * @brief Computes the size of a sorted tree structure.
 * @param m The number of rows in the tree.
 * @return The computed tree size, which is 5 times `m`.
 */
size_t sort_tree_size(size_t m);

/**
 * @brief Sorts a hierarchical clustering tree.
 *
 * This function reorders the hierarchical clustering tree stored in `Z` and outputs the sorted tree in `Z_new`. 
 * Both trees have 5 columns:
 * 			0. label for merged point 1
 * 			1. label for merged point 2
 * 			2. Change in total likelihood for forming this merge
 * 			3. Total cost of points in cluster
 * 			4. Genomic distance between points
 * The sorting is based on columns 2 and 3.
 *
 * @param[out] Z_new Pointer to an array where the sorted tree will be stored.
 * @param[in] scratch Pointer to preallocated scratch space used for temporary storage.
 * @param[in] Z Pointer to the original hierarchical clustering tree data.
 * @param[in] m The number of rows in the tree.
 */
void sort_tree(double *Z_new, void *scratch, const double *Z, size_t m);

/**
 * @brief Performs hierarchical clustering of CpG sites based on genomic distance and negative log-likelihood ratio.
 *
 * This function implements a hierarchical clustering algorithm for CpG sites, using likelihood-based merging. 
 * The resulting clustering tree is stored in `tree` and must be sorted using `sort_tree()` after execution.
 *
 * @param[out] tree Array storing the resulting clustering tree. Each row represents a merge:
 *                     - Column 0: Identifier of first merged point.
 *                     - Column 1: Identifier of second merged point.
 *                     - Column 2: Change in total likelihood for the merge.
 *                     - Column 3: Total likelihood for the cluster.
 *                     - Column 4: Penalty for genomic distance between points.
 * @param[in] counts Structure holding count matrix (num_of_sites × num_of_samples) containing unmethylated (k0) and methylated (k1) count values for each CpG site.
 * @param[in] chr Array of size `m`, indicating the chromosome each CpG site belongs to.
 * @param[in] pos Array of size `m`, indicating the genomic coordinate of each CpG site.
 * @param[in] num_of_sites Number of CpG sites (rows in `K0` and `K1`).
 * @param[in] num_of_samples Number of samples (columns in `K0` and `K1`).
 */
void fuse_cluster(double *tree, data_elem_t *counts, const int *chr, const int *pos, size_t num_of_sites, size_t num_of_samples);

/**
 * @brief Computes Pearson correlation row sums.
 * @param sxy Array to store the sum of products of x and Y rows.
 * @param sy1 Array to store the sum of Y rows.
 * @param sy2 Array to store the sum of squares of Y rows.
 * @param x Array representing the x vector.
 * @param Y Matrix representing the Y vectors.
 * @param m Number of rows in Y.
 * @param ny Number of columns in Y.
 */
void cor_pearson_rowSums(double *sxy, double *sy1, double *sy2, const double *x, const double *Y, size_t m, size_t ny);

#endif /* LIBISLAND_H_ */