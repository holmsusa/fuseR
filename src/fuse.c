/**
 * @file fuse.c
 * @brief Implementation of R-C interface functions.
 * 
 * This file contains functions that take in R variables and return R variables, 
 * but are implemented in C.
 * The functions in this file are for clustering, sorting trees, cutting trees, 
 * and computing Pearson correlation row sums in an R environment.
 */


#include <Rinternals.h>
#include <R.h>
#include <assert.h>
#include "libfuse.h"


/**
 * @brief Checks if the given SEXP argument is a 2D structure.
 * 
 * @param arg An R object to check.
 * @return int Returns nonzero if the argument is a vector or matrix, otherwise 0.
 */
static
__attribute__((unused))
int
is2D(SEXP arg)
{ return isVector(arg) || isMatrix(arg); }


/**
 * @brief Sorts a tree structure in R.
 * 
 * @param tree_R Input matrix representing the tree structure.
 * @return SEXP A newly allocated matrix with the sorted tree.
 */
SEXP
sort_tree_R(SEXP tree_R)
{
	SEXP tree_sorted_R ; 
	size_t nrow;
	double *scratch; // Scratch memory

	tree_sorted_R = PROTECT(duplicate(tree_R));
	nrow = nrows(tree_R);

	scratch = (double*)R_alloc(sizeof(*scratch), sort_tree_size(nrow));

	sort_tree( REAL(tree_sorted_R), scratch, REAL(tree_R), nrow );

	UNPROTECT(1);
	return tree_sorted_R;
}


/**
 * @brief Performs hierarchical clustering on input matrices.
 * 
 * @param K0_R Integer matrix representing the unmethylated counts.
 * @param K1_R Integer matrix representing the methylated counts.
 * @param CHR_R Integer vector of chromosome identifiers.
 * @param POS_R Integer vector of position identifiers.
 * @return SEXP A matrix with 5 columns containing clustering results.
 * 			1. label for merged point 1
 * 			2. label for merged point 2
 * 			3. Change in total likelihood
 * 			4. Total cost of merge
 * 			5. Genomic distance between points
 */
SEXP 
fuse_cluster_R(SEXP K0_R, SEXP K1_R, SEXP CHR_R, SEXP POS_R)
{
    int nrow = nrows(K0_R);
    int ncol = ncols(K0_R);

    /* Validate input */
    if (!isInteger(K0_R) || !isInteger(K1_R))
        error("'K0' and 'K1' must be integer matrices");

    if (nrows(K1_R) != nrow || ncols(K1_R) != ncol)
        error("'K1' must have the same dimensions as 'K0'");

    /* Allocate output */
    SEXP tree_R = PROTECT(allocMatrix(REALSXP, nrow - 1, 5));

    /* Allocate temporary buffer for counts */
    data_elem_t *counts = (data_elem_t *) R_alloc(nrow * ncol, sizeof(data_elem_t));

    /* Fill in counts (transposed and NA-handled) */
    for (int j = 0; j < ncol; ++j) {
        for (int i = 0; i < nrow; ++i) {
            int idx = i + j * nrow;               // R index (column-major)
            int cidx = j + i * ncol;              // Transposed index (row-major)

            int k0_val = INTEGER(K0_R)[idx];
            int k1_val = INTEGER(K1_R)[idx];

            counts[cidx].k0 = (k0_val == NA_INTEGER) ? FUSE_NA_DOUBLE : (double)k0_val;
            counts[cidx].k1 = (k1_val == NA_INTEGER) ? FUSE_NA_DOUBLE : (double)k1_val;
        }
    }

    
    /* Call core clustering function */
    fuse_cluster(REAL(tree_R), counts, INTEGER(CHR_R), INTEGER(POS_R), nrow, ncol);

    UNPROTECT(1);
    return tree_R;
}


/**
 * @brief Cuts a hierarchical tree into clusters.
 * 
 * @param tree_R Input matrix representing the tree structure.
 * @param num_of_clust_R Number of clusters.
 * @return SEXP A vector of cluster labels.
 */
SEXP
cuttree_R(SEXP tree_R, SEXP num_of_clust_R) 
{
	size_t *labels, num_of_clust, num_of_points;
	SEXP labels_R;

	/* Dimensions */
	num_of_points = nrows(tree_R) + 1.; // It takes n steps to merge n+1 points 
	num_of_clust = (size_t) INTEGER(num_of_clust_R)[0];
	
	/* Allocate memory */
	labels_R = PROTECT(allocVector(INTSXP, num_of_points));
	labels = (size_t*)R_alloc(sizeof(*labels), 2*num_of_points);

	/* Cut tree */
	cutree(labels, REAL(tree_R), num_of_clust, num_of_points); 

	/* Copy into output vector */
	for(size_t i = 0; i < num_of_points; ++i) { 
		INTEGER(labels_R)[i] = (int)labels[i] + 1; // C labels start from 0, R labels from 1
	}

	UNPROTECT(1);
	return labels_R;
}











