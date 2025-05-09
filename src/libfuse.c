/**
 * @file libfuse.c
 * @brief Implementation of binomial clustering, sorting, and cutting of tree.
 */

#include "libfuse.h"
#include <assert.h>
#include <math.h>
#include "heap.h"
#include <stdio.h>
#include <stdlib.h>
#include <alloca.h>


#define LOG_PI (1.1447298858494001638774761886452324688434600830078125)

double
bino_div(double x, double y) {
    return (y > 0.) ? x / y : 0.5;
}

double
bino_div_id(int x, double y) {
    return (y != 0) ? (double)x / y : 0.5;
}

double
bino_xlogy(double x, double y) {
    return (y > 0.) ? x * log(y) : 0.;
}

double
bino_xlogy_id(int x, double y) {
    return (y != 0) ? x * log(y) : 0;
}

double
cost(double k0, double k1)
{

	// Return 0 if k0 and or k1 are NaN
	if (k0 < 0 || k1 < 0) {
        return 0;
    }

	// Else compute likelihood
	double p;
	p = bino_div(k1, k0+k1);

	return (-(bino_xlogy(k1, p) + bino_xlogy(k0, 1.-p)));
}

double
nLLR(const data_elem_t *LHC, const data_elem_t *RHC, double LOGL_LHC, double LOGL_RHC, size_t num_of_samples) {
	size_t j;
	double logl = 0;

	// Looping trough LHC and RHC to compute combined likelihoood
	for( j = 0; j < num_of_samples; ++j) {
		logl += cost((is_na_double(LHC[j].k0) ? 0 : LHC[j].k0) + (is_na_double(RHC[j].k0) ? 0 : RHC[j].k0),
    				(is_na_double(LHC[j].k1) ? 0 : LHC[j].k1) + (is_na_double(RHC[j].k1) ? 0 : RHC[j].k1));


	}	

	return logl - (LOGL_LHC + LOGL_RHC) ;
}

static
double 
dist(const int chr1, const int chr2, const int pos1, const int pos2) {
	double x;

	if(chr1 != chr2) 
		return INFINITY;

	x = 0.01 * (pos2 - pos1);
	return log( 1 + x*x ) + LOG_PI;
}

distance_t 
merged_cost(const data_elem_t *LHC, const data_elem_t *RHC, double LOGL_LHC, double LOGL_RHC, const int chrL, const int chrR, const int posL, const int posR, size_t num_of_samples) {
	distance_t distance;
	distance.nllr = nLLR(LHC, RHC, LOGL_LHC, LOGL_RHC, num_of_samples);  
	distance.d = dist(chrL, chrR , posL, posR);
	return distance;
}

void
init_state(distance_t *distance, double *labels, size_t *l_neighbor, size_t *r_neighbor, const data_elem_t *counts, double *logl, double *tracked_cost, const int *chr, const int *pos, size_t num_of_sites, size_t num_of_samples)
{
	size_t I_NIL, i, j;

	/* Assert data is not empty */
	assert(num_of_sites > 0);

	/* Get invalid index */
	I_NIL = num_of_sites;

	/* Compute binomial log-likelihoods for each site */
	for(i = 0; i<num_of_sites; ++i) {
		logl[i] = 0;
		for (j = 0; j < num_of_samples; ++j ) {
			logl[i] += cost((is_na_double(counts[num_of_samples*i + j].k0) ? 0 : counts[num_of_samples*i + j].k0), 
							(is_na_double(counts[num_of_samples*i + j].k1) ? 0 : counts[num_of_samples*i + j].k1));

			
		}
	}

	/* Set up valid distances, fill in left and right neighbors */
	for (i = 0; i < num_of_sites-1; ++i) {
		distance[i] = merged_cost(&counts[num_of_samples*i], &counts[num_of_samples*(i+1)], logl[i], logl[i+1], chr[i], chr[i+1], pos[i], pos[i+1], num_of_samples);
		l_neighbor[i+1] = i;
		r_neighbor[i] = i+1;
	}

	/* Set up extremes */
	{
		/* Nothing right of the rightmost point */
		distance[num_of_sites-1].d = INFINITY; 
		distance[num_of_sites-1].nllr = INFINITY;

		/* Extremes don't have neighbors */
		l_neighbor[0] = I_NIL;
		r_neighbor[num_of_sites-1] = I_NIL;

		/* Placeholder slot to the furthest right */
		distance[num_of_sites].d = INFINITY;
		distance[num_of_sites].nllr = INFINITY; /* nil slot */
	}
	
	/* Set up labels */
	for (i = 0; i < num_of_sites; ++i) {
		labels[i] = -(double)(i+1);
		tracked_cost[i] = 0;
	}
		
}

size_t 
new_pair(size_t current, distance_t *distance, size_t *l_neighbor, size_t *r_neighbor)
{
	double current_pair, best_pair;

	/* Do we go left or right?*/
	if ( (distance[l_neighbor[current]].d + distance[l_neighbor[current]].nllr) < (distance[r_neighbor[current]].d + distance[r_neighbor[current]].nllr) ) {
		/* Scan left */
		current_pair = distance[current].d + distance[current].nllr;
		while ( (best_pair = (distance[l_neighbor[current]].d + distance[l_neighbor[current]].nllr )) < current_pair) {
			current = l_neighbor[current];
			current_pair = best_pair;
		}
	} else {
		/* Scan right */ 
		current_pair = distance[current].d + distance[current].nllr;
		while ( (best_pair = (distance[r_neighbor[current]].d + distance[r_neighbor[current]].nllr)) < current_pair) {
			current = r_neighbor[current];
			current_pair = best_pair;
		}
	}
	return current;
}

void
update_counts(data_elem_t *counts, double *logl, size_t LHC, size_t RHC, size_t num_of_samples) {
	size_t j;

	for (j = 0; j < num_of_samples; j++) {
		double k0, k1, l0, l1;

		/* Get counts for LHC and RHC */
		k0 = is_na_double(counts[num_of_samples * LHC + j].k0) ? 0 : counts[num_of_samples * LHC + j].k0; // Unmethylated counts in LHC
		k1 = is_na_double(counts[num_of_samples * LHC + j].k1) ? 0 : counts[num_of_samples * LHC + j].k1; // Methylated counts in LHC
		l0 = is_na_double(counts[num_of_samples * RHC + j].k0) ? 0 : counts[num_of_samples * RHC + j].k0; // Unmethylated counts in RHC
		l1 = is_na_double(counts[num_of_samples * RHC + j].k1) ? 0 : counts[num_of_samples * RHC + j].k1; // Methylated coutns in RHC


		/* New counts are written into index of the lefthand cluster */
		counts[num_of_samples*LHC + j].k0 = k0 + l0; 
		counts[num_of_samples*LHC + j].k1 = k1 + l1;
	}

	/* Sum of likelihoods for the cluster is updated*/
	logl[LHC] = logl[LHC] + logl[RHC];
}

#define double_LESS(x, y, H)  ( (H)[((x))] < (H)[((y))] )

struct cookie {
	double *H; // Height
	const double *left; // Label of leftmost point in cluster
};

static
int 
larger(size_t x, size_t y, const struct cookie *cookie)
{
	if (cookie->H[x] < cookie->H[y])
		return 1;
	if (cookie->H[x] > cookie->H[y])
		return 0;

	/* Breaking a tie */
	if (cookie->left[x] > cookie->left[y])
		return 1;
	
	return 0;
}

HEAP_DEFINE(heap_size_t, size_t, size_t, larger, const struct cookie *)

size_t 
sort_tree_size(size_t nrow_tree)
{
    return nrow_tree + nrow_tree + nrow_tree + nrow_tree + nrow_tree;
}

void
sort_tree(double *tree_sorted, void *scratch, const double *tree_unsorted, size_t nrow_tree)
{
	size_t *perm, r, size, p, j, *inv_perm, *heap;
	double q1, q2, *left;
	struct cookie cookie;

	/* Allocate memory */
    left = (double *)scratch; // Leftmost point in every cluster
    perm = (size_t *)&left[nrow_tree]; // id of tree_sorted in tree_unsorted
    inv_perm = (size_t *)&perm[nrow_tree]; // id of tree_unsorted in tree_sorted
    heap = (size_t *)&inv_perm[nrow_tree]; // heap for storing points under consideration
	
	/* Get the leftmost data point in every cluster from column 0 */
	for( j = 0; j < nrow_tree; ++j) {
		left[j] = tree_unsorted[0*nrow_tree + j];

		/* If labels are >0, they are row indices. Starting from 1 in R and from 0 in C, so adding -1*/
		if (left[j] > 0) {
			left[j] = left[(size_t)left[j]-1];  
		} 
	}

	/* Height */
	double *H;
	H = malloc(sizeof(*H) * (nrow_tree));

	/* Sorting based on sum of columns 2 and 3*/
	for( j=0; j<nrow_tree; ++j) {
		H[j] = tree_unsorted[3*nrow_tree + j] + tree_unsorted[4*nrow_tree + j];
	}

	cookie.H = &H[0];
	cookie.left = left;
	
	size = 0; // Initial size of heap is 0
	
	/* Put root to heap */
	*heap = nrow_tree-1;
	++size;
	heap_size_t_make_heap(heap, size, &cookie);
	
	for ( r = nrow_tree ; r-- > 0; ) {  
		
		/* pop best from heap */
		p = heap[0];
	   	heap_size_t_pop_heap( heap, size, &cookie);
	    --size;
		
		/* assign location */
		perm[r] = p;

		/* get children */
		q1 = tree_unsorted[0*nrow_tree + p];
		q2 = tree_unsorted[1*nrow_tree + p];
		
		/* if not leaf, add to heap */
		if (q1 > 0) {
			heap[size] = (size_t)q1 - 1;
			++size;
	 	 	heap_size_t_push_heap(heap, size, &cookie);  /* O(log n ) */
		}
			
		if (q2 > 0) {
			heap[size] = (size_t)q2 - 1;
			++size;
	 	 	heap_size_t_push_heap(heap, size, &cookie);  /* O(log n ) */
		}
	}

	/* Filling sorted tree */
	for( j = 0; j < nrow_tree; j++) {
		inv_perm[perm[j]] = j;

		tree_sorted[0*nrow_tree + j] = tree_unsorted[0*nrow_tree + perm[j]];
		tree_sorted[1*nrow_tree + j] = tree_unsorted[1*nrow_tree + perm[j]];
		tree_sorted[2*nrow_tree + j] = tree_unsorted[2*nrow_tree + perm[j]];
		tree_sorted[3*nrow_tree + j] = tree_unsorted[3*nrow_tree + perm[j]];
		tree_sorted[4*nrow_tree + j] = tree_unsorted[4*nrow_tree + perm[j]];
	}
	
	/* If labels are positive they indicate which row the merge happended on. This changes as well. In R from 1, C from 0. */
	for ( j = 0; j < nrow_tree; ++j) {
		if (tree_sorted[0*nrow_tree + j] > 0) {
			tree_sorted[0*nrow_tree + j] = (double)(inv_perm[(size_t)tree_sorted[0*nrow_tree + j] - 1] + 1);
		}
		if (tree_sorted[1*nrow_tree + j] > 0) {
			tree_sorted[1*nrow_tree + j] = (double)(inv_perm[(size_t)tree_sorted[1*nrow_tree + j] - 1] + 1);
		}
	}
}

void
fuse_cluster(double *tree, data_elem_t *counts, const int *chr, const int *pos, size_t num_of_sites, size_t num_of_samples)
{
	/* Initializing helper variables */
	distance_t *distance; // vector for storing distances between CpG-sites
	double *labels; // vector for storing labels of CpG-sites or clusters. L < 0 means single CpG site, L > 0 means cluster and points to row
	double *logl; // vector holding total binomial log-likelihood of each cluster
	double *tracked_cost; // vector holding changes in total likelihood of data when each cluster was formed
	double new_clust_logl; // For holding log-likelihood of new cluster
	size_t *l_neighbor, *r_neighbor; // vectors that point to the index of the left or right neighbor (these change as points merge into clusters).
	size_t i, j; // counters
	size_t I_NO; // invalid index (rightmost)
	size_t LHC, RHC; // Lefthand cluster and righthand cluster to merge
	

	/* Check empty problem */
	if (!(num_of_sites > 0) || !(num_of_samples > 0))
		return;

	/* Invalid index is rightmost */
	I_NO = num_of_sites; 

	/* Allocating memory */
	distance = malloc(sizeof(*distance) * (num_of_sites+1));
	labels = malloc(sizeof(*labels) * num_of_sites);
	logl = malloc(sizeof(*logl) * num_of_sites);
	l_neighbor = malloc(sizeof(*l_neighbor) * num_of_sites);
	r_neighbor = malloc(sizeof(*r_neighbor) * num_of_sites);
	tracked_cost = malloc(sizeof(*tracked_cost) * num_of_sites);
	
	/* Initializing rest of variables */
	init_state( distance, labels, l_neighbor, r_neighbor, counts, logl, tracked_cost, chr, pos, num_of_sites, num_of_samples );
	
	/* Starting from left */
	LHC = 0; 

	/* Cluster loop */
	for (j = 0; j < num_of_sites-1; ++j) {

		// Resetting
		new_clust_logl = 0;

		/* Scanning to find new pair LHC and RHC to merge. RHC is to the right of LHC*/
		LHC = new_pair(LHC, distance, l_neighbor, r_neighbor);
		RHC = r_neighbor[LHC];

		/* Compute the total sum of likelihoods of merge and add to tree */ 
		for (i = 0; i < num_of_samples; ++i ) {
			 new_clust_logl += cost((is_na_double(counts[num_of_samples*LHC + i].k0) ? 0 :counts[num_of_samples*LHC + i].k0) + 
			 						(is_na_double(counts[num_of_samples*RHC + i].k0) ? 0 :counts[num_of_samples*RHC + i].k0), 
									(is_na_double(counts[num_of_samples*LHC + i].k1) ? 0 :counts[num_of_samples*LHC + i].k1) + 
									(is_na_double(counts[num_of_samples*RHC + i].k1) ? 0 :counts[num_of_samples*RHC + i].k1));	
		}

		/* Record merge into tree */
		tree[0*(num_of_sites-1) + j] = labels[LHC]; // label of first merged point
		tree[1*(num_of_sites-1) + j] = labels[RHC]; // label of second merged point
		tree[2*(num_of_sites-1) + j] = distance[LHC].nllr - tracked_cost[LHC] - tracked_cost[RHC]; // change in total likelihood of data
		tree[3*(num_of_sites-1) + j] = new_clust_logl - logl[LHC] - logl[RHC]; // total likelihood of all points in merge
		tree[4*(num_of_sites-1) + j] = distance[LHC].d; // genomic distance penalty

		/* Update tracked cost */
		tracked_cost[LHC] = distance[LHC].nllr;

		/* Update labels. 
		L[LHC] is j+1, that is, the row of the merge (R indexing starting from 1 and not 0). 
		L[RHC] is 0 as it can no longer be merged with its left neighbour. */
		labels[LHC] = j+1;
		labels[RHC] = 0;

		/* Update pointers.
		Left and right neighbors of new cluster are assigned to LHC. 
		Also the neighbours neighbor information needs updating, unless new right neighbor of cluster is end of vector (invalid index I_NO). 
		RHC is used and assigned to be invalid.*/
		r_neighbor[LHC] = r_neighbor[RHC];
		if (r_neighbor[LHC] != I_NO)
			l_neighbor[r_neighbor[LHC]] = LHC;

		l_neighbor[RHC] = I_NO;
		r_neighbor[RHC] = I_NO;

		/* Update data: counts[LHC] now holds the sum of all merged points */
		update_counts(counts, logl, LHC, RHC, num_of_samples);

		/* Update distances */
		if ( l_neighbor[LHC] != I_NO) {
			distance[l_neighbor[LHC]] = merged_cost( &counts[num_of_samples*l_neighbor[LHC]], &counts[num_of_samples*LHC], logl[l_neighbor[LHC]], logl[LHC], chr[l_neighbor[LHC]], chr[LHC], pos[LHC-1], pos[LHC], num_of_samples);
			
		}
		if ( r_neighbor[LHC] != I_NO) {
			distance[LHC] = merged_cost(&counts[num_of_samples*LHC], &counts[num_of_samples*r_neighbor[LHC]], logl[LHC], logl[r_neighbor[LHC]], chr[LHC], chr[r_neighbor[LHC]], pos[r_neighbor[LHC]-1], pos[r_neighbor[LHC]], num_of_samples);
			

		} else {
			distance[LHC].d = INFINITY; 
			distance[LHC].nllr = INFINITY;
		}
			distance[RHC].d = INFINITY;
			distance[RHC].nllr = INFINITY;

		if (r_neighbor[LHC] == I_NO)
			LHC = l_neighbor[LHC];
	}


	/* Freeing memory */
	free(distance);
	free(labels);
	free(logl);
	free(tracked_cost);
	free(l_neighbor);
	free(r_neighbor);
}

static
size_t
rlabel2linear(double L, size_t N)
{
	if (L < 0.)
		/* -1 => 0 , -2 => 1, ... -N => N-1 */
		return (size_t)( -L - 1. );
	else
		/* 1 => n, 2 => n+1 , ... N => 2*N-2 */
		return (size_t)( L - 1. ) + N;
}

// TODO
void
cutree(size_t *L, const double *Z, size_t k, size_t n)
{
	size_t i, seen;

	assert(L != NULL || n < 1);
	assert(Z != NULL || n < 2);
	assert(1 <= k && k <= n);

	/* Label each leaf */ // OK!
	for (i = 0; i < n; ++i) {
		L[i] = i;
	}
		

	/* Label each using the first member */ // OK!
	for (i = 0; i < n-k; ++i)
	{
		L[n+i] = L[ rlabel2linear( Z[i + 0*(n-1)], n) ];
		if (L[ rlabel2linear( Z[i + 1*(n-1)], n)] < L[n+i])
		 	L[n+i] = L[ rlabel2linear(Z[i + 1*(n-1)], n) ];
		
		
	}

	/* Propagate labels down the tree from the cut */ // NOT TESTED
	for (i = n-k; i-- > 0;)
	{
		L[ rlabel2linear(Z[i + 0*(n-1)], n) ] = L[n+i];
		L[ rlabel2linear( Z[i + 1*(n-1)], n) ] = L[n+i];
	}

	/* Compress labels from [0,n) to [0,k) */ // NOT TESTED
	seen = 0;
	for (i = 0; i < n; ++i)
		if (L[i] == i)
			L[i] = seen++;
		else
			L[i] = L[L[i]];
}


void
cor_pearson_rowSums( double *sxy, double *sy1, double *sy2, const double *x, const double *Y, size_t m, size_t ny)
{
	size_t i, j;
#pragma omp parallel for schedule(static)
  for (j = 0 ; j < ny ; ++j)
  {
	double s, t, u;

    s = 0.;
#pragma omp simd reduction(+: s)
    for (i = 0; i < m; ++i)
      s += Y[i + j*m];
      
    t = 0.;
#pragma omp simd reduction(+: t)
    for (i = 0; i < m; ++i)
      t += Y[i + j*m] * Y[i + j*m];
      
    u = 0.;
#pragma omp simd reduction(+: u)
    for (i = 0; i < m; ++i)
      u += x[i] * Y[i + j*m];

	sxy[j] = u;
    sy1[j] = s;
    sy2[j] = t;
  }    
}

