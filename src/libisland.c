
#include "libisland.h"
#include <assert.h>
#include <math.h>
#include "heap.h"
#include <stdio.h>
#include <stdlib.h>



int
add(int a, int b) {
	return a + b;
}

int
subtract(int a, int b) {
	return a - b;
}






double
bino_div(double x, double y)
{
	if (y > 0.)
		return x / y;
	else
		return 0.5;
}

double
bino_div_id(int x, double y)
{
	if (!y == 0)
		return (double)x / y;
	else
		return 0.5;
}



double
bino_xlogy(double x, double y)
{
	if (y > 0.)
		return x * log(y);
	else
		return 0.;
}

double
bino_xlogy_id(int x, double y)
{
	if (!(y == 0))
		return x * log(y);
	else
		return 0;
}

double
cost(double k0, double k1)
{

	// Return 0 if k0 and or k1 are NaN
	if(k0 < 0) {
		return 0;
	}

	if(k1 < 0) {
		return 0;
	}

	// Else compute likelihood
	double p;
	p = bino_div(k1, k0+k1);

	return (-(bino_xlogy(k1, p) + bino_xlogy(k0, 1.-p)));
}


double
nLLR(const data_elem_t *K1, const data_elem_t *K2, double LOG_L1, double LOG_L2, size_t n) {
	size_t j;
	double res = 0;


/* #pragma omp parallel for schedule(static) reduction(+: res) */
	for( j = 0; j < n; ++j) {

		// Replace NaN with 0
		// res += cost(fmax(K1[j].k0,0) + fmax(K2[j].k0,0), fmax(K1[j].k1,0) + fmax(K2[j].k1,0)) - (cost(K1[j].k0, K1[j].k1) + cost(K2[j].k0, K2[j].k1));
		res += cost(fmax(K1[j].k0,0) + fmax(K2[j].k0,0), fmax(K1[j].k1,0) + fmax(K2[j].k1,0)) ;
	}

	return res - (LOG_L1 + LOG_L2) ;
}

#define LOG_PI (1.1447298858494001638774761886452324688434600830078125)

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
merged_cost(const data_elem_t *K1, const data_elem_t *K2, double LOG_L1, double LOG_L2, const int chr1, const int chr2, const int pos1, const int pos2, size_t n) {
	distance_t D;

	D.nllr = nLLR(K1, K2, LOG_L1, LOG_L2, n);  

	D.d = dist(chr1, chr2, pos1, pos2);
	return D;
}


void
init_state(distance_t *D, double *L, size_t *PL, size_t *PR, const data_elem_t *K, double *LOG_L, double *tracked_cost, const int *chr, const int *pos, size_t m, size_t n)
{
	size_t I_NIL, i, j;

	assert(m > 0);

	/* Get invalid index */
	I_NIL = m;

	for(i = 0; i<m; ++i) {
		LOG_L[i] = 0;
		for (j = 0; j < n; ++j ) {
			LOG_L[i] += cost(K[n*i + j].k0, K[n*i + j].k1);
		}
	}
	

	/* Set up valid data */
	for (i = 0; i < m-1; ++i) {
		D[i] = merged_cost(&K[n*i], &K[n*(i+1)], LOG_L[i], LOG_L[i+1], chr[i], chr[i+1], pos[i], pos[i+1], n);
		PL[i+1] = i;
		PR[i] = i+1;
	}

	/* SEt up sentinel */
	{
		D[m-1].d = INFINITY;
		D[m-1].nllr = INFINITY;
		PL[0] = I_NIL;
		PR[m-1] = I_NIL;

		D[m].d = INFINITY;
		D[m].nllr = INFINITY; /* nil slot */
	}
	

	/* Set up labels */
	for (i = 0; i < m; ++i) {
		L[i] = -(double)(i+1);
		tracked_cost[i] = 0;
	}
		
}

size_t 
new_pair(size_t p, distance_t *D, size_t *PL, size_t *PR, size_t m)
{

	double v, w;

	(void)m;
	if ( (D[PL[p]].d + D[PL[p]].nllr) < (D[PR[p]].d + D[PR[p]].nllr) ) {
		/* Scan left */
		v = D[p].d + D[p].nllr;
		while ( (w = (D[PL[p]].d + D[PL[p]].nllr )) < v) {
			p = PL[p];
			v = w;
		}
	} else {
		/* Scan right */ 
		v = D[p].d + D[p].nllr;
		while ( (w = (D[PR[p]].d + D[PR[p]].nllr)) < v) {
			p = PR[p];
			v = w;
		}
	}
	return p;
}

double
max_d(double A, double B) {
	if(A > B) {
		return A;
	}
	if(B > A) {
		return B;
	}
	return A;
}


void
update_K(data_elem_t *K, double *LOG_L, size_t p1, size_t p2, size_t n) {
	size_t j;

	for (j = 0; j < n; j++) {
		double k0, k1, l0, l1;

		/* NaNs are transformed to -2147483648. So, if a k-value is <0, just add 0. Otherwise, add the value itself*/
		
		k0 = max_d(K[n*p1 + j].k0, 0);
		k1 = max_d(K[n*p1 + j].k1, 0);
		l0 = max_d(K[n*p2 + j].k0, 0);
		l1 = max_d(K[n*p2 + j].k1, 0);


		K[n*p1 + j].k0 = k0 + l0;
		K[n*p1 + j].k1 = k1 + l1;
	}
	LOG_L[p1] = LOG_L[p1] + LOG_L[p2];
}




#define double_LESS(x, y, H)  ( (H)[((x))] < (H)[((y))] )

struct cookie {
	double *H;
	const double *left;
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
sort_tree_size(size_t m)
{
    return m + m + m + m;
}

#include <stdio.h>
#include <alloca.h>
#include <stdlib.h>


void
dump_tree(const double *Z, size_t m1)
{
	fprintf(stderr, "%s\t%s\t%s\n", "merge.1", "merge.2", "height");
	for (size_t i = 0; i < m1; ++i)
		fprintf(stderr, "%f" "\t" "%f" "\t" "%g" "\n", Z[0*m1+i], Z[1*m1+i], Z[2*m1+i]);
}

void
sort_tree(double *Z_new, void *scratch, const double *Z, size_t m)
{
	size_t *perm, r, size, p, j, *inv_perm, *heap;
	double q1, q2, *left;
	struct cookie cookie;

	// int * used ;
	// used = (int *)alloca( m * sizeof(*used));
	// {size_t i;
	// for (i = 0 ; i < m ; ++i)
	// 	used[i] = 0;
	// }

	/* Allocate memory */
    left = (double *)scratch;
    perm = (size_t *)&left[m];
    inv_perm = (size_t *)&perm[m];
    heap = (size_t *)&inv_perm[m];
	
	/* Get the leftmost data point in every cluster */
	for( j = 0; j < m; ++j) {
		left[j] = Z[0*m + j];

		if (left[j] > 0) {
			left[j] = left[(size_t)left[j]-1];
		} 
	}

	/* Height */
	double *H;
	H = malloc(sizeof(*H) * (m));

	for( j=0; j<m; ++j) {
		H[j] = Z[2*m + j] + Z[3*m + j];
	}

	cookie.H = &H[0];
	cookie.left = left;
	

	size = 0;
	
	/* Put root to heap */
	*heap = m-1;
	++size;
// used[heap[size]]++;
	heap_size_t_make_heap(heap, size, &cookie);
	

	for ( r = m ; r-- > 0; ) {    

		// if (size == 0)
		// {
		// 		fprintf(stderr, "FALIL\n");
		// 	for (size_t i = 0; i < m ; ++i) {
		// 		if (used[i] != 1) {
		// 			fprintf(stderr, "used[%zu] = %d\n", i, used[i]);

		// 		dump_tree(Z, m);
		// 			abort();
		// 		}
		// 	}

		// }

		/* pop best from heap */
		p = heap[0];
	   	heap_size_t_pop_heap( heap, size, &cookie);
	    --size;
		
		/* assign location */
		perm[r] = p;

		/* get children */
		q1 = Z[0*m + p];
		q2 = Z[1*m + p];
		
		/* if not leaf, add to stack */
		if (q1 > 0) {
			heap[size] = (size_t)q1 - 1;
// used[heap[size]]++;
			++size;
	 	 	heap_size_t_push_heap(heap, size, &cookie);  /* O(log n ) */
		}
			
		if (q2 > 0) {
			heap[size] = (size_t)q2 - 1;
// used[heap[size]]++;
			++size;
	 	 	heap_size_t_push_heap(heap, size, &cookie);  /* O(log n ) */
		}
	}

	for( j = 0; j < m; j++) {
		inv_perm[perm[j]] = j;

		Z_new[0*m + j] = Z[0*m + perm[j]];
		Z_new[1*m + j] = Z[1*m + perm[j]];
		Z_new[2*m + j] = Z[2*m + perm[j]];
		Z_new[3*m + j] = Z[3*m + perm[j]];
	}
	
	for ( j = 0; j < m; ++j) {
		if (Z_new[0*m + j] > 0) {
			Z_new[0*m + j] = (double)(inv_perm[(size_t)Z_new[0*m + j] - 1] + 1);
		}
		if (Z_new[1*m + j] > 0) {
			Z_new[1*m + j] = (double)(inv_perm[(size_t)Z_new[1*m + j] - 1] + 1);
		}
	}
}

size_t 
island_cluster_size(size_t m, size_t n)
{
    return m+1 + m + 2*m*n + m + m;
}

double
logl_for_clust(size_t *clust, size_t value, double *LOG_L, size_t m)
{
	size_t i;
	double res;

	res = 0;

	for(i = 0; i < m; ++i) {
		if (clust[i] == value) {
			res += LOG_L[i];
		}
	}

	return res;
}

// void
// island_cluster(void *scratch, double *Z, const int *K0, const int *K1, const int *chr, const int *pos, size_t m, size_t n) <-- Handling memory allocations inside function
void
island_cluster(double *Z, const int *K0, const int *K1, const int *chr, const int *pos, size_t m, size_t n)
{
	/* Function that takes as input K0 and K1 and performs a hierarchical clustering. Output tree is unsorted, needs to be sorted with sort_tree().
	
	Input: 
	double *Z: vector for storing results, 
	const int* K0, K1: input matrices with count values (samples X CpG-sites), 
	const int* chr: vector telling which chromosome each CpG-site belongs to, 
	const int *pos: vector telling the genomic coordinate of each CpG-site
	size_t m, n: number of CpG-sites (rows) and samples (cols) in K0 and K1
	
	Output:
	double Z*: vector holding results. Every row is a merge, first two columns holds identifiers of merged points, third column holds value of log-likelihood function for the merge.
	*/




	/* Initializing helper variables
	double *D: vector for storing distances between CpG-sites
	double *L: vector for storing labels of CpG-sites or clusters. L < 0 means single CpG site, L > 0 means cluster. 
	size_t *PL, *PR: vectors that points to the index of the left or right neighbour (these change as points merge into clusters).
	size_t i, j: used for indexing 
	size_t I_NO: invalid index
	size_t p1, p2: indices of points currently being clustered
	data_elem_t K: matrix holdinf the values of K0 and K1 at every CpG-site for every sample as doubles.*/
	double *L, *LOG_L, *tracked_cost;
	size_t *PL, *PR, i, j, I_NO, p1, p2;
	data_elem_t *K;
	distance_t *D;

	/* Check empty problem */
	if (!(m > 0) || !(n > 0))
		return;

	/* Invalid index */
	I_NO = m; 

	/* Allocating memory */
	D = malloc(sizeof(*D) * (m+1));
	L = malloc(sizeof(*L) * m);
	K = malloc(sizeof(*K) * m*n);
	PL = malloc(sizeof(*PL) * m);
	PR = malloc(sizeof(*PR) * m);
	LOG_L = malloc(sizeof(*LOG_L) * m);
	tracked_cost = malloc(sizeof(*tracked_cost) * m);

	// TODO:
	// 	- assert pointers are not NULL
	// 	- assert stack is not overflowing/not using too much memory etc


	/* R matrices need to be transposed for indexing to be correct. Transposing K0 and K1 and copying to K. Also avoid manipulation of R object within function. */
	for (j = 0; j < n; ++j) {
		for (i = 0; i < m; ++i) {
			K[ i*n + j ].k0 = K0[ i + j*m ]; 
			K[ i*n + j ].k1 = K1[ i + j*m ];
		}
	}
	
	/* Initializing rest of variables */
	init_state( D, L, PL, PR, K, LOG_L, tracked_cost, chr, pos, m, n );
	
	/* Boot current pointer */
	p1 = 0; 

	/* Cluster loop */
	for (j = 0; j < m-1; ++j) {


		//printf("%lu / %lu\n", j, m-2);
		/* Scanning to find new pair p1 and p2 to merge. p2 is to the right of p1*/
		p1 = new_pair(p1, D, PL, PR, m);
		p2 = PR[p1];

		/* Record merge into Z */
		Z[0*(m-1) + j] = L[p1];
		Z[1*(m-1) + j] = L[p2];
		Z[2*(m-1) + j] = D[p1].nllr - tracked_cost[p1] - tracked_cost[p2];
		Z[3*(m-1) + j] = D[p1].d;

		/* Update tracked cost*/
		tracked_cost[p1] = D[p1].nllr;

		/* Update labels. 
		L[p1] is j+1, that is, the row of the merge (R indexing starting from 1 and not 0). 
		L[p2] is 0 as it can no longer be merged with its left neighbour. */
		L[p1] = j+1;
#ifndef NDEBUG
		L[p2] = 0;
#endif

		/* Update pointers.
		Left and right neighbours of new cluster are assigned to p1. Also the neighbours neighbour information needs updating, unless new right neighbour of cluster is end of vector (invalid index I_NO). 
		p2 is used and assigned to be */
		PR[p1] = PR[p2];
		if (PR[p1] != I_NO)
			PL[PR[p1]] = p1;
#if 1
		PL[p2] = I_NO;
		PR[p2] = I_NO;
#endif

		/* Update data.
		K now holds the sum of all merged points */
		update_K(K, LOG_L, p1, p2, n);

		/* Update distances */
		if ( PL[p1] != I_NO) {
			D[PL[p1]] = merged_cost( &K[n*PL[p1]], &K[n*p1], LOG_L[PL[p1]], LOG_L[p1], chr[PL[p1]], chr[p1], pos[p1-1], pos[p1], n);
			
		}
		if ( PR[p1] != I_NO) {
			D[p1] = merged_cost(&K[n*p1], &K[n*PR[p1]], LOG_L[p1], LOG_L[PR[p1]], chr[p1], chr[PR[p1]], pos[PR[p1]-1], pos[PR[p1]], n);
			

		} else {
			D[p1].d = INFINITY; // TODO: Smoother way to do this?
			D[p1].nllr = INFINITY;
		}
#ifndef NDEBUG
			D[p2].d = INFINITY;
			D[p2].nllr = INFINITY;
#endif

		if (PR[p1] == I_NO)
			p1 = PL[p1];
	}


	/* Freeing memory */
	free(D);
	free(L);
	free(K);
	free(PL);
	free(PR);
	free(LOG_L);
	free(tracked_cost);

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

