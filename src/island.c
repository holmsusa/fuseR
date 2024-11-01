
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <assert.h>
#include "libisland.h"  // Your additional custom header



static
__attribute__((unused))
int
is2D(SEXP arg)
{ return isVector(arg) || isMatrix(arg); }



SEXP
sort_tree_R(SEXP Z_R)
{
	SEXP Z_new_R ; 
	size_t m;
	double *scratch;

	Z_new_R = PROTECT(duplicate(Z_R));
	m = nrows(Z_R);

	scratch = (double*)R_alloc(sizeof(*scratch), sort_tree_size(m));

	sort_tree( REAL(Z_new_R), scratch, REAL(Z_R), m );

	//return Z_new_R;

	UNPROTECT(1);
	return Z_new_R;
}

SEXP 
island_cluster_R(SEXP K0_R, SEXP K1_R, SEXP CHR_R, SEXP POS_R)
{

	
	int m, n;
	SEXP Z_R;

	/* check type for input */
	if (!(isInteger(K0_R)))
		error("'%s' must be an integer array", "K0");
	if (!(isInteger(K1_R)))
		error("'%s' must be an integer array", "K1");


	/* Get dimensions */
	m = nrows(K0_R);
	n = ncols(K0_R);


	/* Check dimensions */
	if (!(nrows(K0_R) == m && ncols(K0_R) == n))
		error("'%s' must be %d-by-%d", "K0", m, n);
	if (!(nrows(K1_R) == m && ncols(K1_R) == n))
		error("'%s' must be %d-by-%d", "K1", m, n);

	/* Allocate memory */
	Z_R = PROTECT(allocMatrix(REALSXP, m-1, 4));

	/* Run clustering */
	island_cluster( REAL(Z_R), INTEGER(K0_R), INTEGER(K1_R), INTEGER(CHR_R), INTEGER(POS_R), m, n );

	UNPROTECT(1);
	return Z_R;
}

SEXP
cuttree_R(SEXP Z_R, SEXP K_R) 
{
	size_t *L, k, m;

	SEXP L_R;

	/* Dimensions */
	m = nrows(Z_R) + 1.;
	k = (size_t) INTEGER(K_R)[0];

 
	/* Allocate size_t *L and SEXP or int *L. Compute size_t L with function, return int L + 1 since R starts from 1 not 0 */
	
	/* Allocate memory */
	L_R = PROTECT(allocVector(INTSXP, m));
	L = (size_t*)R_alloc(sizeof(*L), 2*m);


	/* Cut tree */
	// L: cluster labels, Z: tree, k: number of clusters, m: number of sites
	cutree(L, REAL(Z_R), k, m); // <-- FAILS HERE

	/* Copy into other vector */
	for(size_t i = 0; i < m; ++i) { // for(size_t i = 0; i <=m; ++i) {
		INTEGER(L_R)[i] = (int)L[i] + 1; // <-- Fails here
	}

	UNPROTECT(1);
	return L_R;
	
}

SEXP
cor_pearson_rowSums_R(SEXP x_R, SEXP Y_R) 
{
	int m, nx, ny;
	SEXP sxy_R, sy1_R, sy2_R, res_R;

	/* Dimensions */
	m = nrows(x_R); // Number of samples
	nx = ncols(x_R); // NUmber of segments in X
	ny = ncols(Y_R); // Number of segments in Y

	if (!(nrows(x_R) == m && ncols(x_R) == 1))
		error("x must be vector");
	
	if (!(nrows(Y_R) == m && ncols(Y_R) == ny))
		error("y must have same dimension as x");

	if (!isReal(x_R))
		error("x must be doubles");
	if (!isReal(Y_R))
		error("y must be doubles");

	(void)nx;
	
	/* Allocate memory for output */
	sxy_R = PROTECT(allocVector(REALSXP, ny));
	sy1_R = PROTECT(allocVector(REALSXP, ny));
	sy2_R = PROTECT(allocVector(REALSXP, ny));
	


	/* Correlation */
	cor_pearson_rowSums(REAL(sxy_R), REAL(sy1_R), REAL(sy2_R), REAL(x_R), REAL(Y_R), m, ny) ;


	/* Create output tuple */
	res_R = PROTECT(allocList(3));
	SETCAR(res_R, sy1_R);
	SET_TAG(res_R, install("sy1"));
	SETCAR(CDR(res_R), sy2_R);
	SET_TAG(CDR(res_R), install("sy2"));
	SETCAR(CDDR(res_R), sxy_R);
	SET_TAG(CDDR(res_R), install("sxy"));

	UNPROTECT(4);
	return res_R;
}


// Register routines
static const R_CallMethodDef callMethods[] = {
  {"island_cluster_R", (DL_FUNC) &island_cluster_R, 4},
  {"sort_tree_R",      (DL_FUNC) &sort_tree_R,      1},
  {"cuttree_R",        (DL_FUNC) &cuttree_R,        2},
  {"cor_pearson_rowSums_R", (DL_FUNC) &cor_pearson_rowSums_R, 2},
  {NULL, NULL, 0}
};


void R_init_island(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}





