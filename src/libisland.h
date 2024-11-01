#ifndef LIBISLAND_H_
#define LIBISLAND_H_

#include <stddef.h> 

typedef struct {
			double k0, k1;
		} data_elem_t;

typedef struct {
	double nllr, d;
		} distance_t;

/* NB. size in doubles */
int add(int a, int b);

int subtract(int a, int b);

double bino_div(double x, double y);

double bino_div_id(int x, double y);

double bino_xlogy(double x, double y);

double bino_xlogy_id(int x, double y);

double cost(double k0, double k1);

double nLLR(const data_elem_t *K1, const data_elem_t *K2, double LOG_L1, double LOG_L2, size_t n);

distance_t merged_cost(const data_elem_t *K1, const data_elem_t *K2, double LOG_L1, double LOG_L2, const int chr1, const int chr2, const int pos1, const int pos2, size_t n);

void init_state(distance_t *D, double *L, size_t *PL, size_t *PR, const data_elem_t *K, double *LOG_L, double *tracked_cost, const int *chr, const int *pos, size_t m, size_t n);

size_t new_pair(size_t p, distance_t *D, size_t *PL, size_t *PR, size_t m);

double max_d(double A, double B);

void update_K(data_elem_t *K, double *LOG_L, size_t p1, size_t p2, size_t n);

double logl_for_clust(size_t *clust, size_t value, double *LOG_L, size_t m);

void dump_tree(const double *Z, size_t m1);


size_t sort_tree_size(size_t m);

size_t island_cluster_size(size_t m, size_t n);

void sort_tree(double *Z_new, void *scratch, const double *Z, size_t m);

void island_cluster(double *Z, const int *K0, const int *K1, const int *chr, const int *pos, size_t m, size_t n);

void cutree(size_t *L, const double *Z, size_t k, size_t n);

void cor_pearson_rowSums( double *sxy, double *sy1, double *sy2, const double *x, const double *Y, size_t m, size_t ny);

#endif /* LIBISLAND_H_ */