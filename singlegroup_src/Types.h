#ifndef __TYPES_HEADER_H__
#define __TYPES_HEADER_H__

// Structures
typedef struct  __attribute__((packed)) __attribute__((aligned(64))) {
	double energy;
	double total_xs;
	double elastic_xs;
	double absorbtion_xs;
	double fission_xs;
	double nu_fission_xs;
//	double dummy1;
//	double dummy2;
} NuclideGridPoint;

typedef struct{
	double energy;
	int * xs_ptrs;
} GridPoint;

typedef struct __attribute__((aligned(16))) {
	double data;
	long index;
} CacheData;

typedef struct{
	int nthreads;
	long n_isotopes;
	long n_gridpoints;
	int lookups;
	char * HM;
} Inputs;

#endif

