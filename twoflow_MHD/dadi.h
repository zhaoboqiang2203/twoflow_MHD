#pragma once
/* The actual coefficients of the differencing on the mesh,
	and the influence of the spatial variation of epsilon:  */
#include "MHD.h"

void init_solve();
void tridag(double* a, double* b, double* c, double* r, double* utri, double* gam, int n);
int solve(double u_in[RMAX][ZMAX], double s[RMAX][ZMAX], int itermax, double tol_test);

void electric_field();
void potential_solve();