#pragma once
/* The actual coefficients of the differencing on the mesh,
	and the influence of the spatial variation of epsilon:  */
#include "MHD.h"

void init_solve();
void tridag(float* a, float* b, float* c, float* r, float* utri, float* gam, int n);
int solve(float u_in[RMAX][ZMAX], float s[RMAX][ZMAX], int itermax, float tol_test);