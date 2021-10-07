#pragma once
/* The actual coefficients of the differencing on the mesh,
	and the influence of the spatial variation of epsilon:  */
#include "MHD.h"

float a_x1geom[RMAX][ZMAX],  b_x1geom[RMAX][ZMAX],  c_x1geom[RMAX][ZMAX];
float a_x2geom[RMAX][ZMAX],  b_x2geom[RMAX][ZMAX],  c_x2geom[RMAX][ZMAX];

/*  The arrays used internal to DADI which contain the coefficients
	  for each tridiagonal matrix solution */
float* a_x1, * b_x1, * c_x1;
float* a_x2, * b_x2, * c_x2;
/*  Various copies of the 'answer' we're working toward */
float u[RMAX][ZMAX];
float uwork[RMAX][ZMAX], ustor[RMAX][ZMAX], ustar[RMAX][ZMAX];
float* r_x1, * v_x1, * gam_x1;
float* r_x2, * v_x2, * gam_x2;
/*  Our fictitious time step */
float del_t0;
/*  epsilon at the grid locations */
float** epsi;
/*  The size of the system */


void tridag(float* a, float* b, float* c, float* r, float* utri, float* gam, int n);
int solve(float u_in[RMAX][ZMAX], float s[RMAX][ZMAX], int itermax, float tol_test);