#pragma once
#include <stdio.h>
#include <memory.h>
#include <fstream>
#include <iostream>

#include <string> 
#include <queue>
#include <cmath>

//#define DADI_DEBUG
//#define FLUID_DEBUG
#define BOUNDARY_DEBUG

#define _for(i,a,b) for( i=(a); i<(b); ++i)
#define _feq(i,a,b) for( i=(a); i<=(b); ++i)

//	macros
#define MIN(a,b)			((a<b) ? (a) : (b))
#define MAX(a,b)			((a>b) ? (a) : (b))
#define sqr(a)	      ((a)*(a))
#define cube(a)       ((a)*(a)*(a))


#define RMAX 401
#define ZMAX 1001

#define R_DIR 0
#define Z_DIR 1

#define LEFT 2
#define UP 8
#define RIGHT 32
#define DOWN 128

#define BND_NUM  10


enum PTypes {
	VACCUM_BOUNDARY = 1,                 //真空边界
	DIELECTRIC_SURFACE_BOUNDARY = 2,         //介质边界
	EXTERN_INTERIOR_BOUNDARY = 4,			 //内部边界
	PERIODIC_BOUNDARY = 8,					 //周期边界
	INLET = 16,								 //内部边界
	MIRROR_REFLECTION_BOUNDARY = 32,          //反射边界
	ANODE_BOUNDARY = 64,						 //阳极边界
	CATHODE_BOUNDARY = 128,					 //阴极边界
	CYLINDRICAL_AXIS = 256,			         //对称轴
	CONDUCTING_BOUNDARY = 512			         //导体边界
};



struct Point
{
	double r;
	double z;
};

struct Boundary
{
	PTypes physics_type;

	Point start; //边界起始位置
	Point end;   //边界终止位置

	int bnd_dir;  // 边界的方向 0 r 方向，1 z 方向
	int boundary_type;
};

struct _U
{
	int r, z;
	double u1;
	double u2;
	double u3;
	double u4;
	double u5;
	double u6;
	double u7;
	double u8;
	double u9;
	double u10;
	double u11;
	double u12;
	double u13;
};

struct _F
{
	int r, z;
	double f1;
	double f2;
	double f3;
	double f4;
	double f5;
	double f6;
	double f7;
	double f8;
	double f9;
	double f10;
	double f11;
	double f12;
	double f13;
};

/* node interface with plasma information */

struct  node
{
	double ne;
	double ni;
	double ver;
	double vez;
	double vetheta;
	double vir;
	double viz;
	double vitheta;
	double br;
	double bz;
	double btheta;
	double pe;
	double pi;
	double ee;
	double ei;
	double mu_ie;       /*电子离子碰撞频率*/
	double delta_ei;    /*电子向离子转移能量*/
};




extern Boundary boundary_array[BND_NUM];
extern int nr, nz;
extern double dr, dz;
extern double dt;
extern int scale;
extern int index;

extern double world[ZMAX][RMAX];
extern double btype[ZMAX][RMAX];
extern double ptype[ZMAX][RMAX];
int initial();
void  boundary_condition();
void matrix_to_csv(double** a, int N, int M, int array_size, char* filename);
int fill_plasma(int tr, int tz, int fill_n);
void output();

const double EPS_0 = 8.85418782e-12;  	// C/(V*m), vacuum permittivity
const double QE = 1.602176565e-19;		// C, electron charge
const double AMU = 1.660538921e-27;		// kg, atomic mass unit
const double ME = 9.10938215e-31;		// kg, electron mass
const double MI = 40 * AMU;		// kg, electron mass
const double K = 1.380648e-23;			// J/K, Boltzmann constant
const double PI = 3.141592653;			// pi
const double EvToK = QE / K;				// 1eV in K ~ 11604
const double gamma = 1.4;
const double MU_0 = 4 * PI * 1e-7;
extern double phi[ZMAX][RMAX];
extern double rho[ZMAX][RMAX];
extern double phi1[ZMAX][RMAX];
extern double rou[ZMAX][RMAX];

extern double Er[ZMAX][RMAX];
extern double Ez[ZMAX][RMAX];

extern double btheta[ZMAX][RMAX];

extern double ne[ZMAX][RMAX];
extern double ni[ZMAX][RMAX];

extern struct node MPDT[ZMAX][RMAX];
extern struct _U U[ZMAX][RMAX], U_bar[ZMAX][RMAX], U_bar2[ZMAX][RMAX];
extern struct _F Fr[ZMAX][RMAX], Fr_bar[ZMAX][RMAX], Fr_bar2[ZMAX][RMAX];
extern struct _F Fz[ZMAX][RMAX], Fz_bar[ZMAX][RMAX], Fz_bar2[ZMAX][RMAX];
extern struct _F s[ZMAX][RMAX], s_bar[ZMAX][RMAX], s_bar2[ZMAX][RMAX];
extern double app_Bz[ZMAX][RMAX];
extern double app_Br[ZMAX][RMAX];

void move();

int solveGS();
void test_sor_code();
void outputcsv();
void output_u_all();
void output_u(int n);

void magnetic_field_initial();