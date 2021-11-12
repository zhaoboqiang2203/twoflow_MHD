#pragma once
#include <stdio.h>
#include <memory.h>
#include <fstream>
#include <iostream>

#include <string> 
#include <queue>
#include <cmath>

#define DADI_DEBUG
//#define FLUID_DEBUG
#define BOUNDARY_DEBUG

#define _for(i,a,b) for( i=(a); i<(b); ++i)
#define _feq(i,a,b) for( i=(a); i<=(b); ++i)

//	macros
#define MIN(a,b)			((a<b) ? (a) : (b))
#define MAX(a,b)			((a>b) ? (a) : (b))
#define sqr(a)	      ((a)*(a))
#define cube(a)       ((a)*(a)*(a))


#define RMAX 101
#define ZMAX 101

#define R_DIR 0
#define Z_DIR 1

#define LEFT 3
#define UP 5
#define RIGHT 7
#define DOWN 11

#define BND_NUM  8
//IB() = 1 vaccum， 2 dielectric surface， 3 extern interior  4 periodic  5 inlet 6 mirror reflection 7 anode 8 dielectric but fix phi
//9 轴对称边界条件（粒子弹性反射，电势梯度为0）10 无穷远边界条件（粒子消失，电势梯度为0） 11导体边界（粒子消失，电势给定或通过电容公式计算得到)

//int ib[num_rgn][4] = { {7,  3, 9, 2},
//					{3, 10, 9, 3},
//					{11, 5, 3,10} };	//边界条件

enum PTypes {
	VACCUM_BOUNDARY = 1,                 //真空边界
	DIELECTRIC_SURFACE_BOUNDARY,         //介质边界
	EXTERN_INTERIOR_BOUNDARY,			 //内部边界
	PERIODIC_BOUNDARY,					 //周期边界
	INLET,								 //内部边界
	MIRROR_REFLECTION_BOUNDARY,          //反射边界
	ANODE_BOUNDARY,						 //阳极边界
	CATHODE_BOUNDARY,					 //阴极边界
	CYLINDRICAL_AXIS,			         //对称轴
	CONDUCTING_BOUNDARY			         //导体边界
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
};




extern Boundary boundary_array[BND_NUM];
extern int nr, nz;
extern double dr, dz;
extern double dt;
extern int scale;

extern double world[RMAX][ZMAX];
extern double btype[RMAX][ZMAX];
extern double ptype[RMAX][ZMAX];
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
extern double phi[RMAX][ZMAX];
extern double rho[RMAX][ZMAX];
extern double phi1[RMAX][ZMAX];
extern double rou[RMAX][ZMAX];

extern double Er[RMAX][ZMAX];
extern double Ez[RMAX][ZMAX];

extern double ne[RMAX][ZMAX];
extern double ni[RMAX][ZMAX];

extern struct node MPDT[RMAX][ZMAX];
extern struct _U U[RMAX][ZMAX], U_bar[RMAX][ZMAX], U_bar2[RMAX][ZMAX];
extern struct _F Fr[RMAX][ZMAX], Fr_bar[RMAX][ZMAX], Fr_bar2[RMAX][ZMAX];
extern struct _F Fz[RMAX][ZMAX], Fz_bar[RMAX][ZMAX], Fz_bar2[RMAX][ZMAX];
extern struct _F s[RMAX][ZMAX], s_bar[RMAX][ZMAX], s_bar2[RMAX][ZMAX];
extern double app_Bz[RMAX][ZMAX];
extern double app_Br[RMAX][ZMAX];

void move();

int solveGS();
void test_sor_code();
void outputcsv();
void output_u_all();
void output_u(int n);

void magnetic_field_initial();