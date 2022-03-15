#pragma once
#include <stdio.h>
#include <memory.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string> 
#include <queue>
#include <cmath>

#include <fstream>
#include <string>
#include <sstream>

using namespace std;

//#define DADI_DEBUG
//#define FLUID_DEBUG
//#define BOUNDARY_DEBUG

#define _for(i,a,b) for( i=(a); i<(b); ++i)
#define _feq(i,a,b) for( i=(a); i<=(b); ++i)
#define sgn(a) (a >= 0 ? 1 : -1)

//	macros
#define MIN(a,b)		((a < b) ? (a) : (b))
#define MAX(a,b)		((a > b) ? (a) : (b))
#define sqr(a)	        ((a) * (a))
#define cube(a)         ((a) * (a) * (a))

#define PRE_DIFF_U_R(i,j,x) (-3*U[i][j].u[x]+3*U[i][j+1].u[x]-U[i][j+2].u[x])
#define PRE_DIFF_U_Z(i,j,x) (-3*U[i][j].u[x]+3*U[i+1][j].u[x]-U[i+2][j].u[x])

#define MID_DIFF_U_R(i,j,x) (-U[i][j-1].u[x]+U[i][j+1].u[x])
#define MID_DIFF_U_Z(i,j,x) (-U[i-1][j].u[x]+U[i+1][j].u[x])

#define AFTER_DIFF_U_R(i,j,x) (U[i][j-2].u[x]-4*U[i][j-1].u[x]+3*U[i][j].u[x])
#define AFTER_DIFF_U_Z(i,j,x) (U[i-2][j].u[x]-4*U[i-1][j].u[x]+3*U[i][j].u[x])

#define RMAX 401
#define ZMAX 1001

#define R_DIR 0
#define Z_DIR 1

#define LEFT 2
#define UP 8
#define RIGHT 32
#define DOWN 128

#define BND_NUM  30


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
	double u[13];
};

struct _F
{
	int r, z;
	double f[13];

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
	double sigma_Q;     /*电子离子碰撞截面*/
	double tau_ei;      /*电子离子弛豫时间*/

	double neq;       /*发射电子密度*/
	double vnqr;
	double vnqz;
	double vnqtheta;

	double peq;       /*发射的虚拟电子密度*/
	double vpqr;
	double vpqz;
	double vpqtheta;

	double angle_b_vi;
	
};

struct _particle
{
	double den;
	double vr;
	double vtheta;
	double vz;

	double pressure;
	double eng;
};

struct _AF
{
	int r, z;
	double f[5];

};

struct _pid 
{
	float set_current;//定义设定值
	float actual_current;//定义实际值
	float err;//定义偏差值
	float err_last;//定义上一个偏差值
	float Kp, Ki, Kd;//定义比例、积分、微分系数
	float ne_density;//定义发射电子密度(控制执行器的变量)
	float integral;//定义积分值
};



extern Boundary boundary_array[BND_NUM];
extern int nr, nz;
extern double dr, dz, dtheta;
extern double dt;
extern int scale;
extern int index;

extern int world[ZMAX][RMAX];
extern int btype[ZMAX][RMAX];
extern int ptype[ZMAX][RMAX];

extern double cathod_cell;
extern double anode_cell;
extern double out_e_den;
extern int bnd_size;

int initial();
void  boundary_condition();
void matrix_to_csv(double** a, int N, int M, int array_size, char* filename);
void matrix_int_to_csv(int** a, int N, int M, int array_size, char* filename);
int fill_plasma(int tr, int tz, int fill_n);
void output();


const double EPS_0 = 8.85418782e-12;  	// C/(V*m), vacuum permittivity
//const double EPS_PLA = 1e8;				//相对电导率
const double QE = 1.602176565e-19;		// C, electron charge
const double AMU = 1.660538921e-27;		// kg, atomic mass unit
const double ME = 9.10938215e-31;		// kg, electron mass
//const double REL_MASS = 40;             //相对原子质量
//const double MI = REL_MASS * AMU;		        // kg, electron mass
const double K = 1.380648e-23;			// J/K, Boltzmann constant
const double PI = 3.141592653;			// pi
const double EvToK = QE / K;				// 1eV in K ~ 11604
const double gamma = 1.4;
const double MU_0 = 4 * PI * 1e-7;

const double NA = 6.02e23;             //Avogadro constant
extern double phi[ZMAX][RMAX];
extern double rho[ZMAX][RMAX];
extern double phi1[ZMAX][RMAX];
extern double rou[ZMAX][RMAX];

extern double Er[ZMAX][RMAX];
extern double Ez[ZMAX][RMAX];

extern double btheta[ZMAX][RMAX];

extern double Jz[ZMAX][RMAX];
extern double Jr[ZMAX][RMAX];

extern struct _particle atom[ZMAX][RMAX];
extern struct node MPDT[ZMAX][RMAX];
extern struct _U U[ZMAX][RMAX], U_bar[ZMAX][RMAX], U_bar2[ZMAX][RMAX];
extern struct _F Fr[ZMAX][RMAX], Fr_bar[ZMAX][RMAX], Fr_bar2[ZMAX][RMAX];
extern struct _F Fz[ZMAX][RMAX], Fz_bar[ZMAX][RMAX], Fz_bar2[ZMAX][RMAX];
extern struct _F s[ZMAX][RMAX], s_bar[ZMAX][RMAX], s_bar2[ZMAX][RMAX];


extern struct _AF Uq[ZMAX][RMAX], Uq_bar[ZMAX][RMAX], Uq_bar2[ZMAX][RMAX];
extern struct _AF Fqr[ZMAX][RMAX], Fqr_bar[ZMAX][RMAX], Fqr_bar2[ZMAX][RMAX];
extern struct _AF Fqz[ZMAX][RMAX], Fqz_bar[ZMAX][RMAX], Fqz_bar2[ZMAX][RMAX];
extern struct _AF sq[ZMAX][RMAX], sq_bar[ZMAX][RMAX], sq_bar2[ZMAX][RMAX];

extern struct _F tau_vis[ZMAX][RMAX];

extern double taurr[ZMAX][RMAX];
extern double taurtheta[ZMAX][RMAX];
extern double taurz[ZMAX][RMAX];
extern double tautheta2[ZMAX][RMAX];
extern double tauthetaz[ZMAX][RMAX];
extern double tauzz[ZMAX][RMAX];

extern double pniz[ZMAX][RMAX];
extern double pnir[ZMAX][RMAX];
extern double tau_ni[ZMAX][RMAX];

extern double vez[ZMAX][RMAX];
extern double ver[ZMAX][RMAX];
extern double vethera[ZMAX][RMAX];

extern double viz[ZMAX][RMAX];
extern double vir[ZMAX][RMAX];
extern double vithera[ZMAX][RMAX];

extern double app_Bz[ZMAX][RMAX];
extern double app_Br[ZMAX][RMAX];

extern struct _pid pid;

extern double bg_den;
extern double inter_e_den;
extern double inter_pla_den;
extern double max_phi;
extern double set_phi;

extern double e_half[ZMAX][RMAX];
extern double e_hrho[ZMAX][RMAX], e_htheta[ZMAX][RMAX], e_hz[ZMAX][RMAX], e_h2[ZMAX][RMAX];
extern double e_srho[ZMAX][RMAX], e_stheta[ZMAX][RMAX], e_sz[ZMAX][RMAX];

extern double i_half[ZMAX][RMAX];
extern double i_hrho[ZMAX][RMAX], i_htheta[ZMAX][RMAX], i_hz[ZMAX][RMAX], i_h2[ZMAX][RMAX];
extern double i_srho[ZMAX][RMAX], i_stheta[ZMAX][RMAX], i_sz[ZMAX][RMAX];

extern double cathode_I;
extern double current_I;
extern double last_I;
extern double coil_I;//线圈电流
extern double coil_R;//线圈内径
extern double coil_W;//线圈宽度（内外径差值）
extern double coil_L;//线圈厚度

extern double REL_MASS;             //相对原子质量
extern double MI;		            // kg, electron mass
extern double EPS_PLA;				//相对电导率

extern double para_p, para_i, para_d;

/// <summary>
/// move.cpp 函数声明
/// </summary>
void move();
double  magnetic_vec_angle(double var, double vaz, double vbr, double vbz);
bool is_electron_ion_separation(double angle);
bool is_large_max_speed(double ur, double utheta, double uz, double max_speed);
void ionization_collisions(int i, int j);
void coulomb_collision(int i, int j);
void current_caulate();
void current_control();
void viscosity_collision(int i, int j);

int solveGS();
void test_sor_code();
void outputcsv();
void read_csv();
void out_judge();
void wirte_datfile();
void read_datfile();
bool is_read_datfile();
int parameter_read();
void max_write();

int btype_trans(string ss);
int dir_trans(string ss);
PTypes ptype_trans(string ss);

void judge_bit();
void magnetic_field_initial();
void magnetic_display();

int judge_conner(int i, int j);