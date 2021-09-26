#pragma once
#include <stdio.h>
#include <memory.h>
#include <fstream>
#include <iostream>

#include <string> 
#include <queue>


#define _for(i,a,b) for( i=(a); i<(b); ++i)
#define _feq(i,a,b) for( i=(a); i<=(b); ++i)


#define RMAX 101
#define ZMAX 101

#define R_DIR 0
#define Z_DIR 1

#define LEFT 2
#define RIGHT 3
#define UP 4
#define DOWN 5


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
	float r;
	float z;
};

struct Boundary
{
	PTypes physics_type;

	Point start; //边界起始位置
	Point end;   //边界终止位置

	int bnd_dir;  // 边界的方向 0 r 方向，1 z 方向
	int boundary_type;
};

struct U
{
	float roe;
	float vz;
	float vr;
	float vtheta;

};

struct  node
{
	float ne;
	float ni;
	float ver;
	float vez;
	float vetheta;
	float vir;
	float viz;
	float vitheta;
	float br;
	float bz;
	float btheta;
	float pe;
	float pi;
	float ee;
	float ei;
};




extern Boundary boundary_array[BND_NUM];
extern int nr, nz;


int initial();
void matrix_to_csv(float** a, int N, int M, int array_size, char* filename);
int fill_plasma(int tr, int tz, int fill_n);

void output();