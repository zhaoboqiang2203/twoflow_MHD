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
//IB() = 1 vaccum�� 2 dielectric surface�� 3 extern interior  4 periodic  5 inlet 6 mirror reflection 7 anode 8 dielectric but fix phi
//9 ��ԳƱ߽����������ӵ��Է��䣬�����ݶ�Ϊ0��10 ����Զ�߽�������������ʧ�������ݶ�Ϊ0�� 11����߽磨������ʧ�����Ƹ�����ͨ�����ݹ�ʽ����õ�)

//int ib[num_rgn][4] = { {7,  3, 9, 2},
//					{3, 10, 9, 3},
//					{11, 5, 3,10} };	//�߽�����

enum PTypes {
	VACCUM_BOUNDARY = 1,                 //��ձ߽�
	DIELECTRIC_SURFACE_BOUNDARY,         //���ʱ߽�
	EXTERN_INTERIOR_BOUNDARY,			 //�ڲ��߽�
	PERIODIC_BOUNDARY,					 //���ڱ߽�
	INLET,								 //�ڲ��߽�
	MIRROR_REFLECTION_BOUNDARY,          //����߽�
	ANODE_BOUNDARY,						 //�����߽�
	CATHODE_BOUNDARY,					 //�����߽�
	CYLINDRICAL_AXIS,			         //�Գ���
	CONDUCTING_BOUNDARY			         //����߽�
};



struct Point
{
	float r;
	float z;
};

struct Boundary
{
	PTypes physics_type;

	Point start; //�߽���ʼλ��
	Point end;   //�߽���ֹλ��

	int bnd_dir;  // �߽�ķ��� 0 r ����1 z ����
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