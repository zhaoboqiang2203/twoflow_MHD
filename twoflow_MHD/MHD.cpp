#include "MHD.h"
#include "fluid.h"
#include "dadi.h"
#include "testcode.h"
using namespace std;

int nr;
int nz;
double dr, dz;
double dt;

int index;

struct node MPDT[ZMAX][RMAX];
struct _U U[ZMAX][RMAX], U_bar[ZMAX][RMAX], U_bar2[ZMAX][RMAX];
struct _F Fr[ZMAX][RMAX], Fr_bar[ZMAX][RMAX], Fr_bar2[ZMAX][RMAX];
struct _F Fz[ZMAX][RMAX], Fz_bar[ZMAX][RMAX], Fz_bar2[ZMAX][RMAX];
struct _F s[ZMAX][RMAX], s_bar[ZMAX][RMAX], s_bar2[ZMAX][RMAX];

struct _U Uq[ZMAX][RMAX], Uq_bar[ZMAX][RMAX], Uq_bar2[ZMAX][RMAX];
struct _F Fqr[ZMAX][RMAX], Fqr_bar[ZMAX][RMAX], Fqr_bar2[ZMAX][RMAX];
struct _F Fqz[ZMAX][RMAX], Fqz_bar[ZMAX][RMAX], Fqz_bar2[ZMAX][RMAX];
struct _F sq[ZMAX][RMAX], sq_bar[ZMAX][RMAX], sq_bar2[ZMAX][RMAX];

double phi[ZMAX][RMAX];
double rho[ZMAX][RMAX];

double phi1[ZMAX][RMAX];
double rou[ZMAX][RMAX];

int scale;

double vez[ZMAX][RMAX];
double ver[ZMAX][RMAX];
double vethera[ZMAX][RMAX];

double viz[ZMAX][RMAX];
double vir[ZMAX][RMAX];
double vithera[ZMAX][RMAX];

double app_Bz[ZMAX][RMAX];
double app_Br[ZMAX][RMAX];
double denJ[ZMAX][RMAX];
double res_out[ZMAX][RMAX];

double btheta[ZMAX][RMAX];

double max_phi;
double set_phi;
int main()
{
	nz = ZMAX;
	nr = RMAX;
	index = 40001;
	set_phi = 160;
	scale = ZMAX / 200;

	dr = 0.001 / scale;
	dz = 0.001 / scale;
	dt = dr / 1e7;

	//根据背景压强，通气流量，电流密度计算
	bg_den = 1e-3 / (K * 300);
	inter_e_den = 300 * dt / QE / 360 / 50 * 1e9;
	inter_pla_den = 0.04 * dt / 40 * NA / 360 / 20 * 1e9;


	//dt = 0.05 * ((dr * dr) + (dz * dz));
	printf("dt = %e\n", dt);
	printf("inter_e_den = %e\n", inter_e_den);
	printf("inter_pla_den = %e\n", inter_pla_den);

	initial();
	magnetic_field_initial();
	//dadi initial
	init_solve();
	
	while (index--)
	{
		printf("index %d\n", index);
		boundary_condition();
		electron_flow();
		//ion_flow();

		if (index % 100 == 0)
		{
			Q_fluid();
			potential_solve();
			move_q();
		}
		
		move();
		mag_phi();
		if (index % 100 == 0)
		//if(index < 36000)
		{
			output();
		}
		
	}


	return 0;
}

void matrix_to_csv(double** a, int N, int M, int array_size, char* filename)
{
	//fstream myfile(".\\output\\test_data.csv", ios::out);
	//fstream myfile(filename, ios::app);
	fstream myfile(filename, ios::out);
	if (!myfile.is_open())
	{
		cout << "未成功打开文件" << endl;
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			//cout << *((double*)a + (array_size * i + j)) << endl;
			myfile << *((double*)a + (array_size * i + j)) << ",";
		}
		myfile << endl;
	}
	myfile.close();
}

void output()
{
	// 输出电子密度
	char fname[100];
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ne;
		}
	}
	sprintf_s(fname, (".\\output\\electron density\\electron_density_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	//输出电子速度

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ver;
		}
	}
	sprintf_s(fname, (".\\output\\electron ver\\electron_ver_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vez;
		}
	}
	sprintf_s(fname, (".\\output\\electron vez\\electron_vez_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vetheta;
		}
	}

	sprintf_s(fname, (".\\output\\electron vetheta\\electron_vetheta_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


	//输出离子密度
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ni;
		}
	}
	sprintf_s(fname, (".\\output\\ion density\\ion_density_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	// 输出离子速度
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vir;
		}
	}

	sprintf_s(fname, (".\\output\\ion vir\\ion_vir_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].viz;
		}
	}

	sprintf_s(fname, (".\\output\\ion viz\\ion_viz_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vitheta;
		}
	}

	sprintf_s(fname, (".\\output\\ion vitheta\\ion_vitheta_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	// 输出电势分布

	sprintf_s(fname, (".\\output\\phi\\phi_%d.csv"), index);
	matrix_to_csv((double**)phi, ZMAX, RMAX, RMAX, fname);

	//sprintf_s(fname, (".\\output\\phi\\phi1_%d.csv"), index);
	//matrix_to_csv((double**)phi1, ZMAX, RMAX, RMAX, fname);

	sprintf_s(fname, (".\\output\\rho\\rho_%d.csv"), index);
	matrix_to_csv((double**)rou, ZMAX, RMAX, RMAX, fname);
	// 输出磁场分布

	//输出电子压强分布
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].pe;
		}
	}

	sprintf_s(fname, (".\\output\\electron pe\\electron_pe_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//输出离子压强分布
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].pi;
		}
	}

	sprintf_s(fname, (".\\output\\ion pi\\ion_pi_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	// 输出电子能量分布
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ee / QE;
		}
	}

	sprintf_s(fname, (".\\output\\electron ee\\electron_ee_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	// 输出离子能量分布
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ei / QE;
		}
	}

	sprintf_s(fname, (".\\output\\ion ei\\ion_ei_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


	//Ez分布
	sprintf_s(fname, (".\\output\\Ez\\Ez_%d.csv"), index);
	matrix_to_csv((double**)Ez, ZMAX, RMAX, RMAX, fname);

	//Er分布
	sprintf_s(fname, (".\\output\\Er\\Er_%d.csv"), index);
	matrix_to_csv((double**)Er, ZMAX, RMAX, RMAX, fname);

	// 输出感生磁场分布
	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].br;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\br\\br_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].btheta;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\btheta\\btheta_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].bz;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\bz\\bz_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	 
	sprintf_s(fname, (".\\output\\btheta\\btheta_%d.csv"), index);
	matrix_to_csv((double**)btheta, ZMAX, RMAX, RMAX, fname);

	//电子离子碰撞频率
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].mu_ie;
		}
	}

	sprintf_s(fname, (".\\output\\mu_ie\\mu_ie_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	
	//电流密度r方向
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = (MPDT[i][j].ni * MPDT[i][j].vir - MPDT[i][j].ne * MPDT[i][j].ver) * QE / dt;
		}
	}

	sprintf_s(fname, (".\\output\\Jr\\Jr_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//电流密度z方向
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = (MPDT[i][j].ni * MPDT[i][j].viz - MPDT[i][j].ne * MPDT[i][j].vez) * QE / dt;
		}
	}

	sprintf_s(fname, (".\\output\\Jz\\Jz_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//电子向离子转移能量
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].delta_ei;
		}
	}

	sprintf_s(fname, (".\\output\\delta_ei\\delta_ei_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//电子离子碰撞截面
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].sigma_Q;
		}
	}

	sprintf_s(fname, (".\\output\\sigma\\sigma_Q_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//多余电荷密度
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].peq;
		}
	}

	sprintf_s(fname, (".\\output\\peq\\peq_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].neq;
		}
	}

	sprintf_s(fname, (".\\output\\neq\\neq_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
}
