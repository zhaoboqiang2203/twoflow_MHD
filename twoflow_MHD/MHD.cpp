#include "MHD.h"
#include "fluid.h"
#include "dadi.h"
#include "testcode.h"
using namespace std;

int nr;
int nz;
float dr, dz;
float dt;

int index;

struct node MPDT[RMAX][ZMAX];
struct _U U[RMAX][ZMAX], U_bar[RMAX][ZMAX], U_bar2[RMAX][ZMAX];
struct _F Fr[RMAX][ZMAX], Fr_bar[RMAX][ZMAX], Fr_bar2[RMAX][ZMAX];
struct _F Fz[RMAX][ZMAX], Fz_bar[RMAX][ZMAX], Fz_bar2[RMAX][ZMAX];
struct _F s[RMAX][ZMAX], s_bar[RMAX][ZMAX], s_bar2[RMAX][ZMAX];
float phi[RMAX][ZMAX];
float rho[RMAX][ZMAX];

int scale;

float vez[RMAX][ZMAX];
float ver[RMAX][ZMAX];
float vethera[RMAX][ZMAX];

float viz[RMAX][ZMAX];
float vir[RMAX][ZMAX];
float vithera[RMAX][ZMAX];

float bapp;
float denJ[RMAX][ZMAX];

int main()
{
	char a;
	nz = ZMAX;
	nr = RMAX;

	scale = RMAX / 100;

	dr = 0.001, dz = 0.001;
	dt = 0.1 * ((dr * dr) + (dz * dz));
	printf("dt = %lf\n", dt);
	initial();
	index = 20;
	while (index--)
	{

		printf("index %d\n", index);
		electron_flow();

		potential_solve();
		move();

		boundary_condition();
		output();
	}

	scanf_s("%c", &a);
	return 0;
}

void matrix_to_csv(float** a, int N, int M, int array_size, char* filename)
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
			//cout << *((float*)a + (array_size * i + j)) << endl;
			myfile << *((float*)a + (array_size * i + j)) << ",";
		}
		myfile << endl;
	}
	myfile.close();
}

void output()
{
	float res_out[RMAX][ZMAX];
	//todo 输出电子密度
	char fname[100];
	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ne;
		}
	}
	sprintf_s(fname, (".\\output\\electron density\\electron_density_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);
	//输出电子速度

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ver;
		}
	}
	sprintf_s(fname, (".\\output\\electron ver\\electron_ver_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vez;
		}
	}
	sprintf_s(fname, (".\\output\\electron vez\\electron_vez_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vetheta;
		}
	}

	sprintf_s(fname, (".\\output\\electron vetheta\\electron_vetheta_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);


	//输出离子密度
	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ni;
		}
	}
	sprintf_s(fname, (".\\output\\ion density\\ion_density_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);

	//todo 输出离子速度
	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vir;
		}
	}

	sprintf_s(fname, (".\\output\\ion vir\\ion_vir_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);


	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].viz;
		}
	}

	sprintf_s(fname, (".\\output\\ion viz\\ion_viz_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vitheta;
		}
	}

	sprintf_s(fname, (".\\output\\ion vitheta\\ion_vitheta_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);

	//todo 输出电势分布

	sprintf_s(fname, (".\\output\\phi\\phi_%d.csv"), index);
	matrix_to_csv((float**)phi, ZMAX, RMAX, RMAX, fname);


	sprintf_s(fname, (".\\output\\rho\\rho_%d.csv"), index);
	matrix_to_csv((float**)rho, ZMAX, RMAX, RMAX, fname);
	//todo 输出磁场分布

	//输出电子压强分布
	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].pe;
		}
	}

	sprintf_s(fname, (".\\output\\electron pe\\electron_pe_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);

	//输出离子压强分布
	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].pi;
		}
	}

	sprintf_s(fname, (".\\output\\ion pi\\ion_pi_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);
	//todo 输出电子能量分布
	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ee;
		}
	}

	sprintf_s(fname, (".\\output\\electron ee\\electron_ee_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);
	//todo 输出离子能量分布
	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ei;
		}
	}

	sprintf_s(fname, (".\\output\\ion ei\\ion_ei_%d.csv"), index);
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, fname);

}


void output_u(int n)
{
	float res_out[RMAX][ZMAX];

	if (n == 1)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u1;
			}
		}
	}
	else if (n == 2)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u2;
			}
		}
	}
	else if (n == 3)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u3;
			}
		}
	}
	else if (n == 4)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u4;
			}
		}
	}
	else if (n == 5)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u5;
			}
		}
	}
	else if (n == 6)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u6;
			}
		}
	}
	else if (n == 7)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u7;
			}
		}
	}
	else if (n == 8)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u8;
			}
		}
	}
	else if (n == 9)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u9;
			}
		}
	}
	else if (n == 10)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u10;
			}
		}
	}
	else if (n == 11)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u11;
			}
		}
	}
	else if (n == 12)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u12;
			}
		}
	}
	else if (n == 13)
	{
		for (int i = 0; i < RMAX; i++)
		{
			for (int j = 0; j < ZMAX; j++)
			{
				res_out[i][j] = U[i][j].u13;
			}
		}
	}


	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\test U.csv"));
}

void output_u_all()
{
	float res_out[RMAX][ZMAX];


	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u1;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U1.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u2;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U2.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u3;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U3.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u4;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U4.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u5;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U5.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u6;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U6.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u7;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U7.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u8;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U8.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u9;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U9.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u10;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U10.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u11;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U11.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u12;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U12.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = U[i][j].u13;
		}
	}
	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U13.csv"));
}