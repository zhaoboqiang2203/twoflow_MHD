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
int main()
{
	char a;
	nz = ZMAX;
	nr = RMAX;

	scale = ZMAX / 100;

	dr = 0.001 / scale;
	dz = 0.001 / scale;
	dt = 0.1 * ((dr * dr) + (dz * dz));
	printf("dt = %e\n", dt);
	initial();
	magnetic_field_initial();
	index = 300;
	while (index--)
	{

		printf("index %d\n", index);
		boundary_condition();
		//electron_flow();
		ion_flow();
		potential_solve();
		move();

		
		//if (index % 1000 == 0)
		{
			output();
		}
		
	}

	//scanf_s("%c", &a);
	return 0;
}

void matrix_to_csv(double** a, int N, int M, int array_size, char* filename)
{
	//fstream myfile(".\\output\\test_data.csv", ios::out);
	//fstream myfile(filename, ios::app);
	fstream myfile(filename, ios::out);
	if (!myfile.is_open())
	{
		cout << "δ�ɹ����ļ�" << endl;
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
	// ��������ܶ�
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
	//��������ٶ�

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


	//��������ܶ�
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ni;
		}
	}
	sprintf_s(fname, (".\\output\\ion density\\ion_density_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	// ��������ٶ�
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

	// ������Ʒֲ�

	sprintf_s(fname, (".\\output\\phi\\phi_%d.csv"), index);
	matrix_to_csv((double**)phi, ZMAX, RMAX, RMAX, fname);

	//sprintf_s(fname, (".\\output\\phi\\phi1_%d.csv"), index);
	//matrix_to_csv((double**)phi1, ZMAX, RMAX, RMAX, fname);

	sprintf_s(fname, (".\\output\\rho\\rho_%d.csv"), index);
	matrix_to_csv((double**)rho, ZMAX, RMAX, RMAX, fname);
	// ����ų��ֲ�

	//�������ѹǿ�ֲ�
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].pe;
		}
	}

	sprintf_s(fname, (".\\output\\electron pe\\electron_pe_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//�������ѹǿ�ֲ�
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].pi;
		}
	}

	sprintf_s(fname, (".\\output\\ion pi\\ion_pi_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	// ������������ֲ�
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ee;
		}
	}

	sprintf_s(fname, (".\\output\\electron ee\\electron_ee_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	// ������������ֲ�
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ei;
		}
	}

	sprintf_s(fname, (".\\output\\ion ei\\ion_ei_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


	//Ez�ֲ�
	sprintf_s(fname, (".\\output\\Ez\\Ez_%d.csv"), index);
	matrix_to_csv((double**)Ez, ZMAX, RMAX, RMAX, fname);

	//Er�ֲ�
	sprintf_s(fname, (".\\output\\Er\\Er_%d.csv"), index);
	matrix_to_csv((double**)Er, ZMAX, RMAX, RMAX, fname);

	// ��������ų��ֲ�
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].br;
		}
	}

	sprintf_s(fname, (".\\output\\br\\br_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].br;
		}
	}

	sprintf_s(fname, (".\\output\\btheta\\btheta_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].br;
		}
	}

	sprintf_s(fname, (".\\output\\bz\\bz_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//����������ײƵ��
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].mu_ie;
		}
	}

	sprintf_s(fname, (".\\output\\mu_ie\\mu_ie_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	
	//�����ܶ�r����
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = (MPDT[i][j].ni * MPDT[i][j].vir - MPDT[i][j].ne * MPDT[i][j].ver) * QE / dt;
		}
	}

	sprintf_s(fname, (".\\output\\Jr\\Jr_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//�����ܶ�z����
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = (MPDT[i][j].ni * MPDT[i][j].viz - MPDT[i][j].ne * MPDT[i][j].vez) * QE / dt;
		}
	}

	sprintf_s(fname, (".\\output\\Jz\\Jz_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

}


void output_u(int n)
{

	if (n == 1)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u1;
			}
		}
	}
	else if (n == 2)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u2;
			}
		}
	}
	else if (n == 3)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u3;
			}
		}
	}
	else if (n == 4)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u4;
			}
		}
	}
	else if (n == 5)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u5;
			}
		}
	}
	else if (n == 6)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u6;
			}
		}
	}
	else if (n == 7)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u7;
			}
		}
	}
	else if (n == 8)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u8;
			}
		}
	}
	else if (n == 9)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u9;
			}
		}
	}
	else if (n == 10)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u10;
			}
		}
	}
	else if (n == 11)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u11;
			}
		}
	}
	else if (n == 12)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u12;
			}
		}
	}
	else if (n == 13)
	{
		for (int i = 0; i < ZMAX; i++)
		{
			for (int j = 0; j < RMAX; j++)
			{
				res_out[i][j] = U[i][j].u13;
			}
		}
	}


	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\test U.csv"));
}

void output_u_all()
{
	


	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u1;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U1.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u2;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U2.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u3;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U3.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u4;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U4.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u5;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U5.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u6;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U6.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u7;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U7.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u8;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U8.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u9;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U9.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u10;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U10.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u11;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U11.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u12;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U12.csv"));

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = U[i][j].u13;
		}
	}
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\U13.csv"));
}