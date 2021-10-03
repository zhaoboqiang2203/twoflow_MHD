#include "MHD.h"
#include "fluid.h"
using namespace std;

int nr;
int nz;
float phi[RMAX][ZMAX];

struct node MPDT[RMAX][ZMAX];

int main()
{
	nr = RMAX;
	nz = ZMAX;
	initial();
	
	for(int it = 0;it < 1000; it++)
	{
		electron_flow();
		ion_flow();

		//todo �������

		output();
	}

	return 0;
}


void matrix_to_csv(float** a, int N, int M, int array_size, char* filename)
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
	//todo ��������ܶ�
	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ne;
		}
	}

	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\electron density.csv"));
	//��������ٶ�

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ver;
		}
	}

	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\electron ver.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vez;
		}
	}

	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\electron vez.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vetheta;
		}
	}

	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\electron vetheta.csv"));

	//��������ܶ�
	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ni;
		}
	}

	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\ion density.csv"));
	//todo ��������ٶ�
	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vir;
		}
	}

	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\electron vir.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].viz;
		}
	}

	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\electron viz.csv"));

	for (int i = 0; i < RMAX; i++)
	{
		for (int j = 0; j < ZMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vitheta;
		}
	}

	matrix_to_csv((float**)res_out, ZMAX, RMAX, RMAX, (char*)(".\\output\\electron vitheta.csv"));

	//todo ������Ʒֲ�

	matrix_to_csv((float**)phi, ZMAX, RMAX, RMAX, (char*)(".\\output\\electron vetheta.csv"));
	//todo ����ų��ֲ�

	//todo ��������¶ȷֲ�

	//todo ��������¶ȷֲ�

}