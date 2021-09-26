#include "MHD.h"
#include "fluid.h"
using namespace std;

int nr;
int nz;

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
	//todo ��������ܶ�

	//todo ��������ٶ�

	//todo ��������ܶ�

	//todo ��������ٶ�

	//todo ������Ʒֲ�

	//todo ����ų��ֲ�

	//todo ��������¶ȷֲ�

	//todo ��������¶ȷֲ�

}