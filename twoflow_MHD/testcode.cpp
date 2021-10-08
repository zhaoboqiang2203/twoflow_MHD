#include "testcode.h"

void test_potential()
{
	float res[RMAX][ZMAX];
	float test_s[RMAX][ZMAX];
	//f = cos()

	for(int i=0;i<RMAX;i++)
		for (int j = 0; j < ZMAX; j++)
		{
			res[i][j] = 0;
			test_s[i][j] = 0;
		}

	for (int i = 60; i < 70; i++)
		for (int j = 60; j < 70; j++)
		{
			test_s[i][j] = 1e-16;
		}
	//f = cos();
	matrix_to_csv((float**)test_s, ZMAX, RMAX, RMAX, (char*)(".\\output\\test_s.csv"));
	solve(res, test_s, 1000, 1e-5);
	matrix_to_csv((float**)res, ZMAX, RMAX, RMAX, (char*)(".\\output\\testphi.csv"));
}
