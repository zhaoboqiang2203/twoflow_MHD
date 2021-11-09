#include "fluid.h"
using namespace std;

void electron_flow()
{
	int row = 20, col = 4;

	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nz; j++)
		{
			if (btype[i][j] == 1)
			{

			}
			U[i][j] = cal_u(i, j);
			Fr[i][j] = cal_fr(U[i][j]);
			Fz[i][j] = cal_fz(U[i][j]);
			s[i][j] = cal_s(U[i][j]);
		}
	}

	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nz; j++)
		{
			if (btype[i][j] == 1)
			{
				U_bar[i][j].u1 = U[i][j].u1 - dt / dr * (Fr[i][j + 1].f1 - Fr[i][j].f1) - dt / dz * (Fz[i + 1][j].f1 - Fz[i][j].f1) + dt * s[i][j].f1;
				U_bar[i][j].u2 = U[i][j].u2 - dt / dr * (Fr[i][j + 1].f2 - Fr[i][j].f2) - dt / dz * (Fz[i + 1][j].f2 - Fz[i][j].f2) + dt * s[i][j].f2;
				U_bar[i][j].u3 = U[i][j].u3 - dt / dr * (Fr[i][j + 1].f3 - Fr[i][j].f3) - dt / dz * (Fz[i + 1][j].f3 - Fz[i][j].f3) + dt * s[i][j].f3;
				U_bar[i][j].u4 = U[i][j].u4 - dt / dr * (Fr[i][j + 1].f4 - Fr[i][j].f4) - dt / dz * (Fz[i + 1][j].f4 - Fz[i][j].f4) + dt * s[i][j].f4;
				U_bar[i][j].u5 = U[i][j].u5 - dt / dr * (Fr[i][j + 1].f5 - Fr[i][j].f5) - dt / dz * (Fz[i + 1][j].f5 - Fz[i][j].f5) + dt * s[i][j].f5;
				U_bar[i][j].u6 = U[i][j].u6 - dt / dr * (Fr[i][j + 1].f6 - Fr[i][j].f6) - dt / dz * (Fz[i + 1][j].f6 - Fz[i][j].f6) + dt * s[i][j].f6;
				U_bar[i][j].u7 = U[i][j].u7 - dt / dr * (Fr[i][j + 1].f7 - Fr[i][j].f7) - dt / dz * (Fz[i + 1][j].f7 - Fz[i][j].f7) + dt * s[i][j].f7;
				U_bar[i][j].u8 = U[i][j].u8 - dt / dr * (Fr[i][j + 1].f8 - Fr[i][j].f8) - dt / dz * (Fz[i + 1][j].f8 - Fz[i][j].f8) + dt * s[i][j].f8;
				U_bar[i][j].u9 = U[i][j].u9 - dt / dr * (Fr[i][j + 1].f9 - Fr[i][j].f9) - dt / dz * (Fz[i + 1][j].f9 - Fz[i][j].f9) + dt * s[i][j].f9;
				U_bar[i][j].u10 = U[i][j].u10 - dt / dr * (Fr[i][j + 1].f10 - Fr[i][j].f10) - dt / dz * (Fz[i + 1][j].f10 - Fz[i][j].f10) + dt * s[i][j].f10;
				U_bar[i][j].u11 = U[i][j].u11 - dt / dr * (Fr[i][j + 1].f11 - Fr[i][j].f11) - dt / dz * (Fz[i + 1][j].f11 - Fz[i][j].f11) + dt * s[i][j].f11;
				U_bar[i][j].u12 = U[i][j].u12 - dt / dr * (Fr[i][j + 1].f12 - Fr[i][j].f12) - dt / dz * (Fz[i + 1][j].f12 - Fz[i][j].f12) + dt * s[i][j].f12;
				U_bar[i][j].u13 = U[i][j].u13 - dt / dr * (Fr[i][j + 1].f13 - Fr[i][j].f13) - dt / dz * (Fz[i + 1][j].f13 - Fz[i][j].f13) + dt * s[i][j].f13;


				if (i == row && j == col)
				{
					printf("U[i][j].u1 = %.5e\n", U[i][j].u1);
					printf("Fr[i][j + 1].f1 = %.5e\n", Fr[i][j + 1].f1);
					printf("Fr[i][j].f1 = %.5e\n", Fr[i][j].f1);
					printf("Fz[i + 1][j].f1 = %.5e\n", Fz[i + 1][j].f1);
					printf("Fz[i][j].f1 = %.5e\n", Fz[i][j].f1);
					printf("s[i][j].f1 = %.5e\n", s[i][j].f1);
					printf("U_bar[i][j].u1 = %.5e\n", U_bar[i][j].u1);
					
				}
				//if (U_bar[i][j].u5 != U_bar[i][j].u5)
				//{
				//	printf("U_bar = %.5e\n", U_bar[i][j].u5);
				//	printf("Fz = %.5e\n", Fz[i][j].f5);
				//	printf("Fr = %.5e\n", Fr[i][j].f5);
				//	printf("s = %.5e\n", s[i][j].f5);
				//	printf("U_bar2 = %.5e\n", U_bar2[i][j].u5);
				//	printf("Fz_bar = %.5e\n", Fz_bar[i][j].f5);
				//	printf("Fr_bar = %.5e\n", Fr_bar[i][j].f5);
				//	printf("s_bar = %.5e\n", s_bar[i][j].f5);

				//}
			}
			else if (btype[i][j] == LEFT)//左边的边界，复制右边的参数
			{
				U_bar[i][j].u1 = U_bar[i + 1][j].u1;
				U_bar[i][j].u2 = U_bar[i + 1][j].u2;
				U_bar[i][j].u3 = U_bar[i + 1][j].u3;
				U_bar[i][j].u4 = U_bar[i + 1][j].u4;
				U_bar[i][j].u5 = U_bar[i + 1][j].u5;
				U_bar[i][j].u6 = U_bar[i + 1][j].u6;
				U_bar[i][j].u7 = U_bar[i + 1][j].u7;
				U_bar[i][j].u8 = U_bar[i + 1][j].u8;
				U_bar[i][j].u9 = U_bar[i + 1][j].u9;
				U_bar[i][j].u10 = U_bar[i + 1][j].u10;
				U_bar[i][j].u11 = U_bar[i + 1][j].u11;
				U_bar[i][j].u12 = U_bar[i + 1][j].u12;
				U_bar[i][j].u13 = U_bar[i + 1][j].u13;

			}
			else if (btype[i][j] == RIGHT)
			{
				U_bar[i][j].u1 = U_bar[i - 1][j].u1;
				U_bar[i][j].u2 = U_bar[i - 1][j].u2;
				U_bar[i][j].u3 = U_bar[i - 1][j].u3;
				U_bar[i][j].u4 = U_bar[i - 1][j].u4;
				U_bar[i][j].u5 = U_bar[i - 1][j].u5;
				U_bar[i][j].u6 = U_bar[i - 1][j].u6;
				U_bar[i][j].u7 = U_bar[i - 1][j].u7;
				U_bar[i][j].u8 = U_bar[i - 1][j].u8;
				U_bar[i][j].u9 = U_bar[i - 1][j].u9;
				U_bar[i][j].u10 = U_bar[i - 1][j].u10;
				U_bar[i][j].u11 = U_bar[i - 1][j].u11;
				U_bar[i][j].u12 = U_bar[i - 1][j].u12;
				U_bar[i][j].u13 = U_bar[i - 1][j].u13;
			}
			else if (btype[i][j] == UP)
			{
				U_bar[i][j].u1 = U_bar[i][j - 1].u1;
				U_bar[i][j].u2 = U_bar[i][j - 1].u2;
				U_bar[i][j].u3 = U_bar[i][j - 1].u3;
				U_bar[i][j].u4 = U_bar[i][j - 1].u4;
				U_bar[i][j].u5 = U_bar[i][j - 1].u5;
				U_bar[i][j].u6 = U_bar[i][j - 1].u6;
				U_bar[i][j].u7 = U_bar[i][j - 1].u7;
				U_bar[i][j].u8 = U_bar[i][j - 1].u8;
				U_bar[i][j].u9 = U_bar[i][j - 1].u9;
				U_bar[i][j].u10 = U_bar[i][j - 1].u10;
				U_bar[i][j].u11 = U_bar[i][j - 1].u11;
				U_bar[i][j].u12 = U_bar[i][j - 1].u12;
				U_bar[i][j].u13 = U_bar[i][j - 1].u13;
			}
			else if (btype[i][j] == DOWN)
			{
				U_bar[i][j].u1 = U_bar[i][j + 1].u1;
				U_bar[i][j].u2 = U_bar[i][j + 1].u2;
				U_bar[i][j].u3 = U_bar[i][j + 1].u3;
				U_bar[i][j].u4 = U_bar[i][j + 1].u4;
				U_bar[i][j].u5 = U_bar[i][j + 1].u5;
				U_bar[i][j].u6 = U_bar[i][j + 1].u6;
				U_bar[i][j].u7 = U_bar[i][j + 1].u7;
				U_bar[i][j].u8 = U_bar[i][j + 1].u8;
				U_bar[i][j].u9 = U_bar[i][j + 1].u9;
				U_bar[i][j].u10 = U_bar[i][j + 1].u10;
				U_bar[i][j].u11 = U_bar[i][j + 1].u11;
				U_bar[i][j].u12 = U_bar[i][j + 1].u12;
				U_bar[i][j].u13 = U_bar[i][j + 1].u13;
			}
			else if (btype[i][j] == (LEFT + UP))
			{
				U_bar[i][j].u1 = U_bar[i + 1][j - 1].u1;
				U_bar[i][j].u2 = U_bar[i + 1][j - 1].u2;
				U_bar[i][j].u3 = U_bar[i + 1][j - 1].u3;
				U_bar[i][j].u4 = U_bar[i + 1][j - 1].u4;
				U_bar[i][j].u5 = U_bar[i + 1][j - 1].u5;
				U_bar[i][j].u6 = U_bar[i + 1][j - 1].u6;
				U_bar[i][j].u7 = U_bar[i + 1][j - 1].u7;
				U_bar[i][j].u8 = U_bar[i + 1][j - 1].u8;
				U_bar[i][j].u9 = U_bar[i + 1][j - 1].u9;
				U_bar[i][j].u10 = U_bar[i + 1][j - 1].u10;
				U_bar[i][j].u11 = U_bar[i + 1][j - 1].u11;
				U_bar[i][j].u12 = U_bar[i + 1][j - 1].u12;
				U_bar[i][j].u13 = U_bar[i + 1][j - 1].u13;
			}
			else if (btype[i][j] == (LEFT + DOWN))
			{
				U_bar[i][j].u1 = U_bar[i + 1][j + 1].u1;
				U_bar[i][j].u2 = U_bar[i + 1][j + 1].u2;
				U_bar[i][j].u3 = U_bar[i + 1][j + 1].u3;
				U_bar[i][j].u4 = U_bar[i + 1][j + 1].u4;
				U_bar[i][j].u5 = U_bar[i + 1][j + 1].u5;
				U_bar[i][j].u6 = U_bar[i + 1][j + 1].u6;
				U_bar[i][j].u7 = U_bar[i + 1][j + 1].u7;
				U_bar[i][j].u8 = U_bar[i + 1][j + 1].u8;
				U_bar[i][j].u9 = U_bar[i + 1][j + 1].u9;
				U_bar[i][j].u10 = U_bar[i + 1][j + 1].u10;
				U_bar[i][j].u11 = U_bar[i + 1][j + 1].u11;
				U_bar[i][j].u12 = U_bar[i + 1][j + 1].u12;
				U_bar[i][j].u13 = U_bar[i + 1][j + 1].u13;
			}
			else if (btype[i][j] == (RIGHT + DOWN))
			{
				U_bar[i][j].u1 = U_bar[i - 1][j + 1].u1;
				U_bar[i][j].u2 = U_bar[i - 1][j + 1].u2;
				U_bar[i][j].u3 = U_bar[i - 1][j + 1].u3;
				U_bar[i][j].u4 = U_bar[i - 1][j + 1].u4;
				U_bar[i][j].u5 = U_bar[i - 1][j + 1].u5;
				U_bar[i][j].u6 = U_bar[i - 1][j + 1].u6;
				U_bar[i][j].u7 = U_bar[i - 1][j + 1].u7;
				U_bar[i][j].u8 = U_bar[i - 1][j + 1].u8;
				U_bar[i][j].u9 = U_bar[i - 1][j + 1].u9;
				U_bar[i][j].u10 = U_bar[i - 1][j + 1].u10;
				U_bar[i][j].u11 = U_bar[i - 1][j + 1].u11;
				U_bar[i][j].u12 = U_bar[i - 1][j + 1].u12;
				U_bar[i][j].u13 = U_bar[i - 1][j + 1].u13;
			}
			else if (btype[i][j] == (RIGHT + UP))
			{
				U_bar[i][j].u1 = U_bar[i - 1][j - 1].u1;
				U_bar[i][j].u2 = U_bar[i - 1][j - 1].u2;
				U_bar[i][j].u3 = U_bar[i - 1][j - 1].u3;
				U_bar[i][j].u4 = U_bar[i - 1][j - 1].u4;
				U_bar[i][j].u5 = U_bar[i - 1][j - 1].u5;
				U_bar[i][j].u6 = U_bar[i - 1][j - 1].u6;
				U_bar[i][j].u7 = U_bar[i - 1][j - 1].u7;
				U_bar[i][j].u8 = U_bar[i - 1][j - 1].u8;
				U_bar[i][j].u9 = U_bar[i - 1][j - 1].u9;
				U_bar[i][j].u10 = U_bar[i - 1][j - 1].u10;
				U_bar[i][j].u11 = U_bar[i - 1][j - 1].u11;
				U_bar[i][j].u12 = U_bar[i - 1][j - 1].u12;
				U_bar[i][j].u13 = U_bar[i - 1][j - 1].u13;
			}
			//else if(btype[i][j] == LEFT)//左边的边界，复制右边的参数
			//{

			//	U_bar[i][j].u1 = U[i][j].u1 - dt / dr * (Fr[i][j + 1].f1 - Fr[i][j].f1) - dt / dz * (Fz[i + 1][j].f1 - Fz[i][j].f1) + dt * s[i][j].f1;
			//	U_bar[i][j].u2 = U[i][j].u2 - dt / dr * (Fr[i][j + 1].f2 - Fr[i][j].f2) - dt / dz * (Fz[i + 1][j].f2 - Fz[i][j].f2) + dt * s[i][j].f2;
			//	U_bar[i][j].u3 = U[i][j].u3 - dt / dr * (Fr[i][j + 1].f3 - Fr[i][j].f3) - dt / dz * (Fz[i + 1][j].f3 - Fz[i][j].f3) + dt * s[i][j].f3;
			//	U_bar[i][j].u4 = U[i][j].u4 - dt / dr * (Fr[i][j + 1].f4 - Fr[i][j].f4) - dt / dz * (Fz[i + 1][j].f4 - Fz[i][j].f4) + dt * s[i][j].f4;
			//	U_bar[i][j].u5 = U[i][j].u5 - dt / dr * (Fr[i][j + 1].f5 - Fr[i][j].f5) - dt / dz * (Fz[i + 1][j].f5 - Fz[i][j].f5) + dt * s[i][j].f5;
			//	U_bar[i][j].u6 = U[i][j].u6 - dt / dr * (Fr[i][j + 1].f6 - Fr[i][j].f6) - dt / dz * (Fz[i + 1][j].f6 - Fz[i][j].f6) + dt * s[i][j].f6;
			//	U_bar[i][j].u7 = U[i][j].u7 - dt / dr * (Fr[i][j + 1].f7 - Fr[i][j].f7) - dt / dz * (Fz[i + 1][j].f7 - Fz[i][j].f7) + dt * s[i][j].f7;
			//	U_bar[i][j].u8 = U[i][j].u8 - dt / dr * (Fr[i][j + 1].f8 - Fr[i][j].f8) - dt / dz * (Fz[i + 1][j].f8 - Fz[i][j].f8) + dt * s[i][j].f8;
			//	U_bar[i][j].u9 = U[i][j].u9 - dt / dr * (Fr[i][j + 1].f9 - Fr[i][j].f9) - dt / dz * (Fz[i + 1][j].f9 - Fz[i][j].f9) + dt * s[i][j].f9;
			//	U_bar[i][j].u10 = U[i][j].u10 - dt / dr * (Fr[i][j + 1].f10 - Fr[i][j].f10) - dt / dz * (Fz[i + 1][j].f10 - Fz[i][j].f10) + dt * s[i][j].f10;
			//	U_bar[i][j].u11 = U[i][j].u11 - dt / dr * (Fr[i][j + 1].f11 - Fr[i][j].f11) - dt / dz * (Fz[i + 1][j].f11 - Fz[i][j].f11) + dt * s[i][j].f11;
			//	U_bar[i][j].u12 = U[i][j].u12 - dt / dr * (Fr[i][j + 1].f12 - Fr[i][j].f12) - dt / dz * (Fz[i + 1][j].f12 - Fz[i][j].f12) + dt * s[i][j].f12;
			//	U_bar[i][j].u13 = U[i][j].u13 - dt / dr * (Fr[i][j + 1].f13 - Fr[i][j].f13) - dt / dz * (Fz[i + 1][j].f13 - Fz[i][j].f13) + dt * s[i][j].f13;
			//	
			//}
			//else if(btype[i][j] == RIGHT) 
			//{
			//	
			//	U_bar[i][j].u1 = U[i][j].u1 - dt / dr * (Fr[i][j + 1].f1 - Fr[i][j].f1) - dt / dz * (Fz[i][j].f1 - Fz[i - 1][j].f1) + dt * s[i][j].f1;
			//	U_bar[i][j].u2 = U[i][j].u2 - dt / dr * (Fr[i][j + 1].f2 - Fr[i][j].f2) - dt / dz * (Fz[i][j].f2 - Fz[i - 1][j].f2) + dt * s[i][j].f2;
			//	U_bar[i][j].u3 = U[i][j].u3 - dt / dr * (Fr[i][j + 1].f3 - Fr[i][j].f3) - dt / dz * (Fz[i][j].f3 - Fz[i - 1][j].f3) + dt * s[i][j].f3;
			//	U_bar[i][j].u4 = U[i][j].u4 - dt / dr * (Fr[i][j + 1].f4 - Fr[i][j].f4) - dt / dz * (Fz[i][j].f4 - Fz[i - 1][j].f4) + dt * s[i][j].f4;
			//	U_bar[i][j].u5 = U[i][j].u5 - dt / dr * (Fr[i][j + 1].f5 - Fr[i][j].f5) - dt / dz * (Fz[i][j].f5 - Fz[i - 1][j].f5) + dt * s[i][j].f5;
			//	U_bar[i][j].u6 = U[i][j].u6 - dt / dr * (Fr[i][j + 1].f6 - Fr[i][j].f6) - dt / dz * (Fz[i][j].f6 - Fz[i - 1][j].f6) + dt * s[i][j].f6;
			//	U_bar[i][j].u7 = U[i][j].u7 - dt / dr * (Fr[i][j + 1].f7 - Fr[i][j].f7) - dt / dz * (Fz[i][j].f7 - Fz[i - 1][j].f7) + dt * s[i][j].f7;
			//	U_bar[i][j].u8 = U[i][j].u8 - dt / dr * (Fr[i][j + 1].f8 - Fr[i][j].f8) - dt / dz * (Fz[i][j].f8 - Fz[i - 1][j].f8) + dt * s[i][j].f8;
			//	U_bar[i][j].u9 = U[i][j].u9 - dt / dr * (Fr[i][j + 1].f9 - Fr[i][j].f9) - dt / dz * (Fz[i][j].f9 - Fz[i - 1][j].f9) + dt * s[i][j].f9;
			//	U_bar[i][j].u10 = U[i][j].u10 - dt / dr * (Fr[i][j + 1].f10 - Fr[i][j].f10) - dt / dz * (Fz[i][j].f10 - Fz[i - 1][j].f10) + dt * s[i][j].f10;
			//	U_bar[i][j].u11 = U[i][j].u11 - dt / dr * (Fr[i][j + 1].f11 - Fr[i][j].f11) - dt / dz * (Fz[i][j].f11 - Fz[i - 1][j].f11) + dt * s[i][j].f11;
			//	U_bar[i][j].u12 = U[i][j].u12 - dt / dr * (Fr[i][j + 1].f12 - Fr[i][j].f12) - dt / dz * (Fz[i][j].f12 - Fz[i - 1][j].f12) + dt * s[i][j].f12;
			//	U_bar[i][j].u13 = U[i][j].u13 - dt / dr * (Fr[i][j + 1].f13 - Fr[i][j].f13) - dt / dz * (Fz[i][j].f13 - Fz[i - 1][j].f13) + dt * s[i][j].f13;
			//}
			//else if(btype[i][j] == UP)
			//{

			//	U_bar[i][j].u1 = U[i][j].u1 - dt / dr * (Fr[i][j].f1 - Fr[i][j - 1].f1) - dt / dz * (Fz[i + 1][j].f1 - Fz[i][j].f1) + dt * s[i][j].f1;
			//	U_bar[i][j].u2 = U[i][j].u2 - dt / dr * (Fr[i][j].f2 - Fr[i][j - 1].f2) - dt / dz * (Fz[i + 1][j].f2 - Fz[i][j].f2) + dt * s[i][j].f2;
			//	U_bar[i][j].u3 = U[i][j].u3 - dt / dr * (Fr[i][j].f3 - Fr[i][j - 1].f3) - dt / dz * (Fz[i + 1][j].f3 - Fz[i][j].f3) + dt * s[i][j].f3;
			//	U_bar[i][j].u4 = U[i][j].u4 - dt / dr * (Fr[i][j].f4 - Fr[i][j - 1].f4) - dt / dz * (Fz[i + 1][j].f4 - Fz[i][j].f4) + dt * s[i][j].f4;
			//	U_bar[i][j].u5 = U[i][j].u5 - dt / dr * (Fr[i][j].f5 - Fr[i][j - 1].f5) - dt / dz * (Fz[i + 1][j].f5 - Fz[i][j].f5) + dt * s[i][j].f5;
			//	U_bar[i][j].u6 = U[i][j].u6 - dt / dr * (Fr[i][j].f6 - Fr[i][j - 1].f6) - dt / dz * (Fz[i + 1][j].f6 - Fz[i][j].f6) + dt * s[i][j].f6;
			//	U_bar[i][j].u7 = U[i][j].u7 - dt / dr * (Fr[i][j].f7 - Fr[i][j - 1].f7) - dt / dz * (Fz[i + 1][j].f7 - Fz[i][j].f7) + dt * s[i][j].f7;
			//	U_bar[i][j].u8 = U[i][j].u8 - dt / dr * (Fr[i][j].f8 - Fr[i][j - 1].f8) - dt / dz * (Fz[i + 1][j].f8 - Fz[i][j].f8) + dt * s[i][j].f8;
			//	U_bar[i][j].u9 = U[i][j].u9 - dt / dr * (Fr[i][j].f9 - Fr[i][j - 1].f9) - dt / dz * (Fz[i + 1][j].f9 - Fz[i][j].f9) + dt * s[i][j].f9;
			//	U_bar[i][j].u10 = U[i][j].u10 - dt / dr * (Fr[i][j].f10 - Fr[i][j - 1].f10) - dt / dz * (Fz[i + 1][j].f10 - Fz[i][j].f10) + dt * s[i][j].f10;
			//	U_bar[i][j].u11 = U[i][j].u11 - dt / dr * (Fr[i][j].f11 - Fr[i][j - 1].f11) - dt / dz * (Fz[i + 1][j].f11 - Fz[i][j].f11) + dt * s[i][j].f11;
			//	U_bar[i][j].u12 = U[i][j].u12 - dt / dr * (Fr[i][j].f12 - Fr[i][j - 1].f12) - dt / dz * (Fz[i + 1][j].f12 - Fz[i][j].f12) + dt * s[i][j].f12;
			//	U_bar[i][j].u13 = U[i][j].u13 - dt / dr * (Fr[i][j].f13 - Fr[i][j - 1].f13) - dt / dz * (Fz[i + 1][j].f13 - Fz[i][j].f13) + dt * s[i][j].f13;
			//}
			//else if(btype[i][j] == DOWN)
			//{

			//	U_bar[i][j].u1 = U[i][j].u1 - dt / dr * (Fr[i][j + 1].f1 - Fr[i][j].f1) - dt / dz * (Fz[i + 1][j].f1 - Fz[i][j].f1) + dt * s[i][j].f1;
			//	U_bar[i][j].u2 = U[i][j].u2 - dt / dr * (Fr[i][j + 1].f2 - Fr[i][j].f2) - dt / dz * (Fz[i + 1][j].f2 - Fz[i][j].f2) + dt * s[i][j].f2;
			//	U_bar[i][j].u3 = U[i][j].u3 - dt / dr * (Fr[i][j + 1].f3 - Fr[i][j].f3) - dt / dz * (Fz[i + 1][j].f3 - Fz[i][j].f3) + dt * s[i][j].f3;
			//	U_bar[i][j].u4 = U[i][j].u4 - dt / dr * (Fr[i][j + 1].f4 - Fr[i][j].f4) - dt / dz * (Fz[i + 1][j].f4 - Fz[i][j].f4) + dt * s[i][j].f4;
			//	U_bar[i][j].u5 = U[i][j].u5 - dt / dr * (Fr[i][j + 1].f5 - Fr[i][j].f5) - dt / dz * (Fz[i + 1][j].f5 - Fz[i][j].f5) + dt * s[i][j].f5;
			//	U_bar[i][j].u6 = U[i][j].u6 - dt / dr * (Fr[i][j + 1].f6 - Fr[i][j].f6) - dt / dz * (Fz[i + 1][j].f6 - Fz[i][j].f6) + dt * s[i][j].f6;
			//	U_bar[i][j].u7 = U[i][j].u7 - dt / dr * (Fr[i][j + 1].f7 - Fr[i][j].f7) - dt / dz * (Fz[i + 1][j].f7 - Fz[i][j].f7) + dt * s[i][j].f7;
			//	U_bar[i][j].u8 = U[i][j].u8 - dt / dr * (Fr[i][j + 1].f8 - Fr[i][j].f8) - dt / dz * (Fz[i + 1][j].f8 - Fz[i][j].f8) + dt * s[i][j].f8;
			//	U_bar[i][j].u9 = U[i][j].u9 - dt / dr * (Fr[i][j + 1].f9 - Fr[i][j].f9) - dt / dz * (Fz[i + 1][j].f9 - Fz[i][j].f9) + dt * s[i][j].f9;
			//	U_bar[i][j].u10 = U[i][j].u10 - dt / dr * (Fr[i][j + 1].f10 - Fr[i][j].f10) - dt / dz * (Fz[i + 1][j].f10 - Fz[i][j].f10) + dt * s[i][j].f10;
			//	U_bar[i][j].u11 = U[i][j].u11 - dt / dr * (Fr[i][j + 1].f11 - Fr[i][j].f11) - dt / dz * (Fz[i + 1][j].f11 - Fz[i][j].f11) + dt * s[i][j].f11;
			//	U_bar[i][j].u12 = U[i][j].u12 - dt / dr * (Fr[i][j + 1].f12 - Fr[i][j].f12) - dt / dz * (Fz[i + 1][j].f12 - Fz[i][j].f12) + dt * s[i][j].f12;
			//	U_bar[i][j].u13 = U[i][j].u13 - dt / dr * (Fr[i][j + 1].f13 - Fr[i][j].f13) - dt / dz * (Fz[i + 1][j].f13 - Fz[i][j].f13) + dt * s[i][j].f13;
			//}
			//else if(btype[i][j] == (LEFT + UP))
			//{

			//	U_bar[i][j].u1 = U[i][j].u1 - dt / dr * (Fr[i][j].f1 - Fr[i][j - 1].f1) - dt / dz * (Fz[i + 1][j].f1 - Fz[i][j].f1) + dt * s[i][j].f1;
			//	U_bar[i][j].u2 = U[i][j].u2 - dt / dr * (Fr[i][j].f2 - Fr[i][j - 1].f2) - dt / dz * (Fz[i + 1][j].f2 - Fz[i][j].f2) + dt * s[i][j].f2;
			//	U_bar[i][j].u3 = U[i][j].u3 - dt / dr * (Fr[i][j].f3 - Fr[i][j - 1].f3) - dt / dz * (Fz[i + 1][j].f3 - Fz[i][j].f3) + dt * s[i][j].f3;
			//	U_bar[i][j].u4 = U[i][j].u4 - dt / dr * (Fr[i][j].f4 - Fr[i][j - 1].f4) - dt / dz * (Fz[i + 1][j].f4 - Fz[i][j].f4) + dt * s[i][j].f4;
			//	U_bar[i][j].u5 = U[i][j].u5 - dt / dr * (Fr[i][j].f5 - Fr[i][j - 1].f5) - dt / dz * (Fz[i + 1][j].f5 - Fz[i][j].f5) + dt * s[i][j].f5;
			//	U_bar[i][j].u6 = U[i][j].u6 - dt / dr * (Fr[i][j].f6 - Fr[i][j - 1].f6) - dt / dz * (Fz[i + 1][j].f6 - Fz[i][j].f6) + dt * s[i][j].f6;
			//	U_bar[i][j].u7 = U[i][j].u7 - dt / dr * (Fr[i][j].f7 - Fr[i][j - 1].f7) - dt / dz * (Fz[i + 1][j].f7 - Fz[i][j].f7) + dt * s[i][j].f7;
			//	U_bar[i][j].u8 = U[i][j].u8 - dt / dr * (Fr[i][j].f8 - Fr[i][j - 1].f8) - dt / dz * (Fz[i + 1][j].f8 - Fz[i][j].f8) + dt * s[i][j].f8;
			//	U_bar[i][j].u9 = U[i][j].u9 - dt / dr * (Fr[i][j].f9 - Fr[i][j - 1].f9) - dt / dz * (Fz[i + 1][j].f9 - Fz[i][j].f9) + dt * s[i][j].f9;
			//	U_bar[i][j].u10 = U[i][j].u10 - dt / dr * (Fr[i][j].f10 - Fr[i][j - 1].f10) - dt / dz * (Fz[i + 1][j].f10 - Fz[i][j].f10) + dt * s[i][j].f10;
			//	U_bar[i][j].u11 = U[i][j].u11 - dt / dr * (Fr[i][j].f11 - Fr[i][j - 1].f11) - dt / dz * (Fz[i + 1][j].f11 - Fz[i][j].f11) + dt * s[i][j].f11;
			//	U_bar[i][j].u12 = U[i][j].u12 - dt / dr * (Fr[i][j].f12 - Fr[i][j - 1].f12) - dt / dz * (Fz[i + 1][j].f12 - Fz[i][j].f12) + dt * s[i][j].f12;
			//	U_bar[i][j].u13 = U[i][j].u13 - dt / dr * (Fr[i][j].f13 - Fr[i][j - 1].f13) - dt / dz * (Fz[i + 1][j].f13 - Fz[i][j].f13) + dt * s[i][j].f13;
			//}
			//else if(btype[i][j] == (LEFT + DOWN))
			//{

			//	U_bar[i][j].u1 = U[i][j].u1 - dt / dr * (Fr[i][j + 1].f1 - Fr[i][j].f1) - dt / dz * (Fz[i + 1][j].f1 - Fz[i][j].f1) + dt * s[i][j].f1;
			//	U_bar[i][j].u2 = U[i][j].u2 - dt / dr * (Fr[i][j + 1].f2 - Fr[i][j].f2) - dt / dz * (Fz[i + 1][j].f2 - Fz[i][j].f2) + dt * s[i][j].f2;
			//	U_bar[i][j].u3 = U[i][j].u3 - dt / dr * (Fr[i][j + 1].f3 - Fr[i][j].f3) - dt / dz * (Fz[i + 1][j].f3 - Fz[i][j].f3) + dt * s[i][j].f3;
			//	U_bar[i][j].u4 = U[i][j].u4 - dt / dr * (Fr[i][j + 1].f4 - Fr[i][j].f4) - dt / dz * (Fz[i + 1][j].f4 - Fz[i][j].f4) + dt * s[i][j].f4;
			//	U_bar[i][j].u5 = U[i][j].u5 - dt / dr * (Fr[i][j + 1].f5 - Fr[i][j].f5) - dt / dz * (Fz[i + 1][j].f5 - Fz[i][j].f5) + dt * s[i][j].f5;
			//	U_bar[i][j].u6 = U[i][j].u6 - dt / dr * (Fr[i][j + 1].f6 - Fr[i][j].f6) - dt / dz * (Fz[i + 1][j].f6 - Fz[i][j].f6) + dt * s[i][j].f6;
			//	U_bar[i][j].u7 = U[i][j].u7 - dt / dr * (Fr[i][j + 1].f7 - Fr[i][j].f7) - dt / dz * (Fz[i + 1][j].f7 - Fz[i][j].f7) + dt * s[i][j].f7;
			//	U_bar[i][j].u8 = U[i][j].u8 - dt / dr * (Fr[i][j + 1].f8 - Fr[i][j].f8) - dt / dz * (Fz[i + 1][j].f8 - Fz[i][j].f8) + dt * s[i][j].f8;
			//	U_bar[i][j].u9 = U[i][j].u9 - dt / dr * (Fr[i][j + 1].f9 - Fr[i][j].f9) - dt / dz * (Fz[i + 1][j].f9 - Fz[i][j].f9) + dt * s[i][j].f9;
			//	U_bar[i][j].u10 = U[i][j].u10 - dt / dr * (Fr[i][j + 1].f10 - Fr[i][j].f10) - dt / dz * (Fz[i + 1][j].f10 - Fz[i][j].f10) + dt * s[i][j].f10;
			//	U_bar[i][j].u11 = U[i][j].u11 - dt / dr * (Fr[i][j + 1].f11 - Fr[i][j].f11) - dt / dz * (Fz[i + 1][j].f11 - Fz[i][j].f11) + dt * s[i][j].f11;
			//	U_bar[i][j].u12 = U[i][j].u12 - dt / dr * (Fr[i][j + 1].f12 - Fr[i][j].f12) - dt / dz * (Fz[i + 1][j].f12 - Fz[i][j].f12) + dt * s[i][j].f12;
			//	U_bar[i][j].u13 = U[i][j].u13 - dt / dr * (Fr[i][j + 1].f13 - Fr[i][j].f13) - dt / dz * (Fz[i + 1][j].f13 - Fz[i][j].f13) + dt * s[i][j].f13;
			//}
			//else if(btype[i][j] == (RIGHT + DOWN))
			//{

			//	U_bar[i][j].u1 = U[i][j].u1 - dt / dr * (Fr[i][j + 1].f1 - Fr[i][j].f1) - dt / dz * (Fz[i][j].f1 - Fz[i - 1][j].f1) + dt * s[i][j].f1;
			//	U_bar[i][j].u2 = U[i][j].u2 - dt / dr * (Fr[i][j + 1].f2 - Fr[i][j].f2) - dt / dz * (Fz[i][j].f2 - Fz[i - 1][j].f2) + dt * s[i][j].f2;
			//	U_bar[i][j].u3 = U[i][j].u3 - dt / dr * (Fr[i][j + 1].f3 - Fr[i][j].f3) - dt / dz * (Fz[i][j].f3 - Fz[i - 1][j].f3) + dt * s[i][j].f3;
			//	U_bar[i][j].u4 = U[i][j].u4 - dt / dr * (Fr[i][j + 1].f4 - Fr[i][j].f4) - dt / dz * (Fz[i][j].f4 - Fz[i - 1][j].f4) + dt * s[i][j].f4;
			//	U_bar[i][j].u5 = U[i][j].u5 - dt / dr * (Fr[i][j + 1].f5 - Fr[i][j].f5) - dt / dz * (Fz[i][j].f5 - Fz[i - 1][j].f5) + dt * s[i][j].f5;
			//	U_bar[i][j].u6 = U[i][j].u6 - dt / dr * (Fr[i][j + 1].f6 - Fr[i][j].f6) - dt / dz * (Fz[i][j].f6 - Fz[i - 1][j].f6) + dt * s[i][j].f6;
			//	U_bar[i][j].u7 = U[i][j].u7 - dt / dr * (Fr[i][j + 1].f7 - Fr[i][j].f7) - dt / dz * (Fz[i][j].f7 - Fz[i - 1][j].f7) + dt * s[i][j].f7;
			//	U_bar[i][j].u8 = U[i][j].u8 - dt / dr * (Fr[i][j + 1].f8 - Fr[i][j].f8) - dt / dz * (Fz[i][j].f8 - Fz[i - 1][j].f8) + dt * s[i][j].f8;
			//	U_bar[i][j].u9 = U[i][j].u9 - dt / dr * (Fr[i][j + 1].f9 - Fr[i][j].f9) - dt / dz * (Fz[i][j].f9 - Fz[i - 1][j].f9) + dt * s[i][j].f9;
			//	U_bar[i][j].u10 = U[i][j].u10 - dt / dr * (Fr[i][j + 1].f10 - Fr[i][j].f10) - dt / dz * (Fz[i][j].f10 - Fz[i - 1][j].f10) + dt * s[i][j].f10;
			//	U_bar[i][j].u11 = U[i][j].u11 - dt / dr * (Fr[i][j + 1].f11 - Fr[i][j].f11) - dt / dz * (Fz[i][j].f11 - Fz[i - 1][j].f11) + dt * s[i][j].f11;
			//	U_bar[i][j].u12 = U[i][j].u12 - dt / dr * (Fr[i][j + 1].f12 - Fr[i][j].f12) - dt / dz * (Fz[i][j].f12 - Fz[i - 1][j].f12) + dt * s[i][j].f12;
			//	U_bar[i][j].u13 = U[i][j].u13 - dt / dr * (Fr[i][j + 1].f13 - Fr[i][j].f13) - dt / dz * (Fz[i][j].f13 - Fz[i - 1][j].f13) + dt * s[i][j].f13;
			//}
			//else if(btype[i][j] == (RIGHT + UP))
			//{

			//	U_bar[i][j].u1 = U[i][j].u1 - dt / dr * (Fr[i][j].f1 - Fr[i][j - 1].f1) - dt / dz * (Fz[i][j].f1 - Fz[i - 1][j].f1) + dt * s[i][j].f1;
			//	U_bar[i][j].u2 = U[i][j].u2 - dt / dr * (Fr[i][j].f2 - Fr[i][j - 1].f2) - dt / dz * (Fz[i][j].f2 - Fz[i - 1][j].f2) + dt * s[i][j].f2;
			//	U_bar[i][j].u3 = U[i][j].u3 - dt / dr * (Fr[i][j].f3 - Fr[i][j - 1].f3) - dt / dz * (Fz[i][j].f3 - Fz[i - 1][j].f3) + dt * s[i][j].f3;
			//	U_bar[i][j].u4 = U[i][j].u4 - dt / dr * (Fr[i][j].f4 - Fr[i][j - 1].f4) - dt / dz * (Fz[i][j].f4 - Fz[i - 1][j].f4) + dt * s[i][j].f4;
			//	U_bar[i][j].u5 = U[i][j].u5 - dt / dr * (Fr[i][j].f5 - Fr[i][j - 1].f5) - dt / dz * (Fz[i][j].f5 - Fz[i - 1][j].f5) + dt * s[i][j].f5;
			//	U_bar[i][j].u6 = U[i][j].u6 - dt / dr * (Fr[i][j].f6 - Fr[i][j - 1].f6) - dt / dz * (Fz[i][j].f6 - Fz[i - 1][j].f6) + dt * s[i][j].f6;
			//	U_bar[i][j].u7 = U[i][j].u7 - dt / dr * (Fr[i][j].f7 - Fr[i][j - 1].f7) - dt / dz * (Fz[i][j].f7 - Fz[i - 1][j].f7) + dt * s[i][j].f7;
			//	U_bar[i][j].u8 = U[i][j].u8 - dt / dr * (Fr[i][j].f8 - Fr[i][j - 1].f8) - dt / dz * (Fz[i][j].f8 - Fz[i - 1][j].f8) + dt * s[i][j].f8;
			//	U_bar[i][j].u9 = U[i][j].u9 - dt / dr * (Fr[i][j].f9 - Fr[i][j - 1].f9) - dt / dz * (Fz[i][j].f9 - Fz[i - 1][j].f9) + dt * s[i][j].f9;
			//	U_bar[i][j].u10 = U[i][j].u10 - dt / dr * (Fr[i][j].f10 - Fr[i][j - 1].f10) - dt / dz * (Fz[i][j].f10 - Fz[i - 1][j].f10) + dt * s[i][j].f10;
			//	U_bar[i][j].u11 = U[i][j].u11 - dt / dr * (Fr[i][j].f11 - Fr[i][j - 1].f11) - dt / dz * (Fz[i][j].f11 - Fz[i - 1][j].f11) + dt * s[i][j].f11;
			//	U_bar[i][j].u12 = U[i][j].u12 - dt / dr * (Fr[i][j].f12 - Fr[i][j - 1].f12) - dt / dz * (Fz[i][j].f12 - Fz[i - 1][j].f12) + dt * s[i][j].f12;
			//	U_bar[i][j].u13 = U[i][j].u13 - dt / dr * (Fr[i][j].f13 - Fr[i][j - 1].f13) - dt / dz * (Fz[i][j].f13 - Fz[i - 1][j].f13) + dt * s[i][j].f13;
			//}
			else //if(btype[i][j] == 0)
			{
				U_bar[i][j].u1 = 0;
				U_bar[i][j].u2 = 0;
				U_bar[i][j].u3 = 0;
				U_bar[i][j].u4 = 0;
				U_bar[i][j].u5 = 0;
				U_bar[i][j].u6 = 0;
				U_bar[i][j].u7 = 0;
				U_bar[i][j].u8 = 0;
				U_bar[i][j].u9 = 0;
				U_bar[i][j].u10 = 0;
				U_bar[i][j].u11 = 0;
				U_bar[i][j].u12 = 0;
				U_bar[i][j].u13 = 0;


			}

			Fr_bar[i][j] = cal_fr(U_bar[i][j]);
			Fz_bar[i][j] = cal_fz(U_bar[i][j]);
			s_bar[i][j] = cal_s(U_bar[i][j]);
		}
	}

	//矫正步
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nz; j++)
		{
			if (btype[i][j] == 1)
			{
				U_bar2[i][j].u1 = 0.5 * (U[i][j].u1 + U_bar[i][j].u1 - dt / dr * (Fr_bar[i][j].f1 - Fr_bar[i][j - 1].f1) - dt / dz * (Fz_bar[i][j].f1 - Fz_bar[i - 1][j].f1) + dt * s_bar[i][j].f1);
				U_bar2[i][j].u2 = 0.5 * (U[i][j].u2 + U_bar[i][j].u2 - dt / dr * (Fr_bar[i][j].f2 - Fr_bar[i][j - 1].f2) - dt / dz * (Fz_bar[i][j].f2 - Fz_bar[i - 1][j].f2) + dt * s_bar[i][j].f2);
				U_bar2[i][j].u3 = 0.5 * (U[i][j].u3 + U_bar[i][j].u3 - dt / dr * (Fr_bar[i][j].f3 - Fr_bar[i][j - 1].f3) - dt / dz * (Fz_bar[i][j].f3 - Fz_bar[i - 1][j].f3) + dt * s_bar[i][j].f3);
				U_bar2[i][j].u4 = 0.5 * (U[i][j].u4 + U_bar[i][j].u4 - dt / dr * (Fr_bar[i][j].f4 - Fr_bar[i][j - 1].f4) - dt / dz * (Fz_bar[i][j].f4 - Fz_bar[i - 1][j].f4) + dt * s_bar[i][j].f4);
				U_bar2[i][j].u5 = 0.5 * (U[i][j].u5 + U_bar[i][j].u5 - dt / dr * (Fr_bar[i][j].f5 - Fr_bar[i][j - 1].f5) - dt / dz * (Fz_bar[i][j].f5 - Fz_bar[i - 1][j].f5) + dt * s_bar[i][j].f5);
				U_bar2[i][j].u6 = 0.5 * (U[i][j].u6 + U_bar[i][j].u6 - dt / dr * (Fr_bar[i][j].f6 - Fr_bar[i][j - 1].f6) - dt / dz * (Fz_bar[i][j].f6 - Fz_bar[i - 1][j].f6) + dt * s_bar[i][j].f6);
				U_bar2[i][j].u7 = 0.5 * (U[i][j].u7 + U_bar[i][j].u7 - dt / dr * (Fr_bar[i][j].f7 - Fr_bar[i][j - 1].f7) - dt / dz * (Fz_bar[i][j].f7 - Fz_bar[i - 1][j].f7) + dt * s_bar[i][j].f7);
				U_bar2[i][j].u8 = 0.5 * (U[i][j].u8 + U_bar[i][j].u8 - dt / dr * (Fr_bar[i][j].f8 - Fr_bar[i][j - 1].f8) - dt / dz * (Fz_bar[i][j].f8 - Fz_bar[i - 1][j].f8) + dt * s_bar[i][j].f8);
				U_bar2[i][j].u9 = 0.5 * (U[i][j].u9 + U_bar[i][j].u9 - dt / dr * (Fr_bar[i][j].f9 - Fr_bar[i][j - 1].f9) - dt / dz * (Fz_bar[i][j].f9 - Fz_bar[i - 1][j].f9) + dt * s_bar[i][j].f9);
				U_bar2[i][j].u10 = 0.5 * (U[i][j].u10 + U_bar[i][j].u10 - dt / dr * (Fr_bar[i][j].f10 - Fr_bar[i][j - 1].f10) - dt / dz * (Fz_bar[i][j].f10 - Fz_bar[i - 1][j].f10) + dt * s_bar[i][j].f10);
				U_bar2[i][j].u11 = 0.5 * (U[i][j].u11 + U_bar[i][j].u11 - dt / dr * (Fr_bar[i][j].f11 - Fr_bar[i][j - 1].f11) - dt / dz * (Fz_bar[i][j].f11 - Fz_bar[i - 1][j].f11) + dt * s_bar[i][j].f11);
				U_bar2[i][j].u12 = 0.5 * (U[i][j].u12 + U_bar[i][j].u12 - dt / dr * (Fr_bar[i][j].f12 - Fr_bar[i][j - 1].f12) - dt / dz * (Fz_bar[i][j].f12 - Fz_bar[i - 1][j].f12) + dt * s_bar[i][j].f12);
				U_bar2[i][j].u13 = 0.5 * (U[i][j].u13 + U_bar[i][j].u13 - dt / dr * (Fr_bar[i][j].f13 - Fr_bar[i][j - 1].f13) - dt / dz * (Fz_bar[i][j].f13 - Fz_bar[i - 1][j].f13) + dt * s_bar[i][j].f13);
				
				if (i == row && j == col)
				{
					printf("U[i][j].u1 = %.5e\n", U[i][j].u1);
					printf("Fr_bar[i][j].f1 = %.5e\n", Fr_bar[i][j].f1);
					printf("Fr_bar[i][j - 1].f1 = %.5e\n", Fr_bar[i][j - 1].f1);
					printf("Fz_bar[i][j].f1 = %.5e\n", Fz_bar[i][j].f1);
					printf("Fz_bar[i - 1][j].f1 = %.5e\n", Fz_bar[i - 1][j].f1);
					printf("s_bar[i][j].f1 = %.5e\n", s_bar[i][j].f1);
					printf("U_bar2[i][j].u1 = %.5e\n", U_bar2[i][j].u1);

				}
			}
			else if (btype[i][j] == LEFT)//左边的边界，复制右边的参数
			{
				U_bar2[i][j].u1 = U_bar2[i + 1][j].u1;
				U_bar2[i][j].u2 = U_bar2[i + 1][j].u2;
				U_bar2[i][j].u3 = U_bar2[i + 1][j].u3;
				U_bar2[i][j].u4 = U_bar2[i + 1][j].u4;
				U_bar2[i][j].u5 = U_bar2[i + 1][j].u5;
				U_bar2[i][j].u6 = U_bar2[i + 1][j].u6;
				U_bar2[i][j].u7 = U_bar2[i + 1][j].u7;
				U_bar2[i][j].u8 = U_bar2[i + 1][j].u8;
				U_bar2[i][j].u9 = U_bar2[i + 1][j].u9;
				U_bar2[i][j].u10 = U_bar2[i + 1][j].u10;
				U_bar2[i][j].u11 = U_bar2[i + 1][j].u11;
				U_bar2[i][j].u12 = U_bar2[i + 1][j].u12;
				U_bar2[i][j].u13 = U_bar2[i + 1][j].u13;

			}
			else if (btype[i][j] == RIGHT)
			{
				U_bar2[i][j].u1 = U_bar2[i - 1][j].u1;
				U_bar2[i][j].u2 = U_bar2[i - 1][j].u2;
				U_bar2[i][j].u3 = U_bar2[i - 1][j].u3;
				U_bar2[i][j].u4 = U_bar2[i - 1][j].u4;
				U_bar2[i][j].u5 = U_bar2[i - 1][j].u5;
				U_bar2[i][j].u6 = U_bar2[i - 1][j].u6;
				U_bar2[i][j].u7 = U_bar2[i - 1][j].u7;
				U_bar2[i][j].u8 = U_bar2[i - 1][j].u8;
				U_bar2[i][j].u9 = U_bar2[i - 1][j].u9;
				U_bar2[i][j].u10 = U_bar2[i - 1][j].u10;
				U_bar2[i][j].u11 = U_bar2[i - 1][j].u11;
				U_bar2[i][j].u12 = U_bar2[i - 1][j].u12;
				U_bar2[i][j].u13 = U_bar2[i - 1][j].u13;
			}
			else if (btype[i][j] == UP)
			{
				U_bar2[i][j].u1 = U_bar2[i][j - 1].u1;
				U_bar2[i][j].u2 = U_bar2[i][j - 1].u2;
				U_bar2[i][j].u3 = U_bar2[i][j - 1].u3;
				U_bar2[i][j].u4 = U_bar2[i][j - 1].u4;
				U_bar2[i][j].u5 = U_bar2[i][j - 1].u5;
				U_bar2[i][j].u6 = U_bar2[i][j - 1].u6;
				U_bar2[i][j].u7 = U_bar2[i][j - 1].u7;
				U_bar2[i][j].u8 = U_bar2[i][j - 1].u8;
				U_bar2[i][j].u9 = U_bar2[i][j - 1].u9;
				U_bar2[i][j].u10 = U_bar2[i][j - 1].u10;
				U_bar2[i][j].u11 = U_bar2[i][j - 1].u11;
				U_bar2[i][j].u12 = U_bar2[i][j - 1].u12;
				U_bar2[i][j].u13 = U_bar2[i][j - 1].u13;
			}
			else if (btype[i][j] == DOWN)
			{
				U_bar2[i][j].u1 = U_bar2[i][j + 1].u1;
				U_bar2[i][j].u2 = U_bar2[i][j + 1].u2;
				U_bar2[i][j].u3 = U_bar2[i][j + 1].u3;
				U_bar2[i][j].u4 = U_bar2[i][j + 1].u4;
				U_bar2[i][j].u5 = U_bar2[i][j + 1].u5;
				U_bar2[i][j].u6 = U_bar2[i][j + 1].u6;
				U_bar2[i][j].u7 = U_bar2[i][j + 1].u7;
				U_bar2[i][j].u8 = U_bar2[i][j + 1].u8;
				U_bar2[i][j].u9 = U_bar2[i][j + 1].u9;
				U_bar2[i][j].u10 = U_bar2[i][j + 1].u10;
				U_bar2[i][j].u11 = U_bar2[i][j + 1].u11;
				U_bar2[i][j].u12 = U_bar2[i][j + 1].u12;
				U_bar2[i][j].u13 = U_bar2[i][j + 1].u13;
			}
			else if (btype[i][j] == (LEFT + UP))
			{
				U_bar2[i][j].u1 = U_bar2[i + 1][j - 1].u1;
				U_bar2[i][j].u2 = U_bar2[i + 1][j - 1].u2;
				U_bar2[i][j].u3 = U_bar2[i + 1][j - 1].u3;
				U_bar2[i][j].u4 = U_bar2[i + 1][j - 1].u4;
				U_bar2[i][j].u5 = U_bar2[i + 1][j - 1].u5;
				U_bar2[i][j].u6 = U_bar2[i + 1][j - 1].u6;
				U_bar2[i][j].u7 = U_bar2[i + 1][j - 1].u7;
				U_bar2[i][j].u8 = U_bar2[i + 1][j - 1].u8;
				U_bar2[i][j].u9 = U_bar2[i + 1][j - 1].u9;
				U_bar2[i][j].u10 = U_bar2[i + 1][j - 1].u10;
				U_bar2[i][j].u11 = U_bar2[i + 1][j - 1].u11;
				U_bar2[i][j].u12 = U_bar2[i + 1][j - 1].u12;
				U_bar2[i][j].u13 = U_bar2[i + 1][j - 1].u13;
			}
			else if (btype[i][j] == (LEFT + DOWN))
			{
				U_bar2[i][j].u1 = U_bar2[i + 1][j + 1].u1;
				U_bar2[i][j].u2 = U_bar2[i + 1][j + 1].u2;
				U_bar2[i][j].u3 = U_bar2[i + 1][j + 1].u3;
				U_bar2[i][j].u4 = U_bar2[i + 1][j + 1].u4;
				U_bar2[i][j].u5 = U_bar2[i + 1][j + 1].u5;
				U_bar2[i][j].u6 = U_bar2[i + 1][j + 1].u6;
				U_bar2[i][j].u7 = U_bar2[i + 1][j + 1].u7;
				U_bar2[i][j].u8 = U_bar2[i + 1][j + 1].u8;
				U_bar2[i][j].u9 = U_bar2[i + 1][j + 1].u9;
				U_bar2[i][j].u10 = U_bar2[i + 1][j + 1].u10;
				U_bar2[i][j].u11 = U_bar2[i + 1][j + 1].u11;
				U_bar2[i][j].u12 = U_bar2[i + 1][j + 1].u12;
				U_bar2[i][j].u13 = U_bar2[i + 1][j + 1].u13;
			}
			else if (btype[i][j] == (RIGHT + DOWN))
			{
				U_bar2[i][j].u1 = U_bar2[i - 1][j + 1].u1;
				U_bar2[i][j].u2 = U_bar2[i - 1][j + 1].u2;
				U_bar2[i][j].u3 = U_bar2[i - 1][j + 1].u3;
				U_bar2[i][j].u4 = U_bar2[i - 1][j + 1].u4;
				U_bar2[i][j].u5 = U_bar2[i - 1][j + 1].u5;
				U_bar2[i][j].u6 = U_bar2[i - 1][j + 1].u6;
				U_bar2[i][j].u7 = U_bar2[i - 1][j + 1].u7;
				U_bar2[i][j].u8 = U_bar2[i - 1][j + 1].u8;
				U_bar2[i][j].u9 = U_bar2[i - 1][j + 1].u9;
				U_bar2[i][j].u10 = U_bar2[i - 1][j + 1].u10;
				U_bar2[i][j].u11 = U_bar2[i - 1][j + 1].u11;
				U_bar2[i][j].u12 = U_bar2[i - 1][j + 1].u12;
				U_bar2[i][j].u13 = U_bar2[i - 1][j + 1].u13;
			}
			else if (btype[i][j] == (RIGHT + UP))
			{
				U_bar2[i][j].u1 = U_bar2[i - 1][j - 1].u1;
				U_bar2[i][j].u2 = U_bar2[i - 1][j - 1].u2;
				U_bar2[i][j].u3 = U_bar2[i - 1][j - 1].u3;
				U_bar2[i][j].u4 = U_bar2[i - 1][j - 1].u4;
				U_bar2[i][j].u5 = U_bar2[i - 1][j - 1].u5;
				U_bar2[i][j].u6 = U_bar2[i - 1][j - 1].u6;
				U_bar2[i][j].u7 = U_bar2[i - 1][j - 1].u7;
				U_bar2[i][j].u8 = U_bar2[i - 1][j - 1].u8;
				U_bar2[i][j].u9 = U_bar2[i - 1][j - 1].u9;
				U_bar2[i][j].u10 = U_bar2[i - 1][j - 1].u10;
				U_bar2[i][j].u11 = U_bar2[i - 1][j - 1].u11;
				U_bar2[i][j].u12 = U_bar2[i - 1][j - 1].u12;
				U_bar2[i][j].u13 = U_bar2[i - 1][j - 1].u13;
			}
			//else if(btype[i][j] == LEFT)//左边的边界，复制右边的参数
			//{
			//	
			//	U_bar2[i][j].u1 = 0.5 * (U[i][j].u1 + U_bar[i][j].u1 - dt / dr * (Fr_bar[i][j].f1 - Fr_bar[i][j - 1].f1) - dt / dz * (Fz_bar[i + 1][j].f1 - Fz_bar[i][j].f1) + dt * s_bar[i][j].f1);
			//	U_bar2[i][j].u2 = 0.5 * (U[i][j].u2 + U_bar[i][j].u2 - dt / dr * (Fr_bar[i][j].f2 - Fr_bar[i][j - 1].f2) - dt / dz * (Fz_bar[i + 1][j].f2 - Fz_bar[i][j].f2) + dt * s_bar[i][j].f2);
			//	U_bar2[i][j].u3 = 0.5 * (U[i][j].u3 + U_bar[i][j].u3 - dt / dr * (Fr_bar[i][j].f3 - Fr_bar[i][j - 1].f3) - dt / dz * (Fz_bar[i + 1][j].f3 - Fz_bar[i][j].f3) + dt * s_bar[i][j].f3);
			//	U_bar2[i][j].u4 = 0.5 * (U[i][j].u4 + U_bar[i][j].u4 - dt / dr * (Fr_bar[i][j].f4 - Fr_bar[i][j - 1].f4) - dt / dz * (Fz_bar[i + 1][j].f4 - Fz_bar[i][j].f4) + dt * s_bar[i][j].f4);
			//	U_bar2[i][j].u5 = 0.5 * (U[i][j].u5 + U_bar[i][j].u5 - dt / dr * (Fr_bar[i][j].f5 - Fr_bar[i][j - 1].f5) - dt / dz * (Fz_bar[i + 1][j].f5 - Fz_bar[i][j].f5) + dt * s_bar[i][j].f5);
			//	U_bar2[i][j].u6 = 0.5 * (U[i][j].u6 + U_bar[i][j].u6 - dt / dr * (Fr_bar[i][j].f6 - Fr_bar[i][j - 1].f6) - dt / dz * (Fz_bar[i + 1][j].f6 - Fz_bar[i][j].f6) + dt * s_bar[i][j].f6);
			//	U_bar2[i][j].u7 = 0.5 * (U[i][j].u7 + U_bar[i][j].u7 - dt / dr * (Fr_bar[i][j].f7 - Fr_bar[i][j - 1].f7) - dt / dz * (Fz_bar[i + 1][j].f7 - Fz_bar[i][j].f7) + dt * s_bar[i][j].f7);
			//	U_bar2[i][j].u8 = 0.5 * (U[i][j].u8 + U_bar[i][j].u8 - dt / dr * (Fr_bar[i][j].f8 - Fr_bar[i][j - 1].f8) - dt / dz * (Fz_bar[i + 1][j].f8 - Fz_bar[i][j].f8) + dt * s_bar[i][j].f8);
			//	U_bar2[i][j].u9 = 0.5 * (U[i][j].u9 + U_bar[i][j].u9 - dt / dr * (Fr_bar[i][j].f9 - Fr_bar[i][j - 1].f9) - dt / dz * (Fz_bar[i + 1][j].f9 - Fz_bar[i][j].f9) + dt * s_bar[i][j].f9);
			//	U_bar2[i][j].u10 = 0.5 * (U[i][j].u10 + U_bar[i][j].u10  - dt / dr * (Fr_bar[i][j].f10 - Fr_bar[i][j - 1].f10) - dt / dz * (Fz_bar[i + 1][j].f10 - Fz_bar[i][j].f10) + dt * s_bar[i][j].f10);
			//	U_bar2[i][j].u11 = 0.5 * (U[i][j].u11 + U_bar[i][j].u11  - dt / dr * (Fr_bar[i][j].f11 - Fr_bar[i][j - 1].f11) - dt / dz * (Fz_bar[i + 1][j].f11 - Fz_bar[i][j].f11) + dt * s_bar[i][j].f11);
			//	U_bar2[i][j].u12 = 0.5 * (U[i][j].u12 + U_bar[i][j].u12  - dt / dr * (Fr_bar[i][j].f12 - Fr_bar[i][j - 1].f12) - dt / dz * (Fz_bar[i + 1][j].f12 - Fz_bar[i][j].f12) + dt * s_bar[i][j].f12);
			//	U_bar2[i][j].u13 = 0.5 * (U[i][j].u13 + U_bar[i][j].u13  - dt / dr * (Fr_bar[i][j].f13 - Fr_bar[i][j - 1].f13) - dt / dz * (Fz_bar[i + 1][j].f13 - Fz_bar[i][j].f13) + dt * s_bar[i][j].f13);

			//}
			//else if(btype[i][j] == RIGHT) 
			//{

			//	U_bar2[i][j].u1 = 0.5 * (U[i][j].u1 + U_bar[i][j].u1 - dt / dr * (Fr_bar[i][j].f1 - Fr_bar[i][j - 1].f1) - dt / dz * (Fz_bar[i][j].f1 - Fz_bar[i - 1][j].f1) + dt * s_bar[i][j].f1);
			//	U_bar2[i][j].u2 = 0.5 * (U[i][j].u2 + U_bar[i][j].u2 - dt / dr * (Fr_bar[i][j].f2 - Fr_bar[i][j - 1].f2) - dt / dz * (Fz_bar[i][j].f2 - Fz_bar[i - 1][j].f2) + dt * s_bar[i][j].f2);
			//	U_bar2[i][j].u3 = 0.5 * (U[i][j].u3 + U_bar[i][j].u3 - dt / dr * (Fr_bar[i][j].f3 - Fr_bar[i][j - 1].f3) - dt / dz * (Fz_bar[i][j].f3 - Fz_bar[i - 1][j].f3) + dt * s_bar[i][j].f3);
			//	U_bar2[i][j].u4 = 0.5 * (U[i][j].u4 + U_bar[i][j].u4 - dt / dr * (Fr_bar[i][j].f4 - Fr_bar[i][j - 1].f4) - dt / dz * (Fz_bar[i][j].f4 - Fz_bar[i - 1][j].f4) + dt * s_bar[i][j].f4);
			//	U_bar2[i][j].u5 = 0.5 * (U[i][j].u5 + U_bar[i][j].u5 - dt / dr * (Fr_bar[i][j].f5 - Fr_bar[i][j - 1].f5) - dt / dz * (Fz_bar[i][j].f5 - Fz_bar[i - 1][j].f5) + dt * s_bar[i][j].f5);
			//	U_bar2[i][j].u6 = 0.5 * (U[i][j].u6 + U_bar[i][j].u6 - dt / dr * (Fr_bar[i][j].f6 - Fr_bar[i][j - 1].f6) - dt / dz * (Fz_bar[i][j].f6 - Fz_bar[i - 1][j].f6) + dt * s_bar[i][j].f6);
			//	U_bar2[i][j].u7 = 0.5 * (U[i][j].u7 + U_bar[i][j].u7 - dt / dr * (Fr_bar[i][j].f7 - Fr_bar[i][j - 1].f7) - dt / dz * (Fz_bar[i][j].f7 - Fz_bar[i - 1][j].f7) + dt * s_bar[i][j].f7);
			//	U_bar2[i][j].u8 = 0.5 * (U[i][j].u8 + U_bar[i][j].u8 - dt / dr * (Fr_bar[i][j].f8 - Fr_bar[i][j - 1].f8) - dt / dz * (Fz_bar[i][j].f8 - Fz_bar[i - 1][j].f8) + dt * s_bar[i][j].f8);
			//	U_bar2[i][j].u9 = 0.5 * (U[i][j].u9 + U_bar[i][j].u9 - dt / dr * (Fr_bar[i][j].f9 - Fr_bar[i][j - 1].f9) - dt / dz * (Fz_bar[i][j].f9 - Fz_bar[i - 1][j].f9) + dt * s_bar[i][j].f9);
			//	U_bar2[i][j].u10 = 0.5 * (U[i][j].u10 + U_bar[i][j].u10  - dt / dr * (Fr_bar[i][j].f10 - Fr_bar[i][j - 1].f10) - dt / dz * (Fz_bar[i][j].f10 - Fz_bar[i - 1][j].f10) + dt * s_bar[i][j].f10);
			//	U_bar2[i][j].u11 = 0.5 * (U[i][j].u11 + U_bar[i][j].u11  - dt / dr * (Fr_bar[i][j].f11 - Fr_bar[i][j - 1].f11) - dt / dz * (Fz_bar[i][j].f11 - Fz_bar[i - 1][j].f11) + dt * s_bar[i][j].f11);
			//	U_bar2[i][j].u12 = 0.5 * (U[i][j].u12 + U_bar[i][j].u12  - dt / dr * (Fr_bar[i][j].f12 - Fr_bar[i][j - 1].f12) - dt / dz * (Fz_bar[i][j].f12 - Fz_bar[i - 1][j].f12) + dt * s_bar[i][j].f12);
			//	U_bar2[i][j].u13 = 0.5 * (U[i][j].u13 + U_bar[i][j].u13  - dt / dr * (Fr_bar[i][j].f13 - Fr_bar[i][j - 1].f13) - dt / dz * (Fz_bar[i][j].f13 - Fz_bar[i - 1][j].f13) + dt * s_bar[i][j].f13);
			//}
			//else if(btype[i][j] == UP)
			//{

			//	U_bar2[i][j].u1 = 0.5 * (U[i][j].u1 + U_bar[i][j].u1 - dt / dr * (Fr_bar[i][j].f1 - Fr_bar[i][j - 1].f1) - dt / dz * (Fz_bar[i][j].f1 - Fz_bar[i - 1][j].f1) + dt * s_bar[i][j].f1);
			//	U_bar2[i][j].u2 = 0.5 * (U[i][j].u2 + U_bar[i][j].u2 - dt / dr * (Fr_bar[i][j].f2 - Fr_bar[i][j - 1].f2) - dt / dz * (Fz_bar[i][j].f2 - Fz_bar[i - 1][j].f2) + dt * s_bar[i][j].f2);
			//	U_bar2[i][j].u3 = 0.5 * (U[i][j].u3 + U_bar[i][j].u3 - dt / dr * (Fr_bar[i][j].f3 - Fr_bar[i][j - 1].f3) - dt / dz * (Fz_bar[i][j].f3 - Fz_bar[i - 1][j].f3) + dt * s_bar[i][j].f3);
			//	U_bar2[i][j].u4 = 0.5 * (U[i][j].u4 + U_bar[i][j].u4 - dt / dr * (Fr_bar[i][j].f4 - Fr_bar[i][j - 1].f4) - dt / dz * (Fz_bar[i][j].f4 - Fz_bar[i - 1][j].f4) + dt * s_bar[i][j].f4);
			//	U_bar2[i][j].u5 = 0.5 * (U[i][j].u5 + U_bar[i][j].u5 - dt / dr * (Fr_bar[i][j].f5 - Fr_bar[i][j - 1].f5) - dt / dz * (Fz_bar[i][j].f5 - Fz_bar[i - 1][j].f5) + dt * s_bar[i][j].f5);
			//	U_bar2[i][j].u6 = 0.5 * (U[i][j].u6 + U_bar[i][j].u6 - dt / dr * (Fr_bar[i][j].f6 - Fr_bar[i][j - 1].f6) - dt / dz * (Fz_bar[i][j].f6 - Fz_bar[i - 1][j].f6) + dt * s_bar[i][j].f6);
			//	U_bar2[i][j].u7 = 0.5 * (U[i][j].u7 + U_bar[i][j].u7 - dt / dr * (Fr_bar[i][j].f7 - Fr_bar[i][j - 1].f7) - dt / dz * (Fz_bar[i][j].f7 - Fz_bar[i - 1][j].f7) + dt * s_bar[i][j].f7);
			//	U_bar2[i][j].u8 = 0.5 * (U[i][j].u8 + U_bar[i][j].u8 - dt / dr * (Fr_bar[i][j].f8 - Fr_bar[i][j - 1].f8) - dt / dz * (Fz_bar[i][j].f8 - Fz_bar[i - 1][j].f8) + dt * s_bar[i][j].f8);
			//	U_bar2[i][j].u9 = 0.5 * (U[i][j].u9 + U_bar[i][j].u9 - dt / dr * (Fr_bar[i][j].f9 - Fr_bar[i][j - 1].f9) - dt / dz * (Fz_bar[i][j].f9 - Fz_bar[i - 1][j].f9) + dt * s_bar[i][j].f9);
			//	U_bar2[i][j].u10 = 0.5 * (U[i][j].u10 + U_bar[i][j].u10  - dt / dr * (Fr_bar[i][j].f10 - Fr_bar[i][j - 1].f10) - dt / dz * (Fz_bar[i][j].f10 - Fz_bar[i - 1][j].f10) + dt * s_bar[i][j].f10);
			//	U_bar2[i][j].u11 = 0.5 * (U[i][j].u11 + U_bar[i][j].u11  - dt / dr * (Fr_bar[i][j].f11 - Fr_bar[i][j - 1].f11) - dt / dz * (Fz_bar[i][j].f11 - Fz_bar[i - 1][j].f11) + dt * s_bar[i][j].f11);
			//	U_bar2[i][j].u12 = 0.5 * (U[i][j].u12 + U_bar[i][j].u12  - dt / dr * (Fr_bar[i][j].f12 - Fr_bar[i][j - 1].f12) - dt / dz * (Fz_bar[i][j].f12 - Fz_bar[i - 1][j].f12) + dt * s_bar[i][j].f12);
			//	U_bar2[i][j].u13 = 0.5 * (U[i][j].u13 + U_bar[i][j].u13  - dt / dr * (Fr_bar[i][j].f13 - Fr_bar[i][j - 1].f13) - dt / dz * (Fz_bar[i][j].f13 - Fz_bar[i - 1][j].f13) + dt * s_bar[i][j].f13);
			//}
			//else if(btype[i][j] == DOWN)
			//{


			//	U_bar2[i][j].u1 = 0.5 * (U[i][j].u1 + U_bar[i][j].u1 - dt / dr * (Fr_bar[i][j + 1].f1 - Fr_bar[i][j].f1) - dt / dz * (Fz_bar[i][j].f1 - Fz_bar[i - 1][j].f1) + dt * s_bar[i][j].f1);
			//	U_bar2[i][j].u2 = 0.5 * (U[i][j].u2 + U_bar[i][j].u2 - dt / dr * (Fr_bar[i][j + 1].f2 - Fr_bar[i][j].f2) - dt / dz * (Fz_bar[i][j].f2 - Fz_bar[i - 1][j].f2) + dt * s_bar[i][j].f2);
			//	U_bar2[i][j].u3 = 0.5 * (U[i][j].u3 + U_bar[i][j].u3 - dt / dr * (Fr_bar[i][j + 1].f3 - Fr_bar[i][j].f3) - dt / dz * (Fz_bar[i][j].f3 - Fz_bar[i - 1][j].f3) + dt * s_bar[i][j].f3);
			//	U_bar2[i][j].u4 = 0.5 * (U[i][j].u4 + U_bar[i][j].u4 - dt / dr * (Fr_bar[i][j + 1].f4 - Fr_bar[i][j].f4) - dt / dz * (Fz_bar[i][j].f4 - Fz_bar[i - 1][j].f4) + dt * s_bar[i][j].f4);
			//	U_bar2[i][j].u5 = 0.5 * (U[i][j].u5 + U_bar[i][j].u5 - dt / dr * (Fr_bar[i][j + 1].f5 - Fr_bar[i][j].f5) - dt / dz * (Fz_bar[i][j].f5 - Fz_bar[i - 1][j].f5) + dt * s_bar[i][j].f5);
			//	U_bar2[i][j].u6 = 0.5 * (U[i][j].u6 + U_bar[i][j].u6 - dt / dr * (Fr_bar[i][j + 1].f6 - Fr_bar[i][j].f6) - dt / dz * (Fz_bar[i][j].f6 - Fz_bar[i - 1][j].f6) + dt * s_bar[i][j].f6);
			//	U_bar2[i][j].u7 = 0.5 * (U[i][j].u7 + U_bar[i][j].u7 - dt / dr * (Fr_bar[i][j + 1].f7 - Fr_bar[i][j].f7) - dt / dz * (Fz_bar[i][j].f7 - Fz_bar[i - 1][j].f7) + dt * s_bar[i][j].f7);
			//	U_bar2[i][j].u8 = 0.5 * (U[i][j].u8 + U_bar[i][j].u8 - dt / dr * (Fr_bar[i][j + 1].f8 - Fr_bar[i][j].f8) - dt / dz * (Fz_bar[i][j].f8 - Fz_bar[i - 1][j].f8) + dt * s_bar[i][j].f8);
			//	U_bar2[i][j].u9 = 0.5 * (U[i][j].u9 + U_bar[i][j].u9 - dt / dr * (Fr_bar[i][j + 1].f9 - Fr_bar[i][j].f9) - dt / dz * (Fz_bar[i][j].f9 - Fz_bar[i - 1][j].f9) + dt * s_bar[i][j].f9);
			//	U_bar2[i][j].u10 = 0.5 * (U[i][j].u10 + U_bar[i][j].u10  - dt / dr * (Fr_bar[i][j + 1].f10 - Fr_bar[i][j].f10) - dt / dz * (Fz_bar[i][j].f10 - Fz_bar[i - 1][j].f10) + dt * s_bar[i][j].f10);
			//	U_bar2[i][j].u11 = 0.5 * (U[i][j].u11 + U_bar[i][j].u11  - dt / dr * (Fr_bar[i][j + 1].f11 - Fr_bar[i][j].f11) - dt / dz * (Fz_bar[i][j].f11 - Fz_bar[i - 1][j].f11) + dt * s_bar[i][j].f11);
			//	U_bar2[i][j].u12 = 0.5 * (U[i][j].u12 + U_bar[i][j].u12  - dt / dr * (Fr_bar[i][j + 1].f12 - Fr_bar[i][j].f12) - dt / dz * (Fz_bar[i][j].f12 - Fz_bar[i - 1][j].f12) + dt * s_bar[i][j].f12);
			//	U_bar2[i][j].u13 = 0.5 * (U[i][j].u13 + U_bar[i][j].u13  - dt / dr * (Fr_bar[i][j + 1].f13 - Fr_bar[i][j].f13) - dt / dz * (Fz_bar[i][j].f13 - Fz_bar[i - 1][j].f13) + dt * s_bar[i][j].f13);
			//}
			//else if(btype[i][j] == (LEFT + UP))
			//{


			//	U_bar2[i][j].u1 = 0.5 * (U[i][j].u1 + U_bar[i][j].u1 - dt / dr * (Fr_bar[i][j].f1 - Fr_bar[i][j - 1].f1) - dt / dz * (Fz_bar[i + 1][j].f1 - Fz_bar[i][j].f1) + dt * s_bar[i][j].f1);
			//	U_bar2[i][j].u2 = 0.5 * (U[i][j].u2 + U_bar[i][j].u2 - dt / dr * (Fr_bar[i][j].f2 - Fr_bar[i][j - 1].f2) - dt / dz * (Fz_bar[i + 1][j].f2 - Fz_bar[i][j].f2) + dt * s_bar[i][j].f2);
			//	U_bar2[i][j].u3 = 0.5 * (U[i][j].u3 + U_bar[i][j].u3 - dt / dr * (Fr_bar[i][j].f3 - Fr_bar[i][j - 1].f3) - dt / dz * (Fz_bar[i + 1][j].f3 - Fz_bar[i][j].f3) + dt * s_bar[i][j].f3);
			//	U_bar2[i][j].u4 = 0.5 * (U[i][j].u4 + U_bar[i][j].u4 - dt / dr * (Fr_bar[i][j].f4 - Fr_bar[i][j - 1].f4) - dt / dz * (Fz_bar[i + 1][j].f4 - Fz_bar[i][j].f4) + dt * s_bar[i][j].f4);
			//	U_bar2[i][j].u5 = 0.5 * (U[i][j].u5 + U_bar[i][j].u5 - dt / dr * (Fr_bar[i][j].f5 - Fr_bar[i][j - 1].f5) - dt / dz * (Fz_bar[i + 1][j].f5 - Fz_bar[i][j].f5) + dt * s_bar[i][j].f5);
			//	U_bar2[i][j].u6 = 0.5 * (U[i][j].u6 + U_bar[i][j].u6 - dt / dr * (Fr_bar[i][j].f6 - Fr_bar[i][j - 1].f6) - dt / dz * (Fz_bar[i + 1][j].f6 - Fz_bar[i][j].f6) + dt * s_bar[i][j].f6);
			//	U_bar2[i][j].u7 = 0.5 * (U[i][j].u7 + U_bar[i][j].u7 - dt / dr * (Fr_bar[i][j].f7 - Fr_bar[i][j - 1].f7) - dt / dz * (Fz_bar[i + 1][j].f7 - Fz_bar[i][j].f7) + dt * s_bar[i][j].f7);
			//	U_bar2[i][j].u8 = 0.5 * (U[i][j].u8 + U_bar[i][j].u8 - dt / dr * (Fr_bar[i][j].f8 - Fr_bar[i][j - 1].f8) - dt / dz * (Fz_bar[i + 1][j].f8 - Fz_bar[i][j].f8) + dt * s_bar[i][j].f8);
			//	U_bar2[i][j].u9 = 0.5 * (U[i][j].u9 + U_bar[i][j].u9 - dt / dr * (Fr_bar[i][j].f9 - Fr_bar[i][j - 1].f9) - dt / dz * (Fz_bar[i + 1][j].f9 - Fz_bar[i][j].f9) + dt * s_bar[i][j].f9);
			//	U_bar2[i][j].u10 = 0.5 * (U[i][j].u10 + U_bar[i][j].u10  - dt / dr * (Fr_bar[i][j].f10 - Fr_bar[i][j - 1].f10) - dt / dz * (Fz_bar[i + 1][j].f10 - Fz_bar[i][j].f10) + dt * s_bar[i][j].f10);
			//	U_bar2[i][j].u11 = 0.5 * (U[i][j].u11 + U_bar[i][j].u11  - dt / dr * (Fr_bar[i][j].f11 - Fr_bar[i][j - 1].f11) - dt / dz * (Fz_bar[i + 1][j].f11 - Fz_bar[i][j].f11) + dt * s_bar[i][j].f11);
			//	U_bar2[i][j].u12 = 0.5 * (U[i][j].u12 + U_bar[i][j].u12  - dt / dr * (Fr_bar[i][j].f12 - Fr_bar[i][j - 1].f12) - dt / dz * (Fz_bar[i + 1][j].f12 - Fz_bar[i][j].f12) + dt * s_bar[i][j].f12);
			//	U_bar2[i][j].u13 = 0.5 * (U[i][j].u13 + U_bar[i][j].u13  - dt / dr * (Fr_bar[i][j].f13 - Fr_bar[i][j - 1].f13) - dt / dz * (Fz_bar[i + 1][j].f13 - Fz_bar[i][j].f13) + dt * s_bar[i][j].f13);
			//}
			//else if(btype[i][j] == (LEFT + DOWN))
			//{


			//	U_bar2[i][j].u1 = 0.5 * (U[i][j].u1 + U_bar[i][j].u1 - dt / dr * (Fr_bar[i][j + 1].f1 - Fr_bar[i][j].f1) - dt / dz * (Fz_bar[i + 1][j].f1 - Fz_bar[i][j].f1) + dt * s_bar[i][j].f1);
			//	U_bar2[i][j].u2 = 0.5 * (U[i][j].u2 + U_bar[i][j].u2 - dt / dr * (Fr_bar[i][j + 1].f2 - Fr_bar[i][j].f2) - dt / dz * (Fz_bar[i + 1][j].f2 - Fz_bar[i][j].f2) + dt * s_bar[i][j].f2);
			//	U_bar2[i][j].u3 = 0.5 * (U[i][j].u3 + U_bar[i][j].u3 - dt / dr * (Fr_bar[i][j + 1].f3 - Fr_bar[i][j].f3) - dt / dz * (Fz_bar[i + 1][j].f3 - Fz_bar[i][j].f3) + dt * s_bar[i][j].f3);
			//	U_bar2[i][j].u4 = 0.5 * (U[i][j].u4 + U_bar[i][j].u4 - dt / dr * (Fr_bar[i][j + 1].f4 - Fr_bar[i][j].f4) - dt / dz * (Fz_bar[i + 1][j].f4 - Fz_bar[i][j].f4) + dt * s_bar[i][j].f4);
			//	U_bar2[i][j].u5 = 0.5 * (U[i][j].u5 + U_bar[i][j].u5 - dt / dr * (Fr_bar[i][j + 1].f5 - Fr_bar[i][j].f5) - dt / dz * (Fz_bar[i + 1][j].f5 - Fz_bar[i][j].f5) + dt * s_bar[i][j].f5);
			//	U_bar2[i][j].u6 = 0.5 * (U[i][j].u6 + U_bar[i][j].u6 - dt / dr * (Fr_bar[i][j + 1].f6 - Fr_bar[i][j].f6) - dt / dz * (Fz_bar[i + 1][j].f6 - Fz_bar[i][j].f6) + dt * s_bar[i][j].f6);
			//	U_bar2[i][j].u7 = 0.5 * (U[i][j].u7 + U_bar[i][j].u7 - dt / dr * (Fr_bar[i][j + 1].f7 - Fr_bar[i][j].f7) - dt / dz * (Fz_bar[i + 1][j].f7 - Fz_bar[i][j].f7) + dt * s_bar[i][j].f7);
			//	U_bar2[i][j].u8 = 0.5 * (U[i][j].u8 + U_bar[i][j].u8 - dt / dr * (Fr_bar[i][j + 1].f8 - Fr_bar[i][j].f8) - dt / dz * (Fz_bar[i + 1][j].f8 - Fz_bar[i][j].f8) + dt * s_bar[i][j].f8);
			//	U_bar2[i][j].u9 = 0.5 * (U[i][j].u9 + U_bar[i][j].u9 - dt / dr * (Fr_bar[i][j + 1].f9 - Fr_bar[i][j].f9) - dt / dz * (Fz_bar[i + 1][j].f9 - Fz_bar[i][j].f9) + dt * s_bar[i][j].f9);
			//	U_bar2[i][j].u10 = 0.5 * (U[i][j].u10 + U_bar[i][j].u10  - dt / dr * (Fr_bar[i][j + 1].f10 - Fr_bar[i][j].f10) - dt / dz * (Fz_bar[i + 1][j].f10 - Fz_bar[i][j].f10) + dt * s_bar[i][j].f10);
			//	U_bar2[i][j].u11 = 0.5 * (U[i][j].u11 + U_bar[i][j].u11  - dt / dr * (Fr_bar[i][j + 1].f11 - Fr_bar[i][j].f11) - dt / dz * (Fz_bar[i + 1][j].f11 - Fz_bar[i][j].f11) + dt * s_bar[i][j].f11);
			//	U_bar2[i][j].u12 = 0.5 * (U[i][j].u12 + U_bar[i][j].u12  - dt / dr * (Fr_bar[i][j + 1].f12 - Fr_bar[i][j].f12) - dt / dz * (Fz_bar[i + 1][j].f12 - Fz_bar[i][j].f12) + dt * s_bar[i][j].f12);
			//	U_bar2[i][j].u13 = 0.5 * (U[i][j].u13 + U_bar[i][j].u13  - dt / dr * (Fr_bar[i][j + 1].f13 - Fr_bar[i][j].f13) - dt / dz * (Fz_bar[i + 1][j].f13 - Fz_bar[i][j].f13) + dt * s_bar[i][j].f13);
			//}
			//else if(btype[i][j] == (RIGHT + DOWN))
			//{


			//	U_bar2[i][j].u1 = 0.5 * (U[i][j].u1 + U_bar[i][j].u1 - dt / dr * (Fr_bar[i][j + 1].f1 - Fr_bar[i][j].f1) - dt / dz * (Fz_bar[i][j].f1 - Fz_bar[i - 1][j].f1) + dt * s_bar[i][j].f1);
			//	U_bar2[i][j].u2 = 0.5 * (U[i][j].u2 + U_bar[i][j].u2 - dt / dr * (Fr_bar[i][j + 1].f2 - Fr_bar[i][j].f2) - dt / dz * (Fz_bar[i][j].f2 - Fz_bar[i - 1][j].f2) + dt * s_bar[i][j].f2);
			//	U_bar2[i][j].u3 = 0.5 * (U[i][j].u3 + U_bar[i][j].u3 - dt / dr * (Fr_bar[i][j + 1].f3 - Fr_bar[i][j].f3) - dt / dz * (Fz_bar[i][j].f3 - Fz_bar[i - 1][j].f3) + dt * s_bar[i][j].f3);
			//	U_bar2[i][j].u4 = 0.5 * (U[i][j].u4 + U_bar[i][j].u4 - dt / dr * (Fr_bar[i][j + 1].f4 - Fr_bar[i][j].f4) - dt / dz * (Fz_bar[i][j].f4 - Fz_bar[i - 1][j].f4) + dt * s_bar[i][j].f4);
			//	U_bar2[i][j].u5 = 0.5 * (U[i][j].u5 + U_bar[i][j].u5 - dt / dr * (Fr_bar[i][j + 1].f5 - Fr_bar[i][j].f5) - dt / dz * (Fz_bar[i][j].f5 - Fz_bar[i - 1][j].f5) + dt * s_bar[i][j].f5);
			//	U_bar2[i][j].u6 = 0.5 * (U[i][j].u6 + U_bar[i][j].u6 - dt / dr * (Fr_bar[i][j + 1].f6 - Fr_bar[i][j].f6) - dt / dz * (Fz_bar[i][j].f6 - Fz_bar[i - 1][j].f6) + dt * s_bar[i][j].f6);
			//	U_bar2[i][j].u7 = 0.5 * (U[i][j].u7 + U_bar[i][j].u7 - dt / dr * (Fr_bar[i][j + 1].f7 - Fr_bar[i][j].f7) - dt / dz * (Fz_bar[i][j].f7 - Fz_bar[i - 1][j].f7) + dt * s_bar[i][j].f7);
			//	U_bar2[i][j].u8 = 0.5 * (U[i][j].u8 + U_bar[i][j].u8 - dt / dr * (Fr_bar[i][j + 1].f8 - Fr_bar[i][j].f8) - dt / dz * (Fz_bar[i][j].f8 - Fz_bar[i - 1][j].f8) + dt * s_bar[i][j].f8);
			//	U_bar2[i][j].u9 = 0.5 * (U[i][j].u9 + U_bar[i][j].u9 - dt / dr * (Fr_bar[i][j + 1].f9 - Fr_bar[i][j].f9) - dt / dz * (Fz_bar[i][j].f9 - Fz_bar[i - 1][j].f9) + dt * s_bar[i][j].f9);
			//	U_bar2[i][j].u10 = 0.5 * (U[i][j].u10 + U_bar[i][j].u10  - dt / dr * (Fr_bar[i][j + 1].f10 - Fr_bar[i][j].f10) - dt / dz * (Fz_bar[i][j].f10 - Fz_bar[i - 1][j].f10) + dt * s_bar[i][j].f10);
			//	U_bar2[i][j].u11 = 0.5 * (U[i][j].u11 + U_bar[i][j].u11  - dt / dr * (Fr_bar[i][j + 1].f11 - Fr_bar[i][j].f11) - dt / dz * (Fz_bar[i][j].f11 - Fz_bar[i - 1][j].f11) + dt * s_bar[i][j].f11);
			//	U_bar2[i][j].u12 = 0.5 * (U[i][j].u12 + U_bar[i][j].u12  - dt / dr * (Fr_bar[i][j + 1].f12 - Fr_bar[i][j].f12) - dt / dz * (Fz_bar[i][j].f12 - Fz_bar[i - 1][j].f12) + dt * s_bar[i][j].f12);
			//	U_bar2[i][j].u13 = 0.5 * (U[i][j].u13 + U_bar[i][j].u13  - dt / dr * (Fr_bar[i][j + 1].f13 - Fr_bar[i][j].f13) - dt / dz * (Fz_bar[i][j].f13 - Fz_bar[i - 1][j].f13) + dt * s_bar[i][j].f13);
			//}
			//else if(btype[i][j] == (RIGHT + UP))
			//{

			//	U_bar2[i][j].u1 = 0.5 * (U[i][j].u1 + U_bar[i][j].u1 - dt / dr * (Fr_bar[i][j].f1 - Fr_bar[i][j - 1].f1) - dt / dz * (Fz_bar[i][j].f1 - Fz_bar[i - 1][j].f1) + dt * s_bar[i][j].f1);
			//	U_bar2[i][j].u2 = 0.5 * (U[i][j].u2 + U_bar[i][j].u2 - dt / dr * (Fr_bar[i][j].f2 - Fr_bar[i][j - 1].f2) - dt / dz * (Fz_bar[i][j].f2 - Fz_bar[i - 1][j].f2) + dt * s_bar[i][j].f2);
			//	U_bar2[i][j].u3 = 0.5 * (U[i][j].u3 + U_bar[i][j].u3 - dt / dr * (Fr_bar[i][j].f3 - Fr_bar[i][j - 1].f3) - dt / dz * (Fz_bar[i][j].f3 - Fz_bar[i - 1][j].f3) + dt * s_bar[i][j].f3);
			//	U_bar2[i][j].u4 = 0.5 * (U[i][j].u4 + U_bar[i][j].u4 - dt / dr * (Fr_bar[i][j].f4 - Fr_bar[i][j - 1].f4) - dt / dz * (Fz_bar[i][j].f4 - Fz_bar[i - 1][j].f4) + dt * s_bar[i][j].f4);
			//	U_bar2[i][j].u5 = 0.5 * (U[i][j].u5 + U_bar[i][j].u5 - dt / dr * (Fr_bar[i][j].f5 - Fr_bar[i][j - 1].f5) - dt / dz * (Fz_bar[i][j].f5 - Fz_bar[i - 1][j].f5) + dt * s_bar[i][j].f5);
			//	U_bar2[i][j].u6 = 0.5 * (U[i][j].u6 + U_bar[i][j].u6 - dt / dr * (Fr_bar[i][j].f6 - Fr_bar[i][j - 1].f6) - dt / dz * (Fz_bar[i][j].f6 - Fz_bar[i - 1][j].f6) + dt * s_bar[i][j].f6);
			//	U_bar2[i][j].u7 = 0.5 * (U[i][j].u7 + U_bar[i][j].u7 - dt / dr * (Fr_bar[i][j].f7 - Fr_bar[i][j - 1].f7) - dt / dz * (Fz_bar[i][j].f7 - Fz_bar[i - 1][j].f7) + dt * s_bar[i][j].f7);
			//	U_bar2[i][j].u8 = 0.5 * (U[i][j].u8 + U_bar[i][j].u8 - dt / dr * (Fr_bar[i][j].f8 - Fr_bar[i][j - 1].f8) - dt / dz * (Fz_bar[i][j].f8 - Fz_bar[i - 1][j].f8) + dt * s_bar[i][j].f8);
			//	U_bar2[i][j].u9 = 0.5 * (U[i][j].u9 + U_bar[i][j].u9 - dt / dr * (Fr_bar[i][j].f9 - Fr_bar[i][j - 1].f9) - dt / dz * (Fz_bar[i][j].f9 - Fz_bar[i - 1][j].f9) + dt * s_bar[i][j].f9);
			//	U_bar2[i][j].u10 = 0.5 * (U[i][j].u10 + U_bar[i][j].u10  - dt / dr * (Fr_bar[i][j].f10 - Fr_bar[i][j - 1].f10) - dt / dz * (Fz_bar[i][j].f10 - Fz_bar[i - 1][j].f10) + dt * s_bar[i][j].f10);
			//	U_bar2[i][j].u11 = 0.5 * (U[i][j].u11 + U_bar[i][j].u11  - dt / dr * (Fr_bar[i][j].f11 - Fr_bar[i][j - 1].f11) - dt / dz * (Fz_bar[i][j].f11 - Fz_bar[i - 1][j].f11) + dt * s_bar[i][j].f11);
			//	U_bar2[i][j].u12 = 0.5 * (U[i][j].u12 + U_bar[i][j].u12  - dt / dr * (Fr_bar[i][j].f12 - Fr_bar[i][j - 1].f12) - dt / dz * (Fz_bar[i][j].f12 - Fz_bar[i - 1][j].f12) + dt * s_bar[i][j].f12);
			//	U_bar2[i][j].u13 = 0.5 * (U[i][j].u13 + U_bar[i][j].u13  - dt / dr * (Fr_bar[i][j].f13 - Fr_bar[i][j - 1].f13) - dt / dz * (Fz_bar[i][j].f13 - Fz_bar[i - 1][j].f13) + dt * s_bar[i][j].f13);
			//}
			else //if(btype[i][j] == 0)
			{
				U_bar2[i][j].u1 = 0;
				U_bar2[i][j].u2 = 0;
				U_bar2[i][j].u3 = 0;
				U_bar2[i][j].u4 = 0;
				U_bar2[i][j].u5 = 0;
				U_bar2[i][j].u6 = 0;
				U_bar2[i][j].u7 = 0;
				U_bar2[i][j].u8 = 0;
				U_bar2[i][j].u9 = 0;
				U_bar2[i][j].u10 = 0;
				U_bar2[i][j].u11 = 0;
				U_bar2[i][j].u12 = 0;
				U_bar2[i][j].u13 = 0;

			}



			//Fr_bar2[i][j] = cal_fr(U_bar2[i][j]);
			//Fz_bar2[i][j] = cal_fz(U_bar2[i][j]);
		}
	}

	///计算
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nz; j++)
		{		
			if (btype[i][j] == 1)
			{

				double Qr = 0;
				double Qz = 0;
				if (abs(MPDT[i + 1][j].ne + 2 * MPDT[i][j].ne + MPDT[i - 1][j].ne) != 0)
				{
					Qz = abs(MPDT[i + 1][j].ne - 2 * MPDT[i][j].ne + MPDT[i - 1][j].ne) / abs(MPDT[i + 1][j].ne + 2 * MPDT[i][j].ne + MPDT[i - 1][j].ne);
				}
				if (abs(MPDT[i][j + 1].ne + 2 * MPDT[i][j].ne + MPDT[i][j - 1].ne) != 0)
				{
					Qr = abs(MPDT[i][j + 1].ne - 2 * MPDT[i][j].ne + MPDT[i][j - 1].ne) / abs(MPDT[i][j + 1].ne + 2 * MPDT[i][j].ne + MPDT[i][j - 1].ne);
				}

				if (i == row && j == col)
				{
					printf("Qr = %e\n", Qr);
					printf("Qz = %e\n", Qz);
				}

				//if (Qz > 1)
				{
					Qz = 0;
					//printf("Qz = %e\n", Qz);
				}
				//if (Qr > 1)
				{
					Qr = 0;
					//printf("Qr = %e\n", Qr);
				}

				U[i][j].u1 = U_bar2[i][j].u1 + Qr / 2 * (U_bar2[i][j + 1].u1 - 2 * U_bar2[i][j].u1 + U_bar2[i][j - 1].u1) + Qz / 2 * (U_bar2[i + 1][j].u1 - 2 * U_bar2[i][j].u1 + U_bar2[i - 1][j].u1);
				U[i][j].u2 = U_bar2[i][j].u2 + Qr / 2 * (U_bar2[i][j + 1].u2 - 2 * U_bar2[i][j].u2 + U_bar2[i][j - 1].u2) + Qz / 2 * (U_bar2[i + 1][j].u2 - 2 * U_bar2[i][j].u2 + U_bar2[i - 1][j].u2);
				U[i][j].u3 = U_bar2[i][j].u3 + Qr / 2 * (U_bar2[i][j + 1].u3 - 2 * U_bar2[i][j].u3 + U_bar2[i][j - 1].u3) + Qz / 2 * (U_bar2[i + 1][j].u3 - 2 * U_bar2[i][j].u3 + U_bar2[i - 1][j].u3);
				U[i][j].u4 = U_bar2[i][j].u4 + Qr / 2 * (U_bar2[i][j + 1].u4 - 2 * U_bar2[i][j].u4 + U_bar2[i][j - 1].u4) + Qz / 2 * (U_bar2[i + 1][j].u4 - 2 * U_bar2[i][j].u4 + U_bar2[i - 1][j].u4);
				U[i][j].u5 = U_bar2[i][j].u5 + Qr / 2 * (U_bar2[i][j + 1].u5 - 2 * U_bar2[i][j].u5 + U_bar2[i][j - 1].u5) + Qz / 2 * (U_bar2[i + 1][j].u5 - 2 * U_bar2[i][j].u5 + U_bar2[i - 1][j].u5);
				U[i][j].u6 = U_bar2[i][j].u6 + Qr / 2 * (U_bar2[i][j + 1].u6 - 2 * U_bar2[i][j].u6 + U_bar2[i][j - 1].u6) + Qz / 2 * (U_bar2[i + 1][j].u6 - 2 * U_bar2[i][j].u6 + U_bar2[i - 1][j].u6);
				U[i][j].u7 = U_bar2[i][j].u7 + Qr / 2 * (U_bar2[i][j + 1].u7 - 2 * U_bar2[i][j].u7 + U_bar2[i][j - 1].u7) + Qz / 2 * (U_bar2[i + 1][j].u7 - 2 * U_bar2[i][j].u7 + U_bar2[i - 1][j].u7);
				U[i][j].u8 = U_bar2[i][j].u8 + Qr / 2 * (U_bar2[i][j + 1].u8 - 2 * U_bar2[i][j].u8 + U_bar2[i][j - 1].u8) + Qz / 2 * (U_bar2[i + 1][j].u8 - 2 * U_bar2[i][j].u8 + U_bar2[i - 1][j].u8);
				U[i][j].u9 = U_bar2[i][j].u9 + Qr / 2 * (U_bar2[i][j + 1].u9 - 2 * U_bar2[i][j].u9 + U_bar2[i][j - 1].u9) + Qz / 2 * (U_bar2[i + 1][j].u9 - 2 * U_bar2[i][j].u9 + U_bar2[i - 1][j].u9);
				U[i][j].u10 = U_bar2[i][j].u10 + Qr / 2 * (U_bar2[i][j + 1].u10 - 2 * U_bar2[i][j].u10 + U_bar2[i][j - 1].u10) + Qz / 2 * (U_bar2[i + 1][j].u10 - 2 * U_bar2[i][j].u10 + U_bar2[i - 1][j].u10);
				U[i][j].u11 = U_bar2[i][j].u11 + Qr / 2 * (U_bar2[i][j + 1].u11 - 2 * U_bar2[i][j].u11 + U_bar2[i][j - 1].u11) + Qz / 2 * (U_bar2[i + 1][j].u11 - 2 * U_bar2[i][j].u11 + U_bar2[i - 1][j].u11);
				U[i][j].u12 = U_bar2[i][j].u12 + Qr / 2 * (U_bar2[i][j + 1].u12 - 2 * U_bar2[i][j].u12 + U_bar2[i][j - 1].u12) + Qz / 2 * (U_bar2[i + 1][j].u12 - 2 * U_bar2[i][j].u12 + U_bar2[i - 1][j].u12);
				U[i][j].u13 = U_bar2[i][j].u13 + Qr / 2 * (U_bar2[i][j + 1].u13 - 2 * U_bar2[i][j].u13 + U_bar2[i][j - 1].u13) + Qz / 2 * (U_bar2[i + 1][j].u13 - 2 * U_bar2[i][j].u13 + U_bar2[i - 1][j].u13);
				

				if (i == row && j == col)
				{
					//cout << "final U[i][j].u1 = " << U[i][j].u1 << endl;
					printf("final U[i][j].u1 = %.5e\n", U[i][j].u1);
					//printf("Fr_bar[i][j].f1 = %.5e\n", Fr_bar[i][j].f1);
					//printf("Fr_bar[i][j - 1].f1 = %.5e\n", Fr_bar[i][j - 1].f1);
					//printf("Fz_bar[i][j].f1 = %.5e\n", Fz_bar[i][j].f1);
					//printf("Fz_bar[i - 1][j].f1 = %.5e\n", Fz_bar[i - 1][j].f1);
					//printf("s_bar[i][j].f1 = %.5e\n", s_bar[i][j].f1);
					//printf("U_bar2[i][j].u1 = %.5e\n", U_bar2[i][j].u1);

				}
				//if (U[i][j].u5 != U[i][j].u5)
				//{
				//	printf("U_bar = %.5e\n", U_bar[i][j].u5);
				//	printf("Fz = %.5e\n", Fz[i][j].f5);
				//	printf("Fr = %.5e\n", Fr[i][j].f5);
				//	printf("s = %.5e\n", s[i][j].f5);
				//	printf("U_bar2 = %.5e\n", U_bar2[i][j].u5);
				//	printf("Fz_bar = %.5e\n", Fz_bar[i][j].f5);
				//	printf("Fr_bar = %.5e\n", Fr_bar[i][j].f5);
				//	printf("s_bar = %.5e\n", s_bar[i][j].f5);

				//}

			}
			else if (btype[i][j] == LEFT)//左边的边界，复制右边的参数
			{
				U[i][j].u1 = U[i + 1][j].u1;
				U[i][j].u2 = U[i + 1][j].u2;
				U[i][j].u3 = U[i + 1][j].u3;
				U[i][j].u4 = U[i + 1][j].u4;
				U[i][j].u5 = U[i + 1][j].u5;
				U[i][j].u6 = U[i + 1][j].u6;
				U[i][j].u7 = U[i + 1][j].u7;
				U[i][j].u8 = U[i + 1][j].u8;
				U[i][j].u9 = U[i + 1][j].u9;
				U[i][j].u10 = U[i + 1][j].u10;
				U[i][j].u11 = U[i + 1][j].u11;
				U[i][j].u12 = U[i + 1][j].u12;
				U[i][j].u13 = U[i + 1][j].u13;
			}
			else if (btype[i][j] == RIGHT)
			{
				U[i][j].u1 = U[i - 1][j].u1;
				U[i][j].u2 = U[i - 1][j].u2;
				U[i][j].u3 = U[i - 1][j].u3;
				U[i][j].u4 = U[i - 1][j].u4;
				U[i][j].u5 = U[i - 1][j].u5;
				U[i][j].u6 = U[i - 1][j].u6;
				U[i][j].u7 = U[i - 1][j].u7;
				U[i][j].u8 = U[i - 1][j].u8;
				U[i][j].u9 = U[i - 1][j].u9;
				U[i][j].u10 = U[i - 1][j].u10;
				U[i][j].u11 = U[i - 1][j].u11;
				U[i][j].u12 = U[i - 1][j].u12;
				U[i][j].u13 = U[i - 1][j].u13;
			}
			else if (btype[i][j] == UP)
			{
				U[i][j].u1 = U[i][j - 1].u1;
				U[i][j].u2 = U[i][j - 1].u2;
				U[i][j].u3 = U[i][j - 1].u3;
				U[i][j].u4 = U[i][j - 1].u4;
				U[i][j].u5 = U[i][j - 1].u5;
				U[i][j].u6 = U[i][j - 1].u6;
				U[i][j].u7 = U[i][j - 1].u7;
				U[i][j].u8 = U[i][j - 1].u8;
				U[i][j].u9 = U[i][j - 1].u9;
				U[i][j].u10 = U[i][j - 1].u10;
				U[i][j].u11 = U[i][j - 1].u11;
				U[i][j].u12 = U[i][j - 1].u12;
				U[i][j].u13 = U[i][j - 1].u13;
			}
			else if (btype[i][j] == DOWN)
			{
				U[i][j].u1 = U[i][j + 1].u1;
				U[i][j].u2 = U[i][j + 1].u2;
				U[i][j].u3 = U[i][j + 1].u3;
				U[i][j].u4 = U[i][j + 1].u4;
				U[i][j].u5 = U[i][j + 1].u5;
				U[i][j].u6 = U[i][j + 1].u6;
				U[i][j].u7 = U[i][j + 1].u7;
				U[i][j].u8 = U[i][j + 1].u8;
				U[i][j].u9 = U[i][j + 1].u9;
				U[i][j].u10 = U[i][j + 1].u10;
				U[i][j].u11 = U[i][j + 1].u11;
				U[i][j].u12 = U[i][j + 1].u12;
				U[i][j].u13 = U[i][j + 1].u13;
			}
			else if (btype[i][j] == (LEFT + UP))
			{
				U[i][j].u1 = U[i + 1][j - 1].u1;
				U[i][j].u2 = U[i + 1][j - 1].u2;
				U[i][j].u3 = U[i + 1][j - 1].u3;
				U[i][j].u4 = U[i + 1][j - 1].u4;
				U[i][j].u5 = U[i + 1][j - 1].u5;
				U[i][j].u6 = U[i + 1][j - 1].u6;
				U[i][j].u7 = U[i + 1][j - 1].u7;
				U[i][j].u8 = U[i + 1][j - 1].u8;
				U[i][j].u9 = U[i + 1][j - 1].u9;
				U[i][j].u10 = U[i + 1][j - 1].u10;
				U[i][j].u11 = U[i + 1][j - 1].u11;
				U[i][j].u12 = U[i + 1][j - 1].u12;
				U[i][j].u13 = U[i + 1][j - 1].u13;
			}
			else if (btype[i][j] == (LEFT + DOWN))
			{
				U[i][j].u1 = U[i + 1][j + 1].u1;
				U[i][j].u2 = U[i + 1][j + 1].u2;
				U[i][j].u3 = U[i + 1][j + 1].u3;
				U[i][j].u4 = U[i + 1][j + 1].u4;
				U[i][j].u5 = U[i + 1][j + 1].u5;
				U[i][j].u6 = U[i + 1][j + 1].u6;
				U[i][j].u7 = U[i + 1][j + 1].u7;
				U[i][j].u8 = U[i + 1][j + 1].u8;
				U[i][j].u9 = U[i + 1][j + 1].u9;
				U[i][j].u10 = U[i + 1][j + 1].u10;
				U[i][j].u11 = U[i + 1][j + 1].u11;
				U[i][j].u12 = U[i + 1][j + 1].u12;
				U[i][j].u13 = U[i + 1][j + 1].u13;
			}
			else if (btype[i][j] == (RIGHT + DOWN))
			{
				U[i][j].u1 = U[i - 1][j + 1].u1;
				U[i][j].u2 = U[i - 1][j + 1].u2;
				U[i][j].u3 = U[i - 1][j + 1].u3;
				U[i][j].u4 = U[i - 1][j + 1].u4;
				U[i][j].u5 = U[i - 1][j + 1].u5;
				U[i][j].u6 = U[i - 1][j + 1].u6;
				U[i][j].u7 = U[i - 1][j + 1].u7;
				U[i][j].u8 = U[i - 1][j + 1].u8;
				U[i][j].u9 = U[i - 1][j + 1].u9;
				U[i][j].u10 = U[i - 1][j + 1].u10;
				U[i][j].u11 = U[i - 1][j + 1].u11;
				U[i][j].u12 = U[i - 1][j + 1].u12;
				U[i][j].u13 = U[i - 1][j + 1].u13;
			}
			else if (btype[i][j] == (RIGHT + UP))
			{
				U[i][j].u1 = U[i - 1][j - 1].u1;
				U[i][j].u2 = U[i - 1][j - 1].u2;
				U[i][j].u3 = U[i - 1][j - 1].u3;
				U[i][j].u4 = U[i - 1][j - 1].u4;
				U[i][j].u5 = U[i - 1][j - 1].u5;
				U[i][j].u6 = U[i - 1][j - 1].u6;
				U[i][j].u7 = U[i - 1][j - 1].u7;
				U[i][j].u8 = U[i - 1][j - 1].u8;
				U[i][j].u9 = U[i - 1][j - 1].u9;
				U[i][j].u10 = U[i - 1][j - 1].u10;
				U[i][j].u11 = U[i - 1][j - 1].u11;
				U[i][j].u12 = U[i - 1][j - 1].u12;
				U[i][j].u13 = U[i - 1][j - 1].u13;
			}
			else if (btype[i][j] == 0)
			{
				U[i][j].u1 = 0;
				U[i][j].u2 = 0;
				U[i][j].u3 = 0;
				U[i][j].u4 = 0;
				U[i][j].u5 = 0;
				U[i][j].u6 = 0;
				U[i][j].u7 = 0;
				U[i][j].u8 = 0;
				U[i][j].u9 = 0;
				U[i][j].u10 = 0;
				U[i][j].u11 = 0;
				U[i][j].u12 = 0;
				U[i][j].u13 = 0;

			}
			else 
			{
				U[i][j].u1 = U_bar2[i][j].u1;
				U[i][j].u2 = U_bar2[i][j].u2;
				U[i][j].u3 = U_bar2[i][j].u3;
				U[i][j].u4 = U_bar2[i][j].u4;
				U[i][j].u5 = U_bar2[i][j].u5;
				U[i][j].u6 = U_bar2[i][j].u6;
				U[i][j].u7 = U_bar2[i][j].u7;
				U[i][j].u8 = U_bar2[i][j].u8;
				U[i][j].u9 = U_bar2[i][j].u9;
				U[i][j].u10 = U_bar2[i][j].u10;
				U[i][j].u11 = U_bar2[i][j].u11;
				U[i][j].u12 = U_bar2[i][j].u12;
				U[i][j].u13 = U_bar2[i][j].u13;
			}

		}
	}
	//output_u_all();

	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nz; j++)
		{
			MPDT[i][j].ne = U[i][j].u1;
			MPDT[i][j].ni = U[i][j].u2;

			MPDT[i][j].ver = U[i][j].u3;
			MPDT[i][j].vetheta = U[i][j].u4;
			MPDT[i][j].vez = U[i][j].u5;


			MPDT[i][j].vir = U[i][j].u6;
			MPDT[i][j].vitheta = U[i][j].u7 ;
			MPDT[i][j].viz = U[i][j].u8 ;


			MPDT[i][j].br = U[i][j].u9;
			MPDT[i][j].btheta = U[i][j].u10;
			MPDT[i][j].bz = U[i][j].u11;

			MPDT[i][j].pe = U[i][j].u12;
			MPDT[i][j].pi = U[i][j].u13;

			MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1)  + 0.5 * MPDT[i][j].ne * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
			MPDT[i][j].ei = MPDT[i][j].pi / (gamma - 1)  + 0.5 * MPDT[i][j].ni * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);


		}
	}


	return;
}

void ion_flow()
{
	return;
}


struct _U cal_u(int i, int j)
{
	struct _U uij;
	uij.r = i;
	uij.z = j;
	uij.u1 = MPDT[i][j].ne;
	uij.u2 = MPDT[i][j].ni;
	uij.u3 = MPDT[i][j].ver;
	uij.u4 = MPDT[i][j].vetheta;
	uij.u5 = MPDT[i][j].vez;
	uij.u6 = MPDT[i][j].vir;
	uij.u7 = MPDT[i][j].vitheta;
	uij.u8 = MPDT[i][j].viz;
	uij.u9 = MPDT[i][j].br;
	uij.u10 = MPDT[i][j].btheta;
	uij.u11 = MPDT[i][j].bz;
	uij.u12 = MPDT[i][j].pe;
	uij.u13 = MPDT[i][j].pi;
	return uij;
}

struct _F cal_fr(struct _U uij)
{
	struct _F fij;

	double ne = uij.u1;
	double ni = uij.u2;
	double ver = 0;
	double vetheta = 0;
	double vez = 0;
	double vir = 0;
	double vitheta = 0;
	double viz = 0;

	ver = uij.u3;
	vetheta = uij.u4;
	vez = uij.u5;
	vir = uij.u6;
	vitheta = uij.u7;
	viz = uij.u8;

	double br = uij.u9;
	double btheta = uij.u10;
	double bz = uij.u11;

	double pe = uij.u12;
	double pi = uij.u13;

	//double ee = pe / (gamma - 1) + 0.5 * ne * ME * (ver * ver + vetheta * vetheta + vez * vez);
	//double ei = pi / (gamma - 1) + 0.5 * ni * MI * (vir * vir + vitheta * vitheta + viz * viz);

	fij.r = uij.r;
	fij.z = uij.z;
	fij.f1 = ne * ver;
	fij.f2 = ni * vir;
	if (ne == 0)
	{
		fij.f3 = 0;
		fij.f4 = 0;
		fij.f5 = 0;
	}
	else
	{
		fij.f3 = (pe + (btheta * btheta + bz * bz - br * br) / (2 * MU_0)) / (ne * ME);
		fij.f4 = -(btheta * br) / (MU_0 * ne * ME);
		fij.f5 = -(bz * br) / (MU_0 * ne * ME);
	}

	if (ni == 0)
	{
		fij.f6 = 0;
		fij.f7 = 0;
		fij.f8 = 0;
	}
	else
	{
		fij.f6 = (pi + (btheta * btheta + bz * bz - br * br) / (2 * MU_0)) / (ni * MI);
		fij.f7 = -(btheta * br) / (MU_0 * ni * MI);
		fij.f8 = -(bz * br) / (MU_0 * ni * MI);
	}

	fij.f9 = 0;
	fij.f10 = 0.5 * (vir * btheta - ver * btheta - br * vitheta + br * vetheta);
	fij.f11 = 0.5 * (vir * bz - ver * bz - br * viz + br * vez);
	fij.f12 = pe * ver + (gamma - 1) * ( -ver * (btheta * btheta + bz * bz - br * br) / (2 * MU_0) + vetheta * btheta * br / MU_0 + vez * bz * br / MU_0);
	fij.f13 = pi * vir + (gamma - 1) * ( -vir * (btheta * btheta + bz * bz - br * br) / (2 * MU_0) + vitheta * btheta * br / MU_0 + viz * bz * br / MU_0);
	return fij;
}

struct _F cal_fz(struct _U uij)
{
	struct _F fij;

	double ne = uij.u1;
	double ni = uij.u2;
	double ver = 0;
	double vetheta = 0;
	double vez = 0;
	double vir = 0;
	double vitheta = 0;
	double viz = 0;

	ver = uij.u3;
	vetheta = uij.u4;
	vez = uij.u5;
	vir = uij.u6;
	vitheta = uij.u7;
	viz = uij.u8;

	double br = uij.u9;
	double btheta = uij.u10;
	double bz = uij.u11;

	double pe = uij.u12;
	double pi = uij.u13;

	//double ee = pe / (gamma - 1) + 0.5 * ne * ME * (ver * ver + vetheta * vetheta + vez * vez);
	//double ei = pi / (gamma - 1) + 0.5 * ni * MI * (vir * vir + vitheta * vitheta + viz * viz);


	fij.r = uij.r;
	fij.z = uij.z;

	fij.f1 = ne * vez;
	fij.f2 = ni * viz;
	if (ne == 0)
	{
		fij.f3 = 0;
		fij.f4 = 0;
		fij.f5 = 0;
	}
	else
	{
		fij.f3 = - (bz * br) / (MU_0 * ne * ME);
		fij.f4 = - (btheta * bz) / (MU_0 * ne * ME);
		fij.f5 = (pe + (btheta * btheta + br * br - bz * bz) / (2 * MU_0)) / (ne * ME);
	}

	if (ni == 0)
	{
		fij.f6 = 0;
		fij.f7 = 0;
		fij.f8 = 0;
	}
	else
	{
		fij.f6 = -(bz * br) / (MU_0 * ni * MI);
		fij.f7 = -(btheta * bz) / (MU_0 * ni * MI);
		fij.f8 = (pi + (btheta * btheta + br * br - bz * bz) / (2 * MU_0)) / (ni * MI);
	}

	fij.f9 = 0.5 * (viz * br - vez * br + bz * ver - bz * vir);
	fij.f10 = 0.5 * (viz * btheta - vez * btheta - bz * vitheta + bz * vetheta);
	fij.f11 = 0;
	fij.f12 = pe * vez + (gamma - 1) * (-vez * (btheta * btheta + br * br - bz * bz) / (2 * MU_0) + vetheta * btheta * bz / MU_0 + ver * bz * br / MU_0);
	fij.f13 = pi * viz + (gamma - 1) * (-viz * (btheta * btheta + br * br - bz * bz) / (2 * MU_0) + vitheta * btheta * bz / MU_0 + vir * bz * br / MU_0);
	return fij;
}


struct _F cal_s(struct _U uij)
{
	struct _F fij;

	double ne = uij.u1;
	double ni = uij.u2;
	double ver = 0;
	double vetheta = 0;
	double vez = 0;
	double vir = 0;
	double vitheta = 0;
	double viz = 0;

	ver = uij.u3;
	vetheta = uij.u4;
	vez = uij.u5;
	vir = uij.u6;
	vitheta = uij.u7;
	viz = uij.u8;

	double br = uij.u9;
	double btheta = uij.u10;
	double bz = uij.u11;

	double pe = uij.u12;
	double pi = uij.u13;

	double ee = pe / (gamma - 1) + 0.5 * ne * ME * (ver * ver + vetheta * vetheta + vez * vez);
	double ei = pi / (gamma - 1) + 0.5 * ni * MI * (vir * vir + vitheta * vitheta + viz * viz);

	fij.r = uij.r;
	fij.z = uij.z;
	double r1 = 0;
	if (uij.r == 0)
	{
		r1 = -2 / dr;
	}
	else 
	{
		r1 = -1 / (uij.r * dr);
	}
	fij.f1 = r1 * ne * ver;
	fij.f2 = r1 * ni * vir;

	if (ne == 0)
	{
		fij.f3 = 0;
		fij.f4 = 0;
		fij.f5 = 0;
	}
	else
	{
		fij.f3 = r1 * (pe + (btheta * btheta + bz * bz - br * br) / (2 * MU_0)) / (ne * ME);
		fij.f4 = -r1 * (btheta * br) / (MU_0 * ne * ME);
		fij.f5 = -r1 * (bz * br) / (MU_0 * ne * ME);
	}

	if (ni == 0)
	{
		fij.f6 = 0;
		fij.f7 = 0;
		fij.f8 = 0;
	}
	else
	{
		fij.f6 = r1 * (pi + (btheta * btheta + bz * bz - br * br) / (2 * MU_0)) / (ni * MI);
		fij.f7 = -r1 * (btheta * br) / (MU_0 * ni * MI);
		fij.f8 = -r1 * (bz * br) / (MU_0 * ni * MI);
	}

	fij.f9 = 0;
	fij.f10 = r1 * 0.5 * (vir * btheta - ver * btheta - br * vitheta + br * vetheta);
	fij.f11 = r1 * 0.5 * (vir * bz - ver * bz - br * viz + br * vez);
	fij.f12 = r1 * (pe * ver + (gamma - 1) * (-ver * (btheta * btheta + bz * bz - br * br) / (2 * MU_0) + vetheta * btheta * br / MU_0 + vez * bz * br / MU_0));
	fij.f13 = r1 * (pi * vir + (gamma - 1) * (-vir * (btheta * btheta + bz * bz - br * br) / (2 * MU_0) + vitheta * btheta * br / MU_0 + viz * bz * br / MU_0));
	return fij;
}