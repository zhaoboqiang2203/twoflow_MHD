#include "MHD.h"

/*************************************************

Copyright (c) ,2015-2016,  
File name:
Author:   

Date:2010-08-25

Description:描述主要实现的功能
Fucntion list:...

Histiory : //修改历史记录列表，每条修改记录修改日期，修改者，以及修改内容

1.Data:....
  Author:...
  Midification:...

2.


**************************************************/



/*************************************************

Function:                  // 函数名称

Description:               // hanshu ogngnegn 

Calls:                     // 无

Input:                     // 输入参数说明，包括每个参数的作

						   // 用、取值说明及参数间关系。

Output:                    // 对输出参数的说明。

Return:                    // 函数返回值的说明

Others:                   // 其它说明

*************************************************/

void move()
{
	double q_half;
	double hrho, htheta, hz, h2;
	double srho, stheta, sz;
	double urho, utheta, uz;
	double urho_half, utheta_half, uz_half;
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nz; j++)
		{

			q_half = -dt * QE / ME / 2;
			hrho = q_half * app_Br[i][j];
			htheta = 0;
			hz = q_half * app_Bz[i][j];
			h2 = sqr(hrho) + sqr(hz);
			
			srho = 2 * hrho / (1 + h2);
			stheta = 2 * htheta / (1 + h2); 
			sz = 2 * hz / (1 + h2);
			urho = MPDT[i][j].ver + q_half * Er[i][j];
			utheta = MPDT[i][j].vetheta;
			uz = MPDT[i][j].vez + q_half * Ez[i][j];
			urho_half = urho + (utheta * sz - uz * stheta) + ((uz * hrho - urho * hz) * sz - (urho * htheta - utheta * hrho) * stheta);
			utheta_half = utheta + (uz * srho - urho * sz) + ((urho * htheta - utheta * hrho) * srho - (utheta * hz - uz * htheta) * sz);
			uz_half = uz + (urho * stheta - utheta * srho) + ((utheta * hz - uz * htheta) * stheta - (uz * hrho - urho * hz) * srho);
			MPDT[i][j].ver = urho_half + q_half * Er[i][j];
			MPDT[i][j].vetheta = utheta_half;
			MPDT[i][j].vez = uz_half + q_half * Ez[i][j];
		
			//MPDT[i][j].ver = -QE * Er[i][j] / ME * dt;
			//MPDT[i][j].vez = -QE * Ez[i][j] / ME * dt;

			MPDT[i][j].ee = 0.5 * MPDT[i][j].ne * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
			//MPDT[i][j].pe = MPDT[i][j].ne * MPDT[i][j].ee;


			q_half = dt * QE / MI / 2;
			hrho = q_half * app_Br[i][j];
			htheta = 0;
			hz = q_half * app_Bz[i][j];
			h2 = sqr(hrho) + sqr(hz);

			srho = 2 * hrho / (1 + h2);
			stheta = 2 * htheta / (1 + h2);
			sz = 2 * hz / (1 + h2);
			urho = MPDT[i][j].vir + q_half * Er[i][j];
			utheta = MPDT[i][j].vitheta;
			uz = MPDT[i][j].viz + q_half * Ez[i][j];
			urho_half = urho + (utheta * sz - uz * stheta) + ((uz * hrho - urho * hz) * sz - (urho * htheta - utheta * hrho) * stheta);
			utheta_half = utheta + (uz * srho - urho * sz) + ((urho * htheta - utheta * hrho) * srho - (utheta * hz - uz * htheta) * sz);
			uz_half = uz + (urho * stheta - utheta * srho) + ((utheta * hz - uz * htheta) * stheta - (uz * hrho - urho * hz) * srho);
			MPDT[i][j].vir = urho_half + q_half * Er[i][j];
			MPDT[i][j].vitheta = utheta_half;
			MPDT[i][j].viz = uz_half + q_half * Ez[i][j];

			//MPDT[i][j].vir = QE * Er[i][j] / MI * dt;
			//MPDT[i][j].viz = QE * Ez[i][j] / MI * dt;

			MPDT[i][j].ei = 0.5 * MPDT[i][j].ni * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);
			//MPDT[i][j].pi = MPDT[i][j].ni * MPDT[i][j].ei;


			double Q_ei = 1.4e-20;//电子离子碰撞截面
			MPDT[i][j].mu_ie = MPDT[i][j].ni * Q_ei * ((8 * MPDT[i][j].ee) / (PI * ME));

			double m_eir = 0;        //电子向离子转移的径向（r方向）动量
			double m_eitheta = 0;    //电子向离子转移的角向（theta方向）动量
			double m_eiz = 0;        //电子向离子转移的轴向（z方向）动量

			double delta_ei = 0;     //电子向离子转移的能量
			
			double eta = MPDT[i][j].mu_ie * MPDT[i][j].ne * ME;
			m_eir      = eta * (MPDT[i][j].ver - MPDT[i][j].vir);
			m_eitheta  = eta * (MPDT[i][j].vetheta - MPDT[i][j].vitheta);
			m_eiz      = eta * (MPDT[i][j].vez - MPDT[i][j].viz);

			delta_ei = 3 * MPDT[i][j].ne * ME * MPDT[i][j].mu_ie * (MPDT[i][j].ee - MPDT[i][j].ei);

			double pre_eng = MPDT[i][j].ei + MPDT[i][j].ee;


			MPDT[i][j].ver -= m_eir;
			MPDT[i][j].vir += m_eir;
			MPDT[i][j].vetheta -= m_eitheta;
			MPDT[i][j].vitheta += m_eitheta;
			MPDT[i][j].vez -= m_eiz;
			MPDT[i][j].viz += m_eiz;

			MPDT[i][j].ver /=3;
			MPDT[i][j].vitheta /= 3;
			MPDT[i][j].vez /= 3;

			//能量守恒验证

			MPDT[i][j].ee = 0.5 * MPDT[i][j].ne * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
			MPDT[i][j].ei = 0.5 * MPDT[i][j].ni * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);

			//if (i == 200 && j == 11)
			//{
			//	printf("ss\n");
			//}
			//if (pre_eng - MPDT[i][j].ee - MPDT[i][j].ei < 1e-5)
			//{
			//	printf(" OK \n");
			//}
			//else 
			//{
			//	printf("pre_eng - MPDT[i][j].ee - MPDT[i][j].ei = %e\n", pre_eng - MPDT[i][j].ee - MPDT[i][j].ei);
			//}
		}
	}
}


void plasma_para()
{

}


void collision()
{
	
}