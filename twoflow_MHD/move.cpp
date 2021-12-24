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

double max_q_speed = 1e-5;

void move()
{
	double q_half;
	double hrho, htheta, hz, h2;
	double srho, stheta, sz;
	double urho, utheta, uz;
	double urho_half, utheta_half, uz_half;
	for (int i = 0; i < nz; i++)
	{
		for (int j = 0; j < nr; j++)
		{
			if ((ptype[i][j] & VACCUM_BOUNDARY) != 0) continue;
			if (MPDT[i][j].ne == 0 || MPDT[i][j].ni == 0) continue;
			if ((ptype[i][j] & ANODE_BOUNDARY) != 0 || (ptype[i][j] & CATHODE_BOUNDARY) != 0 ) continue;
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

			//电磁场加速之后电子流体动能密度
			double pre_ee =  0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);

			//if (MPDT[i][j].ver > 1.8e6)
			//{
			//	MPDT[i][j].ver = 1.8e6;
			//}
			//else if (MPDT[i][j].ver < -1.8e6)
			//{
			//	MPDT[i][j].ver = -1.8e6;
			//}

			//if (abs(MPDT[i][j].vetheta) > 1.8e6)
			//{
			//	MPDT[i][j].vetheta = 1.8e6;
			//}
			//else if (MPDT[i][j].vetheta < -1.8e6)
			//{
			//	MPDT[i][j].vetheta = -1.8e6;
			//}

			//if (MPDT[i][j].vez > 1.8e6)
			//{
			//	MPDT[i][j].vez = 1.8e6;

			//}
			//else if (MPDT[i][j].vez < -1.8e6)
			//{
			//	MPDT[i][j].vez = -1.8e6;
			//}


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



			double Q_ei = 1.4e-22;//电子离子碰撞截面
			//double sigma_q = PI * pow(QE, 4) / (128 * sqr(EPS_0) * sqr(ME) * pow((MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz), 4));
			MPDT[i][j].sigma_Q = PI * pow(QE, 4) / (128 * sqr(EPS_0) * sqr(ME) * pow((MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz), 4));
			MPDT[i][j].mu_ie = MPDT[i][j].ni * Q_ei * sqrt((12 * MPDT[i][j].ee) / (PI * ME));

			//if(abs(MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz) < 1e-6)
			//{ 
			//	MPDT[i][j].mu_ie = 0;
			//}
			//else
			//{
			//	MPDT[i][j].mu_ie = (MPDT[i][j].ne * PI * pow(QE, 4)) / (128 * sqr(EPS_0) * sqr(ME) * pow((MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz), 3));
			//}
			double m_eir = 0;        //电子向离子转移的径向（r方向）动量
			double m_eitheta = 0;    //电子向离子转移的角向（theta方向）动量
			double m_eiz = 0;        //电子向离子转移的轴向（z方向）动量

			//double delta_ei = 0;     //电子向离子转移的能量
			
			double eta = MPDT[i][j].mu_ie * MPDT[i][j].ne * ME;
			m_eir      = eta * (MPDT[i][j].ver - MPDT[i][j].vir);
			m_eitheta  = eta * (MPDT[i][j].vetheta - MPDT[i][j].vitheta);
			m_eiz      = eta * (MPDT[i][j].vez - MPDT[i][j].viz);

			//delta_ei = 3 * MPDT[i][j].ne * ME * MPDT[i][j].mu_ie * (MPDT[i][j].ee - MPDT[i][j].ei);

			double pre_eng = MPDT[i][j].ei + MPDT[i][j].ee;

			//if (MPDT[i][j].ne != 0 && MPDT[i][j].ni != 0)
			//{
			//	MPDT[i][j].ver -= m_eir / (MPDT[i][j].ne * ME);
			//	MPDT[i][j].vir += m_eir / (MPDT[i][j].ni * MI);
			//	MPDT[i][j].vetheta -= m_eitheta / (MPDT[i][j].ne * ME);
			//	MPDT[i][j].vitheta += m_eitheta / (MPDT[i][j].ni * MI);
			//	MPDT[i][j].vez -= m_eiz / (MPDT[i][j].ne * ME);
			//	MPDT[i][j].viz += m_eiz / (MPDT[i][j].ni * MI);
			//}

			if (MPDT[i][j].ni != 0)
			{
				MPDT[i][j].vir += m_eir / (MPDT[i][j].ni * MI);
				MPDT[i][j].vitheta += m_eitheta / (MPDT[i][j].ni * MI);
				MPDT[i][j].viz += m_eiz / (MPDT[i][j].ni * MI);
			}

			MPDT[i][j].ver = MPDT[i][j].vir;
			MPDT[i][j].vetheta = MPDT[i][j].vitheta;
			MPDT[i][j].vez = MPDT[i][j].viz;

			//离子速度和磁场夹角
			MPDT[i][j].angle_b_vi = magnetic_vec_angle(app_Br[i][j], app_Bz[i][j], MPDT[i][j].vir, MPDT[i][j].viz);

			//碰撞之后电子流体能量
			double alter_ee = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);

			double pre_ei = MPDT[i][j].ei;
			double ev = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez) ;
			MPDT[i][j].ei = 0.5 * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz) ;

			double dee = pre_ee - alter_ee;
			double dei = MPDT[i][j].ei - pre_ei;

			MPDT[i][j].pe += MPDT[i][j].ne * ME * (pre_ee - alter_ee) * (gamma - 1);
			MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne * ME) + 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
			MPDT[i][j].delta_ei = dei;
			

		}
	}
}

void move_q()
{
	double q_half;
	double hrho, htheta, hz, h2;
	double srho, stheta, sz;
	double urho, utheta, uz;
	double urho_half, utheta_half, uz_half;
	for (int i = 0; i < nz; i++)
	{
		for (int j = 0; j < nr; j++)
		{
			if ((ptype[i][j] & VACCUM_BOUNDARY) != 0) continue;
			if (MPDT[i][j].neq == 0 || MPDT[i][j].peq == 0) continue;

			//判断电子电荷和离子电荷是否分离
			if (is_electron_ion_separation(MPDT[i][j].angle_b_vi) == true)
			{
				//计算多余电子运动
				q_half = -dtq * QE / ME / 2;
				hrho = q_half * app_Br[i][j];
				htheta = 0;
				hz = q_half * app_Bz[i][j];
				h2 = sqr(hrho) + sqr(hz);

				srho = 2 * hrho / (1 + h2);
				stheta = 2 * htheta / (1 + h2);
				sz = 2 * hz / (1 + h2);
				urho = MPDT[i][j].vnqr + q_half * Er[i][j];
				utheta = MPDT[i][j].vnqtheta;
				uz = MPDT[i][j].vnqz + q_half * Ez[i][j];
				urho_half = urho + (utheta * sz - uz * stheta) + ((uz * hrho - urho * hz) * sz - (urho * htheta - utheta * hrho) * stheta);
				utheta_half = utheta + (uz * srho - urho * sz) + ((urho * htheta - utheta * hrho) * srho - (utheta * hz - uz * htheta) * sz);
				uz_half = uz + (urho * stheta - utheta * srho) + ((utheta * hz - uz * htheta) * stheta - (uz * hrho - urho * hz) * srho);
				MPDT[i][j].vnqr = urho_half + q_half * Er[i][j];
				MPDT[i][j].vnqtheta = utheta_half;
				MPDT[i][j].vnqz = uz_half + q_half * Ez[i][j];

				urho = MPDT[i][j].vnqr;
				utheta = MPDT[i][j].vnqtheta;
				uz = MPDT[i][j].vnqz;


				if (is_large_max_speed(urho, utheta, uz, max_q_speed) == true)
				{
					double tu = sqrt(sqr(urho) + sqr(utheta) + sqr(uz));
					MPDT[i][j].vnqr = urho / tu * max_q_speed;
					MPDT[i][j].vnqtheta = utheta / tu * max_q_speed;
					MPDT[i][j].vnqz = uz / tu * max_q_speed;
				}
			}
			else
			{
				if (is_large_max_speed(MPDT[i][j].vir, MPDT[i][j].vitheta, MPDT[i][j].viz, max_q_speed) == true)
				{
					double tu = sqrt(sqr(MPDT[i][j].vir) + sqr(MPDT[i][j].vitheta) + sqr(MPDT[i][j].viz));
					MPDT[i][j].vnqr = MPDT[i][j].vir / tu * max_q_speed;
					MPDT[i][j].vnqtheta = MPDT[i][j].vitheta / tu * max_q_speed;
					MPDT[i][j].vnqz = MPDT[i][j].viz / tu * max_q_speed;
				}
			}

			q_half = dtq * QE / ME / 2;
			hrho = q_half * app_Br[i][j];
			htheta = 0;
			hz = q_half * app_Bz[i][j];
			h2 = sqr(hrho) + sqr(hz);

			srho = 2 * hrho / (1 + h2);
			stheta = 2 * htheta / (1 + h2);
			sz = 2 * hz / (1 + h2);
			urho = MPDT[i][j].vpqr + q_half * Er[i][j];
			utheta = MPDT[i][j].vpqtheta;
			uz = MPDT[i][j].vpqz + q_half * Ez[i][j];
			urho_half = urho + (utheta * sz - uz * stheta) + ((uz * hrho - urho * hz) * sz - (urho * htheta - utheta * hrho) * stheta);
			utheta_half = utheta + (uz * srho - urho * sz) + ((urho * htheta - utheta * hrho) * srho - (utheta * hz - uz * htheta) * sz);
			uz_half = uz + (urho * stheta - utheta * srho) + ((utheta * hz - uz * htheta) * stheta - (uz * hrho - urho * hz) * srho);
			MPDT[i][j].vpqr = urho_half + q_half * Er[i][j];
			MPDT[i][j].vpqtheta = utheta_half;
			MPDT[i][j].vpqz = uz_half + q_half * Ez[i][j];

			
			urho = MPDT[i][j].vpqr;
			utheta = MPDT[i][j].vpqtheta;
			uz = MPDT[i][j].vpqz;

			if (is_large_max_speed(urho, utheta, uz, max_q_speed) == true)
			{
				double tu = sqrt(sqr(urho) + sqr(utheta) + sqr(uz));
				MPDT[i][j].vpqr = urho / tu * max_q_speed;
				MPDT[i][j].vpqtheta = utheta / tu * max_q_speed;
				MPDT[i][j].vpqz = uz / tu * max_q_speed;
			}
			
		}
	}
}

void plasma_para()
{

}


void collision()
{
	
}

double  magnetic_vec_angle(double var,double vaz,double vbr,double vbz)
{
	double tmp = (var * vbr + vaz * vbz) / (sqrt(sqr(var) + sqr(vaz)) * sqrt(sqr(vbr) + sqr(vbz)));
	return acos(tmp) * 180 / PI;
}

bool is_electron_ion_separation(double angle)
{
	
	return angle > 10;
}

bool is_large_max_speed(double ur, double utheta, double uz, double max_speed)
{
	return sqr(ur) + sqr(utheta) + sqr(uz) > sqr(max_speed);
}