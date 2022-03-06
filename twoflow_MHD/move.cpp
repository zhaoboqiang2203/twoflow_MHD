#include "MHD.h"

/*************************************************

Copyright (c) ,2015-2016,  
File name:
Author:   

Date:2010-08-25

Description:������Ҫʵ�ֵĹ���
Fucntion list:...

Histiory : //�޸���ʷ��¼�б�ÿ���޸ļ�¼�޸����ڣ��޸��ߣ��Լ��޸�����

1.Data:....
  Author:...
  Midification:...

2.


**************************************************/



/*************************************************

Function:                  // ��������

Description:               // hanshu ogngnegn 

Calls:                     // ��

Input:                     // �������˵��������ÿ����������

						   // �á�ȡֵ˵�����������ϵ��

Output:                    // �����������˵����

Return:                    // ��������ֵ��˵��

Others:                   // ����˵��

*************************************************/

double max_q_speed = 1e-5;

void move()
{
	double q_half;
	double hrho, htheta, hz, h2;
	double srho, stheta, sz;
	double urho, utheta, uz;
	double urho_half, utheta_half, uz_half;
	double Epr, Epz;
	double t0;
	double ct;

	register int i, j, k;

	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			//if (i == 390 && j == 0)
			//{
			//	printf("ptype[i][j] = %d\n", ptype[i][j]);
			//}
			if ((ptype[i][j] & VACCUM_BOUNDARY) != 0) continue;
			if (MPDT[i][j].ne == 0 || MPDT[i][j].ni == 0) continue;
			if ((ptype[i][j] & ANODE_BOUNDARY) != 0 || (ptype[i][j] & CATHODE_BOUNDARY) != 0 ) continue;

			if (sqr(app_Br[i][j]) + sqr(app_Bz[i][j]) == 0)
			{
				ct = 0;
			}
			else
			{
				ct = (app_Br[i][j] * Er[i][j] + Ez[i][j] * app_Bz[i][j]) / (sqr(app_Br[i][j]) + sqr(app_Bz[i][j]));
			}
			
			Epr = Er[i][j] - ct * app_Br[i][j];
			Epz = Ez[i][j] - ct * app_Bz[i][j];

			t0 = ME / (QE * sqrt(sqr(app_Br[i][j]) + sqr(app_Bz[i][j])));

			//fmod(dt, t0);
			//q_half = -dt * QE / ME / 2;

			q_half = -fmod(dt, t0) * QE / ME / 2;
			hrho = q_half * app_Br[i][j];
			htheta = 0;
			hz = q_half * app_Bz[i][j];
			h2 = sqr(hrho) + sqr(hz);
			
			srho = 2 * hrho / (1 + h2);
			stheta = 2 * htheta / (1 + h2); 
			sz = 2 * hz / (1 + h2);
			urho = MPDT[i][j].ver + q_half * Epr;
			utheta = MPDT[i][j].vetheta;
			uz = MPDT[i][j].vez + q_half * Epz;
			urho_half = urho + (utheta * sz - uz * stheta) + ((uz * hrho - urho * hz) * sz - (urho * htheta - utheta * hrho) * stheta);
			utheta_half = utheta + (uz * srho - urho * sz) + ((urho * htheta - utheta * hrho) * srho - (utheta * hz - uz * htheta) * sz);
			uz_half = uz + (urho * stheta - utheta * srho) + ((utheta * hz - uz * htheta) * stheta - (uz * hrho - urho * hz) * srho);
			MPDT[i][j].ver = urho_half + q_half * Epr;
			MPDT[i][j].vetheta = utheta_half;
			MPDT[i][j].vez = uz_half + q_half * Epz;
		
			MPDT[i][j].ver += -QE * ct * app_Br[i][j] / ME * dt;
			MPDT[i][j].vez += -QE * ct * app_Bz[i][j] / ME * dt;






			t0 = MI / (QE * sqrt(sqr(app_Br[i][j]) + sqr(app_Bz[i][j])));
			
			q_half = fmod(dt, t0) * QE / MI / 2;
			hrho = q_half * app_Br[i][j];
			htheta = 0;
			hz = q_half * app_Bz[i][j];
			h2 = sqr(hrho) + sqr(hz);

			srho = 2 * hrho / (1 + h2);
			stheta = 2 * htheta / (1 + h2);
			sz = 2 * hz / (1 + h2);
			urho = MPDT[i][j].vir + q_half * Epr;
			utheta = MPDT[i][j].vitheta;
			uz = MPDT[i][j].viz + q_half * Epz;
			urho_half = urho + (utheta * sz - uz * stheta) + ((uz * hrho - urho * hz) * sz - (urho * htheta - utheta * hrho) * stheta);
			utheta_half = utheta + (uz * srho - urho * sz) + ((urho * htheta - utheta * hrho) * srho - (utheta * hz - uz * htheta) * sz);
			uz_half = uz + (urho * stheta - utheta * srho) + ((utheta * hz - uz * htheta) * stheta - (uz * hrho - urho * hz) * srho);
			MPDT[i][j].vir = urho_half + q_half * Epr;
			MPDT[i][j].vitheta = utheta_half;
			MPDT[i][j].viz = uz_half + q_half * Epz;
			
			MPDT[i][j].vir += QE * ct * app_Br[i][j] / MI * dt;
			MPDT[i][j].viz += QE * ct * app_Bz[i][j] / MI * dt;

			MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne * ME) + 0.5 * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
			MPDT[i][j].ei = MPDT[i][j].pi / (gamma - 1) / (MPDT[i][j].ni * MI) + 0.5 * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);


			ionization_collisions(i, j);
			coulomb_collision(i, j);
			
			//viscosity_collision(i, j);
		//	double pre_ei = 0.5 * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);



		//	double Q_ei = 1.4e-22;//����������ײ����
		//	MPDT[i][j].tau_ei = pow(pre_ee,1.5) * sqrt(ME) / (11.313708 * PI * MPDT[i][j].ni * pow(QE, 4) * 10);
		//	//double sigma_q = PI * pow(QE, 4) / (128 * sqr(EPS_0) * sqr(ME) * pow((MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz), 4));
		//	MPDT[i][j].sigma_Q = PI * pow(QE, 4) / (128 * sqr(EPS_0) * sqr(ME) * pow((MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz), 4));
		//	MPDT[i][j].mu_ie = MPDT[i][j].ni * Q_ei * sqrt((12 * MPDT[i][j].ee) / (PI * ME));

		//	//if(abs(MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz) < 1e-6)
		//	//{ 
		//	//	MPDT[i][j].mu_ie = 0;
		//	//}
		//	//else
		//	//{
		//	//	MPDT[i][j].mu_ie = (MPDT[i][j].ne * PI * pow(QE, 4)) / (128 * sqr(EPS_0) * sqr(ME) * pow((MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz), 3));
		//	//}
		//	double m_eir = 0;        //����������ת�Ƶľ���r���򣩶���
		//	double m_eitheta = 0;    //����������ת�ƵĽ���theta���򣩶���
		//	double m_eiz = 0;        //����������ת�Ƶ�����z���򣩶���

		//	MPDT[i][j].delta_ei = dt / MPDT[i][j].tau_ei * (pre_ee - pre_ei);     //����������ת�Ƶ�����
		//	
		//	//double eta = MPDT[i][j].mu_ie * MPDT[i][j].ne * ME;
		//	//m_eir      = eta * (MPDT[i][j].ver - MPDT[i][j].vir);
		//	//m_eitheta  = eta * (MPDT[i][j].vetheta - MPDT[i][j].vitheta);
		//	//m_eiz      = eta * (MPDT[i][j].vez - MPDT[i][j].viz);

		//	////delta_ei = 3 * MPDT[i][j].ne * ME * MPDT[i][j].mu_ie * (MPDT[i][j].ee - MPDT[i][j].ei);

		//	//double pre_eng = MPDT[i][j].ei + MPDT[i][j].ee;

		//	//if (MPDT[i][j].ne != 0 && MPDT[i][j].ni != 0)
		//	//{
		//	//	MPDT[i][j].ver -= m_eir / (MPDT[i][j].ne * ME);
		//	//	MPDT[i][j].vir += m_eir / (MPDT[i][j].ni * MI);
		//	//	MPDT[i][j].vetheta -= m_eitheta / (MPDT[i][j].ne * ME);
		//	//	MPDT[i][j].vitheta += m_eitheta / (MPDT[i][j].ni * MI);
		//	//	MPDT[i][j].vez -= m_eiz / (MPDT[i][j].ne * ME);
		//	//	MPDT[i][j].viz += m_eiz / (MPDT[i][j].ni * MI);
		//	//}
		//	tur = sqr(MPDT[i][j].ver - MPDT[i][j].vir);
		//	tutheta = sqr(MPDT[i][j].vetheta - MPDT[i][j].vitheta);
		//	tuz = sqr(MPDT[i][j].vez - MPDT[i][j].viz);

		//	tu2 = tuz + tutheta + tuz;

		//	MPDT[i][j].vir += sqrt(2 * tur / tu2 * MPDT[i][j].delta_ei / (MPDT[i][j].ni * MI));
		//	MPDT[i][j].vitheta += sqrt(2 * tutheta / tu2 * MPDT[i][j].delta_ei / (MPDT[i][j].ni * MI));
		//	MPDT[i][j].viz += sqrt(2 * tuz / tu2 * MPDT[i][j].delta_ei / (MPDT[i][j].ni * MI));

		//	MPDT[i][j].ver -= sqrt(2 * tur / tu2 * 11 * MPDT[i][j].delta_ei / (MPDT[i][j].ne * ME));
		//	MPDT[i][j].vetheta -= sqrt(2 * tutheta / tu2 * 11 * MPDT[i][j].delta_ei / (MPDT[i][j].ne * ME));
		//	MPDT[i][j].vez -= sqrt(2 * tuz / tu2 * 11 * MPDT[i][j].delta_ei / (MPDT[i][j].ne * ME));

		//	//MPDT[i][j].ver = MPDT[i][j].vir;
		//	//MPDT[i][j].vetheta = MPDT[i][j].vitheta;
		//	//MPDT[i][j].vez = MPDT[i][j].viz;

		//	//�����ٶȺʹų��н�
	MPDT[i][j].angle_b_vi = magnetic_vec_angle(app_Br[i][j], app_Bz[i][j], MPDT[i][j].vir, MPDT[i][j].viz);

		//	////��ײ֮�������������
		//	//double alter_ee = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);

		//	//double pre_ei = MPDT[i][j].ei;
		//	//double ev = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez) ;
		//	MPDT[i][j].ei = 0.5 * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz) ;

		//	//double dee = pre_ee - alter_ee;
		//	//double dei = MPDT[i][j].ei - pre_ei;


		//	MPDT[i][j].pe += MPDT[i][j].ne * 10 * MPDT[i][j].delta_ei * (gamma - 1);
		//	MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne) + 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);




		//	//���벿�ִ���ɣ������������ڵ���ֵ���ӵ��������ܶ�

		//	if (MPDT[i][j].ee > 15.6 * QE)
		//	{
		//		if (MPDT[i][j].ne < 1e21)
		//		{
		//			MPDT[i][j].ni += MPDT[i][j].ne;
		//			MPDT[i][j].ne *= 2;
		//		}

		//		double tep = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne);
		//		double teu = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);

		//		double tep1 = tep / (tep + teu) * 15.6 * QE;
		//		double teu2 = teu / (tep + teu) * 15.6 * QE;


		//		tur = sqr(MPDT[i][j].ver);
		//		tutheta = sqr(MPDT[i][j].vetheta);
		//		tuz = sqr(MPDT[i][j].vez);

		//		tu2 = tuz + tutheta + tuz;

		//		MPDT[i][j].ver -= sqrt(2 * tur / tu2 * teu2 /(MPDT[i][j].ne * ME));
		//		MPDT[i][j].vetheta -= sqrt(2 * tutheta / tu2 * teu2 / (MPDT[i][j].ne * ME));
		//		MPDT[i][j].vez -= sqrt(2 * tuz / tu2 * teu2 / (MPDT[i][j].ne * ME));

		//		MPDT[i][j].pe -= MPDT[i][j].ne * tep1 * (gamma - 1);
		//		MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne) + 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
		//	}
		}
	}


	//for (i = 0; i < nz; i++)
	//{
	//	for (j = 0; j < nr; j++)
	//	{

	//		if (btype[i][j] == 1)
	//		{
	//			MPDT[i][j].vir = vir[i][j];
	//			MPDT[i][j].viz = viz[i][j];

	//			MPDT[i][j].ver = ver[i][j];
	//			MPDT[i][j].vez = vez[i][j];
	//		}
	//		else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
	//		{
	//			MPDT[i][j].vir = vir[i][j];
	//			MPDT[i][j].viz = viz[i][j];

	//			MPDT[i][j].ver = ver[i][j];
	//			MPDT[i][j].vez = vez[i][j];
	//		}
	//		else
	//		{

	//		}
	//	}
	//}
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
			if (MPDT[i][j].neq != 0 ) 
			{
				//�жϵ��ӵ�ɺ����ӵ���Ƿ����
				if (is_electron_ion_separation(MPDT[i][j].angle_b_vi) == true)
				{
					//�����������˶�
					q_half = -dt * QE / ME / 2;
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
					else
					{
						MPDT[i][j].vnqr = MPDT[i][j].vir;
						MPDT[i][j].vnqtheta = MPDT[i][j].vitheta;
						MPDT[i][j].vnqz = MPDT[i][j].viz;
					}
				}
			}

			if (MPDT[i][j].peq != 0)
			{
				q_half = dt * QE / ME / 2;
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
	return true;
	//return angle > 10;
}

bool is_large_max_speed(double ur, double utheta, double uz, double max_speed)
{
	return sqr(ur) + sqr(utheta) + sqr(uz) > sqr(max_speed);
}

/*************************************************

Function:                  // ionization_collisions

Description:               // ԭ�ӵ��뺯��

Calls:                     // ��

Input:                     // int i int j �������������λ��

						   // �á�ȡֵ˵�����������ϵ��

Output:                    // �����������˵����

Return:                    // ��������ֵ��˵��

Others:                    // ����˵��

*************************************************/


void ionization_collisions(int i, int j)
{
	double tur;
	double tutheta;
	double tuz;

	double tu2, teu2;
	//double nn = (1e21 - MPDT[i][j].ne - MPDT[i][j].neq + MPDT[i][j].peq);
	double nn = (1e21 - MPDT[i][j].ni);
	if (nn < 0)
	{
		return;
	}

	if (MPDT[i][j].ne <= 0)
	{
		return;
	}

	double nee = MPDT[i][j].ee;

	if (nee < 0)
	{
		return;
	}

	double p = 1 - exp( -nn * sqrt(2 * nee ) * 1e-20 * dt);
	MPDT[i][j].mu_ie = p;
	if (MPDT[i][j].ee * ME - 15.755 * QE * p > 1e-25)
	{
		//MPDT[i][j].ni += nn * p;
		//MPDT[i][j].ne += nn * p;

		//MPDT[i][j].ee -= 15.755 * QE * p;
		teu2 = 15.755 * QE * p;
		tur = sqr(MPDT[i][j].ver);
		tutheta = sqr(MPDT[i][j].vetheta);
		tuz = sqr(MPDT[i][j].vez);

		tu2 = tuz + tutheta + tuz;

		MPDT[i][j].ver -= sqrt(2 * tur / tu2 * teu2 / ME);
		MPDT[i][j].vetheta -= sqrt(2 * tutheta / tu2 * teu2 / ME);
		MPDT[i][j].vez -= sqrt(2 * tuz / tu2 * teu2 / ME);

		//double add_n = nn * p;

		//MPDT[i][j].ver = MPDT[i][j].ver * MPDT[i][j].ne / (MPDT[i][j].ne + add_n);
		//MPDT[i][j].vetheta = MPDT[i][j].vetheta * MPDT[i][j].ne / (MPDT[i][j].ne + add_n);
		//MPDT[i][j].vez = MPDT[i][j].vez * MPDT[i][j].ne / (MPDT[i][j].ne + add_n);

		//MPDT[i][j].vir = MPDT[i][j].vir * MPDT[i][j].ne / (MPDT[i][j].ni + add_n);
		//MPDT[i][j].vitheta = MPDT[i][j].vitheta * MPDT[i][j].ne / (MPDT[i][j].ni + add_n);
		//MPDT[i][j].viz = MPDT[i][j].viz * MPDT[i][j].ne / (MPDT[i][j].ni + add_n);

		MPDT[i][j].ni += nn * p;
		MPDT[i][j].ne += nn * p;

		//MPDT[i][j].pe -= MPDT[i][j].ne * tep1 * (gamma - 1);
		MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne * ME) + 0.5 * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
		//printf("ee pre = %e,after = %e\n", nee, MPDT[i][j].ee);
		MPDT[i][j].ei = MPDT[i][j].pi / (gamma - 1) / (MPDT[i][j].ni * MI) + 0.5 * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);
	}

	/*if (MPDT[i][j].ee > 15.755 * QE)
	{

		MPDT[i][j].ni += MPDT[i][j].ne;
		MPDT[i][j].ne *= 2;

		double tep = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne);
		double teu = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);

		double tep1 = tep / (tep + teu) * 15.6 * QE;
		double teu2 = teu / (tep + teu) * 15.6 * QE;


		tur = sqr(MPDT[i][j].ver);
		tutheta = sqr(MPDT[i][j].vetheta);
		tuz = sqr(MPDT[i][j].vez);

		tu2 = tuz + tutheta + tuz;

		MPDT[i][j].ver -= sqrt(2 * tur / tu2 * teu2 / (MPDT[i][j].ne * ME));
		MPDT[i][j].vetheta -= sqrt(2 * tutheta / tu2 * teu2 / (MPDT[i][j].ne * ME));
		MPDT[i][j].vez -= sqrt(2 * tuz / tu2 * teu2 / (MPDT[i][j].ne * ME));

		MPDT[i][j].pe -= MPDT[i][j].ne * tep1 * (gamma - 1);
		MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne) + 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
	}*/

	return;
}

/// <summary>
/// ������ײ�����е������Ӷ����Ͷ���ת��
/// </summary>
/// <param name="i"></param>
/// <param name="j"></param>
/// <returns></returns>
void coulomb_collision(int i, int j)
{
	double ur, utheta, uz;
	double ur1, utheta1, uz1;
	double pre_ee, pre_ei;
	double ep;
	double ne = MPDT[i][j].ne;

	pre_ee = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
	pre_ei = 0.5 * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);

	ur = (MPDT[i][j].ver * ne * ME + MPDT[i][j].vir * MPDT[i][j].ni * MI) / (ne * ME + MPDT[i][j].ni * MI);
	utheta = (MPDT[i][j].vetheta * ne * ME + MPDT[i][j].vitheta * MPDT[i][j].ni * MI) / (ne * ME + MPDT[i][j].ni * MI);
	uz = (MPDT[i][j].vez * ne * ME + MPDT[i][j].viz * MPDT[i][j].ni * MI) / (ne * ME + MPDT[i][j].ni * MI);

	ep = pre_ee + pre_ei - 0.5 * (ME + MI) * (sqr(ur) + sqr(utheta) + sqr(uz));

	if (ep < 0)
	{
		ep = 0;
	}
	

	MPDT[i][j].pe = ne * ep * (gamma - 1);

	//ur1 = sqrt((sqr(ur) + sqr(utheta) + sqr(uz)) / 3.0);

	//MPDT[i][j].vir = ur + 0.01 * (ur1 - ur);
	//MPDT[i][j].vitheta = utheta + 0.01 * (ur1 - utheta);
	//MPDT[i][j].viz = uz + 0.01 * (ur1 - uz);

	//MPDT[i][j].ver = MPDT[i][j].vir;
	//MPDT[i][j].vetheta = MPDT[i][j].vitheta;
	//MPDT[i][j].vez = MPDT[i][j].viz;

	MPDT[i][j].vir = ur;
	MPDT[i][j].vitheta = utheta;
	MPDT[i][j].viz = uz;

	MPDT[i][j].ver = ur;
	MPDT[i][j].vetheta = utheta;
	MPDT[i][j].vez = uz;
	//
	//MPDT[i][j].ee += ep;


	MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne * ME) + 0.5 * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
	if (MPDT[i][j].ee > 1807 * QE / ME)
	{
		MPDT[i][j].ee = 1807 * QE / ME;
	}
	MPDT[i][j].ei = MPDT[i][j].pi / (gamma - 1) / (MPDT[i][j].ni * MI) + 0.5 * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);

	return;
}

void viscosity_collision(int i, int j)
{
	if (btype[i][j] == 1)
	{
		vir[i][j] = MPDT[i][j].vir + 0.1 * ((MPDT[i + 1][j].vir + MPDT[i][j].vir + MPDT[i - 1][j].vir) / 3 - MPDT[i][j].vir);
		viz[i][j] = MPDT[i][j].viz + 0.1 * ((MPDT[i][j + 1].viz + MPDT[i][j].viz + MPDT[i][j - 1].viz) / 3 - MPDT[i][j].viz);

		ver[i][j] = MPDT[i][j].ver + 0.1 * ((MPDT[i + 1][j].ver + MPDT[i][j].ver + MPDT[i - 1][j].ver) / 3 - MPDT[i][j].ver);
		vez[i][j] = MPDT[i][j].vez + 0.1 * ((MPDT[i][j + 1].vez + MPDT[i][j].vez + MPDT[i][j - 1].vez) / 3 - MPDT[i][j].vez);
	}
	else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
	{
		vir[i][j] = MPDT[i][j].vir + 0.1 * ((MPDT[i][j + 1].vir + MPDT[i][j].vir + MPDT[i][j - 1].vir) / 3 - MPDT[i][j].vir);
		viz[i][j] = MPDT[i][j].viz + 0.1 * ((MPDT[i + 1][j].viz + MPDT[i][j].viz ) / 2 - MPDT[i][j].viz);

		ver[i][j] = MPDT[i][j].ver + 0.1 * ((MPDT[i][j + 1].ver + MPDT[i][j].ver + MPDT[i][j - 1].ver) / 3 - MPDT[i][j].ver);
		vez[i][j] = MPDT[i][j].vez + 0.1 * ((MPDT[i + 1][j].vez + MPDT[i][j].vez ) / 2 - MPDT[i][j].vez);
	}
	else 
	{
		
	}
}

void recombination_collision()
{
	
}

void radiation()
{

}

//void current_caulate()
//{
//	int start_r = 0, end_r = 0, pos_z = 0;
//	int start_z = 0, end_z = 0, pos_r = 0;
//	register int i, j;
//
//	pos_r = 34;
//	start_z = 0;
//	end_z = 382;
//
//	pos_z = 382;
//	start_r = 0;
//	end_r = 34;
//
//	current_I = 0;
//
//	j = pos_r;
//
//	for(i = start_z; i < end_z; i++)
//	{
//		current_I += (MPDT[i][j].peq * MPDT[i][j].vpqr - MPDT[i][j].neq * MPDT[i][j].vnqr) * QE  * dr * 2 * PI * j * dz ;
//		//current_I += (MPDT[i][j].peq * MPDT[i][j].vpqz - MPDT[i][j].neq * MPDT[i][j].vnqz) * QE ;
//	}
//	i = pos_z;
//
//	for (j = start_r; j < end_r; j++)
//	{
//		//current_I += (MPDT[i][j].peq * MPDT[i][j].vpqr - MPDT[i][j].neq * MPDT[i][j].vnqr) * QE 
//		if (j != 0) 
//		{
//			current_I += (MPDT[i][j].peq * MPDT[i][j].vpqz - MPDT[i][j].neq * MPDT[i][j].vnqz) * QE  * sqr(dr) * dtheta * (2 * j - 1) ;
//		}
//		
//	}
//
//	printf("current = %lf\n", current_I);
//}
//
//
//void current_control()
//{
//	pid.err = pid.set_current - (-current_I);
//	pid.integral += pid.err;
//
//	if (pid.err > 50 || pid.err < -50)
//	{
//		pid.ne_density = pid.Kp * pid.err + pid.Kd * (pid.err - pid.err_last);
//	}
//	else
//	{
//		pid.ne_density = pid.Kp * pid.err + pid.Ki * pid.integral + pid.Kd * (pid.err - pid.err_last);
//	}
//
//	pid.err_last = pid.err;
//	inter_e_den += pid.ne_density;
//
//	if (inter_e_den < 0)
//	{
//		inter_e_den = 0;
//	}
//	out_e_den = inter_e_den * cathod_cell / anode_cell;
//
//}

void current_caulate()
{
	int start_r = 0, end_r = 0, pos_z = 0;
	int start_z = 0, end_z = 0, pos_r = 0;
	register int i, j;

	pos_r = 34;
	start_z = 0;
	end_z = 382;

	pos_z = 382;
	start_r = 0;
	end_r = 34;

	current_I = 0;

	j = pos_r;

	for (i = start_z; i < end_z; i++)
	{
		current_I += (MPDT[i][j].ni * MPDT[i][j].vir - MPDT[i][j].ne * MPDT[i][j].ver) * QE * dr * 2 * PI * j * dz;
	}
	i = pos_z;

	for (j = start_r; j < end_r; j++)
	{
		if (j != 0)
		{
			current_I += (MPDT[i][j].ni * MPDT[i][j].viz - MPDT[i][j].ne * MPDT[i][j].vez) * QE * sqr(dr) * PI * (2 * j - 1);
		}

	}

	printf("current = %lf\n", current_I);
}


void current_control()
{
	pid.err = pid.set_current - (-current_I);
	pid.integral += pid.err;

	if (pid.err > 50 || pid.err < -50)
	{
		pid.ne_density = pid.Kp * pid.err + pid.Kd * (pid.err - pid.err_last);
	}
	else
	{
		pid.ne_density = pid.Kp * pid.err + pid.Ki * pid.integral + pid.Kd * (pid.err - pid.err_last);
	}

	pid.err_last = pid.err;
	inter_e_den += pid.ne_density;

	if (inter_e_den < 0)
	{
		inter_e_den = 0;
	}
	out_e_den = inter_e_den * cathod_cell / anode_cell;

}

