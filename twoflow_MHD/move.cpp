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

double e_half[ZMAX][RMAX];
double e_hrho[ZMAX][RMAX], e_htheta[ZMAX][RMAX], e_hz[ZMAX][RMAX], e_h2[ZMAX][RMAX];
double e_srho[ZMAX][RMAX], e_stheta[ZMAX][RMAX], e_sz[ZMAX][RMAX];

double i_half[ZMAX][RMAX];
double i_hrho[ZMAX][RMAX], i_htheta[ZMAX][RMAX], i_hz[ZMAX][RMAX], i_h2[ZMAX][RMAX];
double i_srho[ZMAX][RMAX], i_stheta[ZMAX][RMAX], i_sz[ZMAX][RMAX];

int cur_pos_r;
int cur_pos_z;

int ion_area_rmin;
int ion_area_rmaxn;
int ion_area_zmin;
int ion_area_zmax;

void move()
{
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

			//t0 = 2 * PI * ME / (QE * sqrt(sqr(app_Br[i][j]) + sqr(app_Bz[i][j])));

			//fmod(dt, t0);
			//q_half = -dt * QE / ME / 2;
			//MPDT[i][j].angle_b_vi = fmod(dt, t0);
			//q_half = -fmod(dt, t0) * QE / ME / 2;
			//hrho = q_half * app_Br[i][j];
			//htheta = 0;
			//hz = q_half * app_Bz[i][j];
			//h2 = sqr(hrho) + sqr(hz);
			//
			//srho = 2 * hrho / (1 + h2);
			//stheta = 2 * htheta / (1 + h2); 
			//sz = 2 * hz / (1 + h2);
			urho = MPDT[i][j].ver + e_half[i][j] * Epr;
			utheta = MPDT[i][j].vetheta;
			uz = MPDT[i][j].vez + e_half[i][j] * Epz;
			urho_half = urho + (utheta * e_sz[i][j] - uz * e_stheta[i][j]) + ((uz * e_hrho[i][j] - urho * e_hz[i][j]) * e_sz[i][j] - (urho * e_htheta[i][j] - utheta * e_hrho[i][j]) * e_stheta[i][j]);
			utheta_half = utheta + (uz * e_srho[i][j] - urho * e_sz[i][j]) + ((urho * e_htheta[i][j] - utheta * e_hrho[i][j]) * e_srho[i][j] - (utheta * e_hz[i][j] - uz * e_htheta[i][j]) * e_sz[i][j]);
			uz_half = uz + (urho * e_stheta[i][j] - utheta * e_srho[i][j]) + ((utheta * e_hz[i][j] - uz * e_htheta[i][j]) * e_stheta[i][j] - (uz * e_hrho[i][j] - urho * e_hz[i][j]) * e_srho[i][j]);
			MPDT[i][j].ver = urho_half + e_half[i][j] * Epr;
			MPDT[i][j].vetheta = utheta_half;
			MPDT[i][j].vez = uz_half + e_half[i][j] * Epz;
		
			MPDT[i][j].ver += -QE * ct * app_Br[i][j] / ME * dt;
			MPDT[i][j].vez += -QE * ct * app_Bz[i][j] / ME * dt;



			//t0 = 2 * PI * MI / (QE * sqrt(sqr(app_Br[i][j]) + sqr(app_Bz[i][j])));
			//
			//q_half = fmod(dt, t0) * QE / MI / 2;
			//hrho = q_half * app_Br[i][j];
			//htheta = 0;
			//hz = q_half * app_Bz[i][j];
			//h2 = sqr(hrho) + sqr(hz);

			//srho = 2 * hrho / (1 + h2);
			//stheta = 2 * htheta / (1 + h2);
			//sz = 2 * hz / (1 + h2);
			urho = MPDT[i][j].vir + i_half[i][j] * Epr;
			utheta = MPDT[i][j].vitheta;
			uz = MPDT[i][j].viz + i_half[i][j] * Epz;
			urho_half = urho + (utheta * i_sz[i][j] - uz * i_stheta[i][j]) + ((uz * i_hrho[i][j] - urho * i_hz[i][j]) * i_sz[i][j] - (urho * i_htheta[i][j] - utheta * i_hrho[i][j]) * i_stheta[i][j]);
			utheta_half = utheta + (uz * i_srho[i][j] - urho * i_sz[i][j]) + ((urho * i_htheta[i][j] - utheta * i_hrho[i][j]) * i_srho[i][j] - (utheta * i_hz[i][j] - uz * i_htheta[i][j]) * i_sz[i][j]);
			uz_half = uz + (urho * i_stheta[i][j] - utheta * i_srho[i][j]) + ((utheta * i_hz[i][j] - uz * i_htheta[i][j]) * i_stheta[i][j] - (uz * i_hrho[i][j] - urho * i_hz[i][j]) * i_srho[i][j]);
			MPDT[i][j].vir = urho_half + i_half[i][j] * Epr;
			MPDT[i][j].vitheta = utheta_half;
			MPDT[i][j].viz = uz_half + i_half[i][j] * Epz;
			
			MPDT[i][j].vir += QE * ct * app_Br[i][j] / MI * dt;
			MPDT[i][j].viz += QE * ct * app_Bz[i][j] / MI * dt;

			MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne * ME) + 0.5 * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
			MPDT[i][j].ei = MPDT[i][j].pi / (gamma - 1) / (MPDT[i][j].ni * MI) + 0.5 * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);


			ionization_collisions(i, j);
			coulomb_collision(i, j);
			
			//viscosity_collision(i, j);

		//	//离子速度和磁场夹角
	//MPDT[i][j].angle_b_vi = magnetic_vec_angle(app_Br[i][j], app_Bz[i][j], MPDT[i][j].vir, MPDT[i][j].viz);


		}
	}


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

Description:               // 原子电离函数

Calls:                     // 无

Input:                     // int i int j 在流场中网格点位置

						   // 用、取值说明及参数间关系。

Output:                    // 对输出参数的说明。

Return:                    // 函数返回值的说明

Others:                    // 其它说明

*************************************************/


void ionization_collisions(int i, int j)
{
	double tur;
	double tutheta;
	double tuz;

	double tu2, teu2;
	//double nn = (1e21 - MPDT[i][j].ne - MPDT[i][j].neq + MPDT[i][j].peq);
	double nn = (atom[i][j].den - MPDT[i][j].ni);
	double tden = 0;
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

	double  ndot = MPDT[i][j].ne * ME * MPDT[i][j].ee * 0.16 / (Ionization_Energy * QE);

	if (ndot > nn)
	{
		ndot = nn;
	}
	
	tur = 0.4 * MPDT[i][j].ver;
	tutheta = 0.4 * MPDT[i][j].vetheta;
	tuz = 0.4 * MPDT[i][j].vez;


	MPDT[i][j].ver -= tur;
	MPDT[i][j].vetheta -= tutheta;
	MPDT[i][j].vez -= tuz;

	tden = ndot;
	MPDT[i][j].ver = MPDT[i][j].ver * MPDT[i][j].ne / (MPDT[i][j].ne + tden);
	MPDT[i][j].vetheta = MPDT[i][j].vetheta * MPDT[i][j].ne / (MPDT[i][j].ne + tden);
	MPDT[i][j].vez = MPDT[i][j].vez * MPDT[i][j].ne / (MPDT[i][j].ne + tden);

	MPDT[i][j].vir = MPDT[i][j].vir * MPDT[i][j].ni / (MPDT[i][j].ni + tden);
	MPDT[i][j].vitheta = MPDT[i][j].vitheta * MPDT[i][j].ni / (MPDT[i][j].ni + tden);
	MPDT[i][j].viz = MPDT[i][j].viz * MPDT[i][j].ni / (MPDT[i][j].ni + tden);

		
		

	MPDT[i][j].ni += tden;
	MPDT[i][j].ne += tden;
		

		

	//MPDT[i][j].pe -= MPDT[i][j].ne * tep1 * (gamma - 1);
	MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne * ME) + 0.5 * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
	//printf("ee pre = %e,after = %e\n", nee, MPDT[i][j].ee);
	MPDT[i][j].ei = MPDT[i][j].pi / (gamma - 1) / (MPDT[i][j].ni * MI) + 0.5 * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);


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
/// 库仑碰撞，进行电子离子动量和动能转移
/// </summary>
/// <param name="i"></param>
/// <param name="j"></param>
/// <returns></returns>
void coulomb_collision(int i, int j)
{
	double ur, utheta, uz;
	double ur2, utheta2, uz2;
	double upre, upre1, upre2;
	double pre_ee, pre_ei;
	double ep;
	double ne = 0;

	if (MPDT[i][j].ne > MPDT[i][j].ni * 10)
	{
		ne = MPDT[i][j].ni * 10;
	}
	else
	{
		ne = MPDT[i][j].ne;
	}


	pre_ee = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
	pre_ei = 0.5 * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);

	ur = (MPDT[i][j].ver * ne * ME + MPDT[i][j].vir * MPDT[i][j].ni * MI) / (ne * ME + MPDT[i][j].ni * MI);
	utheta = (MPDT[i][j].vetheta * ne * ME + MPDT[i][j].vitheta * MPDT[i][j].ni * MI) / (ne * ME + MPDT[i][j].ni * MI);
	uz = (MPDT[i][j].vez * ne * ME + MPDT[i][j].viz * MPDT[i][j].ni * MI) / (ne * ME + MPDT[i][j].ni * MI);

	ur2 = sqr(ur);
	utheta2 = sqr(utheta);
	uz2 = sqr(uz);

	upre2 = ur2 + utheta2 + uz2;
	upre = sqrt(upre2);

	if (upre > 40000)
	{
		upre1 = upre / pow((upre / 40000.0), 0.66);
		ur = upre1 * ur / upre;
		utheta = upre1 * utheta / upre;
		uz = upre1 * uz / upre;
	}



	MPDT[i][j].vir = ur;
	MPDT[i][j].vitheta = utheta;
	MPDT[i][j].viz = uz;

	MPDT[i][j].ver = ur;
	MPDT[i][j].vetheta = utheta;
	MPDT[i][j].vez = uz;



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

	pos_r = 15;
	start_z = 0;
	end_z = 160;

	pos_z = 160;
	start_r = 0;
	end_r = 15;

	current_I = 0;

	//j = pos_r;

	//for (i = start_z; i < end_z; i++)
	//{
	//	current_I += (MPDT[i][j].ni * MPDT[i][j].vir - MPDT[i][j].ne * MPDT[i][j].ver) * QE * dr * 2 * PI * j * dz;
	//}
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

//void current_caulate()
//{
//	int i, j;
//	current_I = 0;
//	for (int k = 0; k < cath_num; k++)
//	{
//		i = cath_cell_z[k];
//		j = cath_cell_r[k];
//		if (cath_cell_dir[k] == LEFT || cath_cell_dir[k] == RIGHT)
//		{
//			if (j != 0)
//			{
//				current_I += (MPDT[i][j].ni * MPDT[i][j].viz - MPDT[i][j].ne * MPDT[i][j].vez) * QE * sqr(dr) * PI * (2 * j - 1);
//			}
//		}
//
//
//		if (cath_cell_dir[k] == DOWN || cath_cell_dir[k] == UP)
//		{
//			current_I += (MPDT[i][j].ni * MPDT[i][j].vir - MPDT[i][j].ne * MPDT[i][j].ver) * QE * dr * 2 * PI * j * dz;
//		}
//
//	}
//	printf("current = %lf\n", current_I);
//}

void current_control()
{
	pid.err = pid.set_current - (-current_I);
	pid.integral += pid.err;

	//if (pid.err > 50 || pid.err < -50)
	//{
	//	pid.ne_density = pid.Kp * pid.err + pid.Kd * (pid.err - pid.err_last);
	//}
	//else
	//{
		pid.ne_density = pid.Kp * pid.err + pid.Ki * pid.integral + pid.Kd * (pid.err - pid.err_last);
	//}

	pid.err_last = pid.err;
	inter_e_den += pid.ne_density;

	if (inter_e_den < 0)
	{
		inter_e_den = 0;
	}
	out_e_den = inter_e_den * cathod_cell / anode_cell;

}

