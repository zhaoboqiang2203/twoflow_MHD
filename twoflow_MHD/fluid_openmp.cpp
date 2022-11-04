#include "fluid.h"

void electron_flow_openmp()
{
	struct _U Qr;
	struct _U Qz;
	register int i, j, k;
	omp_set_num_threads(16);
#pragma omp	parallel for
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			//printf("i = %d, I am Thread %d\n", i, omp_get_thread_num());
			U[i][j] = cal_u(i, j);
			Fr[i][j] = cal_fr(U[i][j]);
			Fz[i][j] = cal_fz(U[i][j]);
			s[i][j] = cal_s(U[i][j]);

			U_bar[i][j].z = i;
			U_bar[i][j].r = j;
			U_bar2[i][j].z = i;
			U_bar2[i][j].r = j;
		}
	}

	cal_tau_omp();

#pragma omp	parallel for
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				for (k = 0; k < 13; k++)
				{
					U_bar[i][j].u[k] = U[i][j].u[k] - dt / dr * (Fr[i][j + 1].f[k] - Fr[i][j].f[k]) - dt / dz * (Fz[i + 1][j].f[k] - Fz[i][j].f[k]) + dt * s[i][j].f[k];
				}
#ifdef FLUID_DEBUG
				if (i == row && j == col)
				{
					printf("U[i][j].u1 = %.5e\n", U[i][j].u[0]);
					printf("Fr[i][j + 1].f1 = %.5e\n", Fr[i][j + 1].f[0]);
					printf("Fr[i][j].f1 = %.5e\n", Fr[i][j].f[0]);
					printf("Fz[i + 1][j].f1 = %.5e\n", Fz[i + 1][j].f[0]);
					printf("Fz[i][j].f1 = %.5e\n", Fz[i][j].f[0]);
					printf("s[i][j].f1 = %.5e\n", s[i][j].f[0]);
					printf("U_bar[i][j].u1 = %.5e\n", U_bar[i][j].u[0]);

					printf("U[i][j].u3 = %.5e\n", U[i][j].u[2]);
					printf("Fr[i][j + 1].f3 = %.5e\n", Fr[i][j + 1].f[2]);
					printf("Fr[i][j].f3 = %.5e\n", Fr[i][j].f[2]);
					printf("Fz[i + 1][j].f3 = %.5e\n", Fz[i + 1][j].f[2]);
					printf("Fz[i][j].f3 = %.5e\n", Fz[i][j].f[2]);
					printf("s[i][j].f3 = %.5e\n", s[i][j].f[2]);
					printf("U_bar[i][j].u3 = %.5e\n", U_bar[i][j].u[2]);
				}
#endif// FLUID_DEBUG
			}
			else if (btype[i][j] == 0)
			{
				for (k = 0; k < 13; k++)
				{
					U_bar[i][j].u[k] = 0;
				}
			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				for (k = 0; k < 13; k++)
				{
					//U_bar[i][j].u[k] = U[i][j].u[k] - dt / dr * (Fr[i][j + 1].f[k] - Fr[i][j].f[k]) - dt / dz * (Fz[i + 1][j].f[k] - Fz[i][j].f[k]) + dt * s[i][j].f[k];
					U_bar[i][j].u[k] = U[i][j].u[k] - dt / dz * (Fz[i + 1][j].f[k] - Fz[i][j].f[k]) + dt * s[i][j].f[k];
				}
			}
			else
			{
				for (k = 0; k < 13; k++)
				{
					U_bar[i][j].u[k] = U[i][j].u[k];
				}
			}

			Fr_bar[i][j] = cal_fr(U_bar[i][j]);
			Fz_bar[i][j] = cal_fz(U_bar[i][j]);
			s_bar[i][j] = cal_s(U_bar[i][j]);

		}
	}
	cal_tau_bar_omp();

	//矫正步
#pragma omp	parallel for
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				for (k = 0; k < 13; k++)
				{
					U_bar2[i][j].u[k] = 0.5 * (U[i][j].u[k] + U_bar[i][j].u[k] - dt / dr * (Fr_bar[i][j].f[k] - Fr_bar[i][j - 1].f[k]) - dt / dz * (Fz_bar[i][j].f[k] - Fz_bar[i - 1][j].f[k]) + dt * s_bar[i][j].f[k]);
				}
#ifdef FLUID_DEBUG
				if (i == row && j == col)
				{
					printf("Fr_bar[i][j].f[0] = %.5e\n", Fr_bar[i][j].f[0]);
					printf("Fr_bar[i][j - 1].f[0] = %.5e\n", Fr_bar[i][j - 1].f[0]);
					printf("Fz_bar[i][j].f[0] = %.5e\n", Fz_bar[i][j].f[0]);
					printf("Fz_bar[i - 1][j].f[0] = %.5e\n", Fz_bar[i - 1][j].f[0]);
					printf("s_bar[i][j].f[0] = %.5e\n", s_bar[i][j].f[0]);
					printf("U_bar2[i][j].u[0] = %.5e\n", U_bar2[i][j].u[0]);

					printf("Fr_bar[i][j].f3 = %.5e\n", Fr_bar[i][j].f[2]);
					printf("Fr_bar[i][j - 1].f3 = %.5e\n", Fr_bar[i][j - 1].f[2]);
					printf("Fz_bar[i][j].f3 = %.5e\n", Fz_bar[i][j].f[2]);
					printf("Fz_bar[i - 1][j].f3 = %.5e\n", Fz_bar[i - 1][j].f[2]);
					printf("s_bar[i][j].f3 = %.5e\n", s_bar[i][j].f[2]);
					printf("U_bar2[i][j].u3 = %.5e\n", U_bar2[i][j].u[2]);

				}
#endif// FLUID_DEBUG				

			}
			else if (btype[i][j] == 0)
			{
				for (k = 0; k < 13; k++)
				{
					U_bar2[i][j].u[k] = 0;
				}
			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				for (k = 0; k < 13; k++)
				{
					//U_bar2[i][j].u[k] = 0.5 * (U[i][j].u[k] + U_bar[i][j].u[k] - dt / dr * (Fr_bar[i][j + 1].f[k] - Fr_bar[i][j].f[k]) - dt / dz * (Fz_bar[i][j].f[k] - Fz_bar[i - 1][j].f[k]) + dt * s_bar[i][j].f[k]);
					U_bar2[i][j].u[k] = 0.5 * (U[i][j].u[k] + U_bar[i][j].u[k] - dt / dz * (Fz_bar[i][j].f[k] - Fz_bar[i - 1][j].f[k]) + dt * s_bar[i][j].f[k]);
				}
			}
			else
			{
				for (k = 0; k < 13; k++)
				{
					U_bar2[i][j].u[k] = U_bar[i][j].u[k];
				}
			}

		}
	}

	//计算				

#pragma omp	parallel for
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				Qz = arti_vis(MPDT[i + 1][j], MPDT[i][j], MPDT[i - 1][j]);
				Qr = arti_vis(MPDT[i][j + 1], MPDT[i][j], MPDT[i][j - 1]);

				U[i][j].u[0] = U_bar2[i][j].u[0] + Qr.u[0] / 4 * (MPDT[i][j + 1].ne - 2 * MPDT[i][j].ne + MPDT[i][j - 1].ne) + Qz.u[0] / 4 * (MPDT[i + 1][j].ne - 2 * MPDT[i][j].ne + MPDT[i - 1][j].ne);
				U[i][j].u[1] = U_bar2[i][j].u[1] + Qr.u[1] / 4 * (MPDT[i][j + 1].ni - 2 * MPDT[i][j].ni + MPDT[i][j - 1].ni) + Qz.u[1] / 4 * (MPDT[i + 1][j].ni - 2 * MPDT[i][j].ni + MPDT[i - 1][j].ni);
				U[i][j].u[2] = U_bar2[i][j].u[2] + Qr.u[2] / 4 * (MPDT[i][j + 1].ver - 2 * MPDT[i][j].ver + MPDT[i][j - 1].ver) + Qz.u[2] / 4 * (MPDT[i + 1][j].ver - 2 * MPDT[i][j].ver + MPDT[i - 1][j].ver);
				U[i][j].u[3] = U_bar2[i][j].u[3] + Qr.u[3] / 4 * (MPDT[i][j + 1].vetheta - 2 * MPDT[i][j].vetheta + MPDT[i][j - 1].vetheta) + Qz.u[3] / 4 * (MPDT[i + 1][j].vetheta - 2 * MPDT[i][j].vetheta + MPDT[i - 1][j].vetheta);
				U[i][j].u[4] = U_bar2[i][j].u[4] + Qr.u[4] / 4 * (MPDT[i][j + 1].vez - 2 * MPDT[i][j].vez + MPDT[i][j - 1].vez) + Qz.u[4] / 4 * (MPDT[i + 1][j].vez - 2 * MPDT[i][j].vez + MPDT[i - 1][j].vez);
				U[i][j].u[5] = U_bar2[i][j].u[5] + Qr.u[5] / 4 * (MPDT[i][j + 1].vir - 2 * MPDT[i][j].vir + MPDT[i][j - 1].vir) + Qz.u[5] / 4 * (MPDT[i + 1][j].vir - 2 * MPDT[i][j].vir + MPDT[i - 1][j].vir);
				U[i][j].u[6] = U_bar2[i][j].u[6] + Qr.u[6] / 4 * (MPDT[i][j + 1].vitheta - 2 * MPDT[i][j].vitheta + MPDT[i][j - 1].vitheta) + Qz.u[6] / 4 * (MPDT[i + 1][j].vitheta - 2 * MPDT[i][j].vitheta + MPDT[i - 1][j].vitheta);
				U[i][j].u[7] = U_bar2[i][j].u[7] + Qr.u[7] / 4 * (MPDT[i][j + 1].viz - 2 * MPDT[i][j].viz + MPDT[i][j - 1].viz) + Qz.u[7] / 4 * (MPDT[i + 1][j].viz - 2 * MPDT[i][j].viz + MPDT[i - 1][j].viz);
				U[i][j].u[8] = U_bar2[i][j].u[8] + Qr.u[8] / 4 * (MPDT[i][j + 1].br - 2 * MPDT[i][j].br + MPDT[i][j - 1].br) + Qz.u[8] / 4 * (MPDT[i + 1][j].br - 2 * MPDT[i][j].br + MPDT[i - 1][j].br);
				U[i][j].u[9] = U_bar2[i][j].u[9] + Qr.u[9] / 4 * (MPDT[i][j + 1].btheta - 2 * MPDT[i][j].btheta + MPDT[i][j - 1].btheta) + Qz.u[9] / 4 * (MPDT[i + 1][j].btheta - 2 * MPDT[i][j].btheta + MPDT[i - 1][j].btheta);
				U[i][j].u[10] = U_bar2[i][j].u[10] + Qr.u[10] / 4 * (MPDT[i][j + 1].bz - 2 * MPDT[i][j].bz + MPDT[i][j - 1].bz) + Qz.u[10] / 4 * (MPDT[i + 1][j].bz - 2 * MPDT[i][j].bz + MPDT[i - 1][j].bz);
				U[i][j].u[11] = U_bar2[i][j].u[11] + Qr.u[11] / 4 * (MPDT[i][j + 1].ee - 2 * MPDT[i][j].ee + MPDT[i][j - 1].ee) + Qz.u[11] / 4 * (MPDT[i + 1][j].ee - 2 * MPDT[i][j].ee + MPDT[i - 1][j].ee);
				U[i][j].u[12] = U_bar2[i][j].u[12] + Qr.u[12] / 4 * (MPDT[i][j + 1].ei - 2 * MPDT[i][j].ei + MPDT[i][j - 1].ei) + Qz.u[12] / 4 * (MPDT[i + 1][j].ei - 2 * MPDT[i][j].ei + MPDT[i - 1][j].ei);


#ifdef FLUID_DEBUG
				if (i == row && j == col)
				{
					//cout << "final U[i][j].u1 = " << U[i][j].u1 << endl;
					printf("final U[i][j].u1 = %.5e\n", U[i][j].u[0]);
					printf("U_bar2[i][j + 1].u1 = %.5e\n", U_bar2[i][j + 1].u[0]);
					printf("U_bar2[i][j - 1].u1 = %.5e\n", U_bar2[i][j - 1].u[0]);
					printf("U_bar2[i + 1][j].u1 = %.5e\n", U_bar2[i + 1][j].u[0]);
					printf("U_bar2[i - 1][j].u1 = %.5e\n", U_bar2[i - 1][j].u[0]);

					printf("final U[i][j].u3 = %.5e\n", U[i][j].u[2]);
					printf("U_bar2[i][j + 1].u3 = %.5e\n", U_bar2[i][j + 1].u[2]);
					printf("U_bar2[i][j - 1].u3 = %.5e\n", U_bar2[i][j - 1].u[2]);
					printf("U_bar2[i + 1][j].u3 = %.5e\n", U_bar2[i + 1][j].u[2]);
					printf("U_bar2[i - 1][j].u3 = %.5e\n", U_bar2[i - 1][j].u[2]);

				}
#endif// FLUID_DEBUG	

			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				Qz = arti_vis(MPDT[i + 1][j], MPDT[i][j], MPDT[i - 1][j]);

				U[i][j].u[0] = U_bar2[i][j].u[0] + Qz.u[0] / 2 * (MPDT[i + 1][j].ne - 2 * MPDT[i][j].ne + MPDT[i - 1][j].ne);
				U[i][j].u[1] = U_bar2[i][j].u[1] + Qz.u[1] / 2 * (MPDT[i + 1][j].ni - 2 * MPDT[i][j].ni + MPDT[i - 1][j].ni);
				U[i][j].u[2] = U_bar2[i][j].u[2] + Qz.u[2] / 2 * (MPDT[i + 1][j].ver - 2 * MPDT[i][j].ver + MPDT[i - 1][j].ver);
				U[i][j].u[3] = U_bar2[i][j].u[3] + Qz.u[3] / 2 * (MPDT[i + 1][j].vetheta - 2 * MPDT[i][j].vetheta + MPDT[i - 1][j].vetheta);
				U[i][j].u[4] = U_bar2[i][j].u[4] + Qz.u[4] / 2 * (MPDT[i + 1][j].vez - 2 * MPDT[i][j].vez + MPDT[i - 1][j].vez);
				U[i][j].u[5] = U_bar2[i][j].u[5] + Qz.u[5] / 2 * (MPDT[i + 1][j].vir - 2 * MPDT[i][j].vir + MPDT[i - 1][j].vir);
				U[i][j].u[6] = U_bar2[i][j].u[6] + Qz.u[6] / 2 * (MPDT[i + 1][j].vitheta - 2 * MPDT[i][j].vitheta + MPDT[i - 1][j].vitheta);
				U[i][j].u[7] = U_bar2[i][j].u[7] + Qz.u[7] / 2 * (MPDT[i + 1][j].viz - 2 * MPDT[i][j].viz + MPDT[i - 1][j].viz);
				U[i][j].u[8] = U_bar2[i][j].u[8] + Qz.u[8] / 2 * (MPDT[i + 1][j].br - 2 * MPDT[i][j].br + MPDT[i - 1][j].br);
				U[i][j].u[9] = U_bar2[i][j].u[9] + Qz.u[9] / 2 * (MPDT[i + 1][j].btheta - 2 * MPDT[i][j].btheta + MPDT[i - 1][j].btheta);
				U[i][j].u[10] = U_bar2[i][j].u[10] + Qz.u[10] / 2 * (MPDT[i + 1][j].bz - 2 * MPDT[i][j].bz + MPDT[i - 1][j].bz);
				U[i][j].u[11] = U_bar2[i][j].u[11] + Qz.u[11] / 2 * (MPDT[i + 1][j].ee - 2 * MPDT[i][j].ee + MPDT[i - 1][j].ee);
				U[i][j].u[12] = U_bar2[i][j].u[12] + Qz.u[12] / 2 * (MPDT[i + 1][j].ei - 2 * MPDT[i][j].ei + MPDT[i - 1][j].ei);

			}
			else //if (btype[i][j] == 0)
			{
				for (k = 0; k < 13; k++)
				{
					U[i][j].u[k] = U_bar2[i][j].u[k];
				}
			}
		}
	}
	//更新边界信息



#pragma omp	parallel for
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{

			MPDT[i][j].ne = U[i][j].u[0];
			MPDT[i][j].ni = U[i][j].u[1];
			MPDT[i][j].ver = U[i][j].u[2];
			MPDT[i][j].vetheta = U[i][j].u[3];
			MPDT[i][j].vez = U[i][j].u[4];


			MPDT[i][j].vir = U[i][j].u[5];
			MPDT[i][j].vitheta = U[i][j].u[6];
			MPDT[i][j].viz = U[i][j].u[7];


			MPDT[i][j].br = U[i][j].u[8];
			MPDT[i][j].btheta = U[i][j].u[9];
			MPDT[i][j].bz = U[i][j].u[10];

			MPDT[i][j].ee = U[i][j].u[11];
			MPDT[i][j].ei = U[i][j].u[12];

			MPDT[i][j].pe = (gamma - 1) * (MPDT[i][j].ne * ME * MPDT[i][j].ee - 0.5 * MPDT[i][j].ne * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez));
			MPDT[i][j].pi = (gamma - 1) * (MPDT[i][j].ni * MI * MPDT[i][j].ei - 0.5 * MPDT[i][j].ni * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz));

		}

	}

	return;
}


void cal_tau_omp()
{
	double mu_vis = 0;
	register int i = 0, j = 0, k, l;
	int nt = nz * nr;

#pragma omp	parallel for
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			mu_vis = dy_vis(i, j);
			MPDT[i][j].sigma_Q = mu_vis;
			tau_ni[i][j] = U[i][j].u[1];
			if (btype[i][j] == 1)
			{
				taurr[i][j] = 2 * mu_vis * (U[i][j + 1].u[5] - U[i][j - 1].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U[i][j + 1].u[6] - U[i][j - 1].u[6]) / (2 * dr) - U[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U[i][j + 1].u[7] - U[i][j - 1].u[7]) / (2 * dr) + (U[i + 1][j].u[5] - U[i - 1][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U[i + 1][j].u[6] - U[i - 1][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U[i + 1][j].u[7] - U[i - 1][j].u[7]) / (2 * dz);
				pniz[i][j] = (U[i + 1][j].u[1] - U[i - 1][j].u[1]) / (2 * dz);
				pnir[i][j] = (U[i][j + 1].u[1] - U[i][j - 1].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == 0)
			{
				taurr[i][j] = 0;
				taurtheta[i][j] = 0;
				taurz[i][j] = 0;
				tautheta2[i][j] = 0;
				tauthetaz[i][j] = 0;
				tauzz[i][j] = 0;
			}
			else if (btype[i][j] == LEFT)
			{
				taurr[i][j] = 2 * mu_vis * (U[i][j + 1].u[5] - U[i][j - 1].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U[i][j + 1].u[6] - U[i][j - 1].u[6]) / (2 * dr) - U[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U[i][j + 1].u[7] - U[i][j - 1].u[7]) / (2 * dr) + (-3 * U[i][j].u[5] + 3 * U[i + 1][j].u[5] - U[i + 2][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((-3 * U[i][j].u[6] + 3 * U[i + 1][j].u[6] - U[i + 2][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (-3 * U[i][j].u[7] + 3 * U[i + 1][j].u[7] - U[i + 2][j].u[7]) / (2 * dz);
				pniz[i][j] = (-3 * U[i][j].u[1] + 3 * U[i + 1][j].u[1] - U[i + 2][j].u[1]) / (2 * dz);
				pnir[i][j] = (U[i][j + 1].u[1] - U[i][j - 1].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == RIGHT)
			{
				taurr[i][j] = 2 * mu_vis * (U[i][j + 1].u[5] - U[i][j - 1].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U[i][j + 1].u[6] - U[i][j - 1].u[6]) / (2 * dr) - U[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U[i][j + 1].u[7] - U[i][j - 1].u[7]) / (2 * dr) + (U[i - 2][j].u[5] - 4 * U[i - 1][j].u[5] + 3 * U[i][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U[i - 2][j].u[6] - 4 * U[i - 1][j].u[6] + 3 * U[i][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U[i - 2][j].u[7] - 4 * U[i - 1][j].u[7] + 3 * U[i][j].u[7]) / (2 * dz);
				pniz[i][j] = (U[i - 2][j].u[1] - 4 * U[i - 1][j].u[1] + 3 * U[i][j].u[1]) / (2 * dz);
				pnir[i][j] = (U[i][j + 1].u[1] - U[i][j - 1].u[1]) / (2 * dr);

			}
			else if (btype[i][j] == UP)
			{
				taurr[i][j] = 2 * mu_vis * (U[i][j - 2].u[5] - 4 * U[i][j - 1].u[5] + 3 * U[i][j].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U[i][j - 2].u[6] - 4 * U[i][j - 1].u[6] + 3 * U[i][j].u[6]) / (2 * dr) - U[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U[i][j - 2].u[7] - 4 * U[i][j - 1].u[7] + 3 * U[i][j].u[7]) / (2 * dr) + (U[i + 1][j].u[5] - U[i - 1][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U[i + 1][j].u[6] - U[i - 1][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U[i + 1][j].u[7] - U[i][j + 1].u[7]) / (2 * dz);
				pniz[i][j] = (U[i + 1][j].u[1] - U[i - 1][j].u[1]) / (2 * dz);
				pnir[i][j] = (U[i][j - 2].u[1] - 4 * U[i][j - 1].u[1] + 3 * U[i][j].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == DOWN)
			{
				taurr[i][j] = 2 * mu_vis * (-3 * U[i][j].u[5] + 3 * U[i][j + 1].u[5] - U[i][j + 2].u[5]) / (2 * dr);
				if (j == 0)
				{
					taurtheta[i][j] = mu_vis * (-3 * U[i][j].u[6] + 3 * U[i][j + 1].u[6] - U[i][j + 2].u[6]) / (2 * dr);
				}
				else
				{
					taurtheta[i][j] = mu_vis * ((-3 * U[i][j].u[6] + 3 * U[i][j + 1].u[6] - U[i][j + 2].u[6]) / (2 * dr) - U[i][j].u[6] / (j * dr));
				}
				taurz[i][j] = mu_vis * ((-3 * U[i][j].u[7] + 3 * U[i][j + 1].u[7] - U[i][j + 2].u[7]) / (2 * dr) + (U[i + 1][j].u[5] - U[i - 1][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U[i + 1][j].u[6] - U[i - 1][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U[i + 1][j].u[7] - U[i][j + 1].u[7]) / (2 * dz);
				pniz[i][j] = (U[i + 1][j].u[1] - U[i - 1][j].u[1]) / (2 * dz);
				pnir[i][j] = (-3 * U[i][j].u[1] + 3 * U[i][j + 1].u[1] - U[i][j + 2].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == (LEFT + UP))
			{
				taurr[i][j] = 2 * mu_vis * (U[i][j - 2].u[5] - 4 * U[i][j - 1].u[5] + 3 * U[i][j].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U[i][j - 2].u[6] - 4 * U[i][j - 1].u[6] + 3 * U[i][j].u[6]) / (2 * dr) - U[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U[i][j - 2].u[7] - 4 * U[i][j - 1].u[7] + 3 * U[i][j].u[7]) / (2 * dr) + (-3 * U[i][j].u[5] + 3 * U[i + 1][j].u[5] - U[i + 2][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((-3 * U[i][j].u[6] + 3 * U[i + 1][j].u[6] - U[i + 2][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (-3 * U[i][j].u[7] + 3 * U[i + 1][j].u[7] - U[i + 2][j].u[7]) / (2 * dz);
				pniz[i][j] = (-3 * U[i][j].u[1] + 3 * U[i + 1][j].u[1] - U[i + 2][j].u[1]) / (2 * dz);
				pnir[i][j] = (U[i][j - 2].u[1] - 4 * U[i][j - 1].u[1] + 3 * U[i][j].u[1]) / (2 * dr);

			}
			else if (btype[i][j] == (LEFT + DOWN))
			{
				taurr[i][j] = 2 * mu_vis * (-3 * U[i][j].u[5] + 3 * U[i][j + 1].u[5] - U[i][j + 2].u[5]) / (2 * dr);
				if (j == 0)
				{
					taurtheta[i][j] = mu_vis * (-3 * U[i][j].u[6] + 3 * U[i][j + 1].u[6] - U[i][j + 2].u[6]) / (2 * dr);
				}
				else
				{
					taurtheta[i][j] = mu_vis * ((-3 * U[i][j].u[6] + 3 * U[i][j + 1].u[6] - U[i][j + 2].u[6]) / (2 * dr) - U[i][j].u[6] / (j * dr));
				}
				taurz[i][j] = mu_vis * ((-3 * U[i][j].u[7] + 3 * U[i][j + 1].u[7] - U[i][j + 2].u[7]) / (2 * dr) + (-3 * U[i][j].u[5] + 3 * U[i + 1][j].u[5] - U[i + 2][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((-3 * U[i][j].u[6] + 3 * U[i + 1][j].u[6] - U[i + 2][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (-3 * U[i][j].u[7] + 3 * U[i + 1][j].u[7] - U[i + 2][j].u[7]) / (2 * dz);
				pniz[i][j] = (-3 * U[i][j].u[1] + 3 * U[i + 1][j].u[1] - U[i + 2][j].u[1]) / (2 * dz);
				pnir[i][j] = (-3 * U[i][j].u[1] + 3 * U[i][j + 1].u[1] - U[i][j + 2].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == (RIGHT + DOWN))
			{
				taurr[i][j] = 2 * mu_vis * (-3 * U[i][j].u[5] + 3 * U[i][j + 1].u[5] - U[i][j + 2].u[5]) / (2 * dr);
				if (j == 0)
				{
					taurtheta[i][j] = mu_vis * (-3 * U[i][j].u[6] + 3 * U[i][j + 1].u[6] - U[i][j + 2].u[6]) / (2 * dr);
				}
				else
				{
					taurtheta[i][j] = mu_vis * ((-3 * U[i][j].u[6] + 3 * U[i][j + 1].u[6] - U[i][j + 2].u[6]) / (2 * dr) - U[i][j].u[6] / (j * dr));
				}
				taurz[i][j] = mu_vis * ((-3 * U[i][j].u[7] + 3 * U[i][j + 1].u[7] - U[i][j + 2].u[7]) / (2 * dr) + (U[i - 2][j].u[5] - 4 * U[i - 1][j].u[5] + 3 * U[i][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U[i - 2][j].u[6] - 4 * U[i - 1][j].u[6] + 3 * U[i][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U[i - 2][j].u[7] - 4 * U[i - 1][j].u[7] + 3 * U[i][j].u[7]) / (2 * dz);
				pniz[i][j] = (U[i - 2][j].u[1] - 4 * U[i - 1][j].u[1] + 3 * U[i][j].u[1]) / (2 * dz);
				pnir[i][j] = (-3 * U[i][j].u[1] + 3 * U[i][j + 1].u[1] - U[i][j + 2].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == (RIGHT + UP))
			{
				taurr[i][j] = 2 * mu_vis * (U[i][j - 2].u[5] - 4 * U[i][j - 1].u[5] + 3 * U[i][j].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U[i][j - 2].u[6] - 4 * U[i][j - 1].u[6] + 3 * U[i][j].u[6]) / (2 * dr) - U[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U[i][j - 2].u[7] - 4 * U[i][j - 1].u[7] + 3 * U[i][j].u[7]) / (2 * dr) + (U[i - 2][j].u[5] - 4 * U[i - 1][j].u[5] + 3 * U[i][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U[i - 2][j].u[6] - 4 * U[i - 1][j].u[6] + 3 * U[i][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U[i - 2][j].u[7] - 4 * U[i - 1][j].u[7] + 3 * U[i][j].u[7]) / (2 * dz);
				pniz[i][j] = (U[i - 2][j].u[1] - 4 * U[i - 1][j].u[1] + 3 * U[i][j].u[1]) / (2 * dz);
				pnir[i][j] = (U[i][j - 2].u[1] - 4 * U[i][j - 1].u[1] + 3 * U[i][j].u[1]) / (2 * dr);
			}
			else
			{

			}

		}
	}
#pragma omp	parallel for
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			for (k = 0; k < 13; k++)
			{
				tau_vis[i][j].f[k] = 0;
			}

			if (btype[i][j] == 1)
			{
				if (tau_ni[i][j] != 0)
				{
					tau_vis[i][j].f[5] = (taurr[i][j + 1] - taurr[i][j - 1]) / (2 * dr) + taurr[i][j] * pnir[i][j] / tau_ni[i][j] + taurr[i][j] / (dr * j) + taurz[i][j] * pniz[i][j] / tau_ni[i][j] + (taurz[i + 1][j] - taurz[i - 1][j]) / (2 * dz);
					tau_vis[i][j].f[6] = (taurtheta[i][j + 1] - taurtheta[i][j - 1]) / (2 * dr) + taurtheta[i][j] * pnir[i][j] / tau_ni[i][j] + taurtheta[i][j] / (dr * j) + tauthetaz[i][j] * pniz[i][j] / tau_ni[i][j] + (tauthetaz[i + 1][j] - tauthetaz[i - 1][j]) / (2 * dz);
					tau_vis[i][j].f[7] = (taurz[i][j + 1] - taurz[i][j - 1]) / (2 * dr) + taurz[i][j] * pnir[i][j] / tau_ni[i][j] + taurz[i][j] / (dr * j) + tauzz[i][j] * pniz[i][j] / tau_ni[i][j] + (tauzz[i + 1][j] - tauzz[i - 1][j]) / (2 * dz);
					tau_vis[i][j].f[12] = tau_vis[i][j].f[5] * U[i][j].u[5] + tau_vis[i][j].f[6] * U[i][j].u[6] + tau_vis[i][j].f[7] * U[i][j].u[7];


				}

			}
			else if (btype[i][j] == 0)
			{
				tau_vis[i][j].f[5] = 0;
				tau_vis[i][j].f[6] = 0;
				tau_vis[i][j].f[7] = 0;
			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				tau_vis[i][j].f[5] = taurr[i][j] / tau_ni[i][j] * pnir[i][j] + taurz[i][j] / tau_ni[i][j] * pniz[i][j] + (taurz[i + 1][j] - taurz[i - 1][j]) / (2 * dz);
				tau_vis[i][j].f[6] = taurtheta[i][j] / tau_ni[i][j] * pnir[i][j] + tauthetaz[i][j] / tau_ni[i][j] * pniz[i][j] + (tauthetaz[i + 1][j] - tauthetaz[i - 1][j]) / (2 * dz);
				tau_vis[i][j].f[7] = taurz[i][j] / tau_ni[i][j] * pnir[i][j] + tauzz[i][j] / tau_ni[i][j] * pniz[i][j] + (tauzz[i + 1][j] - tauzz[i - 1][j]) / (2 * dz);
				tau_vis[i][j].f[12] = tau_vis[i][j].f[5] * U[i][j].u[5] + tau_vis[i][j].f[6] * U[i][j].u[6] + tau_vis[i][j].f[7] * U[i][j].u[7];

			}
			else
			{
				tau_vis[i][j].f[5] = 0;
				tau_vis[i][j].f[6] = 0;
				tau_vis[i][j].f[7] = 0;
			}

			//MPDT[i][j].tau_vise = ;
			//MPDT[i][j].tau_visr;
			//MPDT[i][j].tau_visz;
			//MPDT[i][j].tau_vistheta;

			//MPDT[i][j].neq = tau_vis[i][j].f[12];
			//MPDT[i][j].vnqr = tau_vis[i][j].f[5];
			//MPDT[i][j].vnqz = tau_vis[i][j].f[7];
			//MPDT[i][j].vnqtheta = tau_vis[i][j].f[6];
		}
	}
#pragma omp	parallel for
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			for (k = 0; k < 13; k++)
			{
				s[i][j].f[k] += tau_vis[i][j].f[k];
			}
		}
	}
}

void cal_tau_bar_omp()
{
	double mu_vis = 0;
	register int i = 0, j = 0, k = 0, l;
	int nt = nz * nr;


#pragma omp	parallel for
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			mu_vis = dy_vis(i, j);
			tau_ni[i][j] = U_bar[i][j].u[1];
			if (btype[i][j] == 1)
			{
				taurr[i][j] = 2 * mu_vis * (U_bar[i][j + 1].u[5] - U_bar[i][j - 1].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U_bar[i][j + 1].u[6] - U_bar[i][j - 1].u[6]) / (2 * dr) - U_bar[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U_bar[i][j + 1].u[7] - U_bar[i][j - 1].u[7]) / (2 * dr) + (U_bar[i + 1][j].u[5] - U_bar[i - 1][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U_bar[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U_bar[i + 1][j].u[6] - U_bar[i - 1][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U_bar[i + 1][j].u[7] - U_bar[i - 1][j].u[7]) / (2 * dz);
				pniz[i][j] = (U_bar[i + 1][j].u[1] - U_bar[i - 1][j].u[1]) / (2 * dz);
				pnir[i][j] = (U_bar[i][j + 1].u[1] - U_bar[i][j - 1].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == 0)
			{
				taurr[i][j] = 0;
				taurtheta[i][j] = 0;
				taurz[i][j] = 0;
				tautheta2[i][j] = 0;
				tauthetaz[i][j] = 0;
				tauzz[i][j] = 0;
				pniz[i][j] = 0;
				pnir[i][j] = 0;
			}
			else if (btype[i][j] == LEFT)
			{
				taurr[i][j] = 2 * mu_vis * (U_bar[i][j + 1].u[5] - U_bar[i][j - 1].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U_bar[i][j + 1].u[6] - U_bar[i][j - 1].u[6]) / (2 * dr) - U_bar[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U_bar[i][j + 1].u[7] - U_bar[i][j - 1].u[7]) / (2 * dr) + (-3 * U_bar[i][j].u[5] + 3 * U_bar[i + 1][j].u[5] - U_bar[i + 2][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U_bar[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((-3 * U_bar[i][j].u[6] + 3 * U_bar[i + 1][j].u[6] - U_bar[i + 2][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (-3 * U_bar[i][j].u[7] + 3 * U_bar[i + 1][j].u[7] - U_bar[i + 2][j].u[7]) / (2 * dz);
				pniz[i][j] = (-3 * U_bar[i][j].u[1] + 3 * U_bar[i + 1][j].u[1] - U_bar[i + 2][j].u[1]) / (2 * dz);
				pnir[i][j] = (U_bar[i][j + 1].u[1] - U_bar[i][j - 1].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == RIGHT)
			{
				taurr[i][j] = 2 * mu_vis * (U_bar[i][j + 1].u[5] - U_bar[i][j - 1].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U_bar[i][j + 1].u[6] - U_bar[i][j - 1].u[6]) / (2 * dr) - U_bar[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U_bar[i][j + 1].u[7] - U_bar[i][j - 1].u[7]) / (2 * dr) + (U_bar[i - 2][j].u[5] - 4 * U_bar[i - 1][j].u[5] + 3 * U_bar[i][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U_bar[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U_bar[i - 2][j].u[6] - 4 * U_bar[i - 1][j].u[6] + 3 * U_bar[i][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U_bar[i - 2][j].u[7] - 4 * U_bar[i - 1][j].u[7] + 3 * U_bar[i][j].u[7]) / (2 * dz);
				pniz[i][j] = (U_bar[i - 2][j].u[1] - 4 * U_bar[i - 1][j].u[1] + 3 * U_bar[i][j].u[1]) / (2 * dz);
				pnir[i][j] = (U_bar[i][j + 1].u[1] - U_bar[i][j - 1].u[1]) / (2 * dr);

			}
			else if (btype[i][j] == UP)
			{
				taurr[i][j] = 2 * mu_vis * (U_bar[i][j - 2].u[5] - 4 * U_bar[i][j - 1].u[5] + 3 * U_bar[i][j].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U_bar[i][j - 2].u[6] - 4 * U_bar[i][j - 1].u[6] + 3 * U_bar[i][j].u[6]) / (2 * dr) - U_bar[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U_bar[i][j - 2].u[7] - 4 * U_bar[i][j - 1].u[7] + 3 * U_bar[i][j].u[7]) / (2 * dr) + (U_bar[i + 1][j].u[5] - U_bar[i - 1][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U_bar[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U_bar[i + 1][j].u[6] - U_bar[i - 1][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U_bar[i + 1][j].u[7] - U_bar[i][j + 1].u[7]) / (2 * dz);
				pniz[i][j] = (U_bar[i + 1][j].u[1] - U_bar[i - 1][j].u[1]) / (2 * dz);
				pnir[i][j] = (U_bar[i][j - 2].u[1] - 4 * U_bar[i][j - 1].u[1] + 3 * U_bar[i][j].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == DOWN)
			{
				taurr[i][j] = 2 * mu_vis * (-3 * U_bar[i][j].u[5] + 3 * U_bar[i][j + 1].u[5] - U_bar[i][j + 2].u[5]) / (2 * dr);
				if (j == 0)
				{
					taurtheta[i][j] = mu_vis * (-3 * U_bar[i][j].u[6] + 3 * U_bar[i][j + 1].u[6] - U_bar[i][j + 2].u[6]) / (2 * dr);
				}
				else
				{
					taurtheta[i][j] = mu_vis * ((-3 * U_bar[i][j].u[6] + 3 * U_bar[i][j + 1].u[6] - U_bar[i][j + 2].u[6]) / (2 * dr) - U_bar[i][j].u[6] / (j * dr));
				}
				taurz[i][j] = mu_vis * ((-3 * U_bar[i][j].u[7] + 3 * U_bar[i][j + 1].u[7] - U_bar[i][j + 2].u[7]) / (2 * dr) + (U_bar[i + 1][j].u[5] - U_bar[i - 1][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U_bar[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U_bar[i + 1][j].u[6] - U_bar[i - 1][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U_bar[i + 1][j].u[7] - U_bar[i][j + 1].u[7]) / (2 * dz);
				pniz[i][j] = (U_bar[i + 1][j].u[1] - U_bar[i - 1][j].u[1]) / (2 * dz);
				pnir[i][j] = (-3 * U_bar[i][j].u[1] + 3 * U_bar[i][j + 1].u[1] - U_bar[i][j + 2].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == (LEFT + UP))
			{
				taurr[i][j] = 2 * mu_vis * (U_bar[i][j - 2].u[5] - 4 * U_bar[i][j - 1].u[5] + 3 * U_bar[i][j].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U_bar[i][j - 2].u[6] - 4 * U_bar[i][j - 1].u[6] + 3 * U_bar[i][j].u[6]) / (2 * dr) - U_bar[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U_bar[i][j - 2].u[7] - 4 * U_bar[i][j - 1].u[7] + 3 * U_bar[i][j].u[7]) / (2 * dr) + (-3 * U_bar[i][j].u[5] + 3 * U_bar[i + 1][j].u[5] - U_bar[i + 2][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U_bar[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((-3 * U_bar[i][j].u[6] + 3 * U_bar[i + 1][j].u[6] - U_bar[i + 2][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (-3 * U_bar[i][j].u[7] + 3 * U_bar[i + 1][j].u[7] - U_bar[i + 2][j].u[7]) / (2 * dz);
				pniz[i][j] = (-3 * U_bar[i][j].u[1] + 3 * U_bar[i + 1][j].u[1] - U_bar[i + 2][j].u[1]) / (2 * dz);
				pnir[i][j] = (U_bar[i][j - 2].u[1] - 4 * U_bar[i][j - 1].u[1] + 3 * U_bar[i][j].u[1]) / (2 * dr);

			}
			else if (btype[i][j] == (LEFT + DOWN))
			{
				taurr[i][j] = 2 * mu_vis * (-3 * U_bar[i][j].u[5] + 3 * U_bar[i][j + 1].u[5] - U_bar[i][j + 2].u[5]) / (2 * dr);
				if (j == 0)
				{
					taurtheta[i][j] = mu_vis * (-3 * U_bar[i][j].u[6] + 3 * U_bar[i][j + 1].u[6] - U_bar[i][j + 2].u[6]) / (2 * dr);
				}
				else
				{
					taurtheta[i][j] = mu_vis * ((-3 * U_bar[i][j].u[6] + 3 * U_bar[i][j + 1].u[6] - U_bar[i][j + 2].u[6]) / (2 * dr) - U_bar[i][j].u[6] / (j * dr));
				}
				taurz[i][j] = mu_vis * ((-3 * U_bar[i][j].u[7] + 3 * U_bar[i][j + 1].u[7] - U_bar[i][j + 2].u[7]) / (2 * dr) + (-3 * U_bar[i][j].u[5] + 3 * U_bar[i + 1][j].u[5] - U_bar[i + 2][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U_bar[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((-3 * U_bar[i][j].u[6] + 3 * U_bar[i + 1][j].u[6] - U_bar[i + 2][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (-3 * U_bar[i][j].u[7] + 3 * U_bar[i + 1][j].u[7] - U_bar[i + 2][j].u[7]) / (2 * dz);
				pniz[i][j] = (-3 * U_bar[i][j].u[1] + 3 * U_bar[i + 1][j].u[1] - U_bar[i + 2][j].u[1]) / (2 * dz);
				pnir[i][j] = (-3 * U_bar[i][j].u[1] + 3 * U_bar[i][j + 1].u[1] - U_bar[i][j + 2].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == (RIGHT + DOWN))
			{
				taurr[i][j] = 2 * mu_vis * (-3 * U_bar[i][j].u[5] + 3 * U_bar[i][j + 1].u[5] - U_bar[i][j + 2].u[5]) / (2 * dr);
				if (j == 0)
				{
					taurtheta[i][j] = mu_vis * (-3 * U_bar[i][j].u[6] + 3 * U_bar[i][j + 1].u[6] - U_bar[i][j + 2].u[6]) / (2 * dr);
				}
				else
				{
					taurtheta[i][j] = mu_vis * ((-3 * U_bar[i][j].u[6] + 3 * U_bar[i][j + 1].u[6] - U_bar[i][j + 2].u[6]) / (2 * dr) - U_bar[i][j].u[6] / (j * dr));
				}
				taurz[i][j] = mu_vis * ((-3 * U_bar[i][j].u[7] + 3 * U_bar[i][j + 1].u[7] - U_bar[i][j + 2].u[7]) / (2 * dr) + (U_bar[i - 2][j].u[5] - 4 * U_bar[i - 1][j].u[5] + 3 * U_bar[i][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U_bar[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U_bar[i - 2][j].u[6] - 4 * U_bar[i - 1][j].u[6] + 3 * U_bar[i][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U_bar[i - 2][j].u[7] - 4 * U_bar[i - 1][j].u[7] + 3 * U_bar[i][j].u[7]) / (2 * dz);
				pniz[i][j] = (U_bar[i - 2][j].u[1] - 4 * U_bar[i - 1][j].u[1] + 3 * U_bar[i][j].u[1]) / (2 * dz);
				pnir[i][j] = (-3 * U_bar[i][j].u[1] + 3 * U_bar[i][j + 1].u[1] - U_bar[i][j + 2].u[1]) / (2 * dr);
			}
			else if (btype[i][j] == (RIGHT + UP))
			{
				taurr[i][j] = 2 * mu_vis * (U_bar[i][j - 2].u[5] - 4 * U_bar[i][j - 1].u[5] + 3 * U_bar[i][j].u[5]) / (2 * dr);
				taurtheta[i][j] = mu_vis * ((U_bar[i][j - 2].u[6] - 4 * U_bar[i][j - 1].u[6] + 3 * U_bar[i][j].u[6]) / (2 * dr) - U_bar[i][j].u[6] / (j * dr));
				taurz[i][j] = mu_vis * ((U_bar[i][j - 2].u[7] - 4 * U_bar[i][j - 1].u[7] + 3 * U_bar[i][j].u[7]) / (2 * dr) + (U_bar[i - 2][j].u[5] - 4 * U_bar[i - 1][j].u[5] + 3 * U_bar[i][j].u[5]) / (2 * dz));
				tautheta2[i][j] = 2 * mu_vis * U_bar[i][j].u[6] / (j * dr);
				tauthetaz[i][j] = mu_vis * ((U_bar[i - 2][j].u[6] - 4 * U_bar[i - 1][j].u[6] + 3 * U_bar[i][j].u[6]) / (2 * dz));
				tauzz[i][j] = 2 * mu_vis * (U_bar[i - 2][j].u[7] - 4 * U_bar[i - 1][j].u[7] + 3 * U_bar[i][j].u[7]) / (2 * dz);
				pniz[i][j] = (U_bar[i - 2][j].u[1] - 4 * U_bar[i - 1][j].u[1] + 3 * U_bar[i][j].u[1]) / (2 * dz);
				pnir[i][j] = (U_bar[i][j - 2].u[1] - 4 * U_bar[i][j - 1].u[1] + 3 * U_bar[i][j].u[1]) / (2 * dr);
			}
			else
			{

			}
		}
	}

#pragma omp	parallel for
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			for (k = 0; k < 13; k++)
			{
				tau_vis[i][j].f[k] = 0;
			}

			if (btype[i][j] == 1)
			{
				if (tau_ni[i][j] != 0)
				{
					tau_vis[i][j].f[5] = (taurr[i][j + 1] - taurr[i][j - 1]) / (2 * dr) + taurr[i][j] * pnir[i][j] / tau_ni[i][j] + taurr[i][j] / (dr * j) + taurz[i][j] * pniz[i][j] / tau_ni[i][j] + (taurz[i + 1][j] - taurz[i - 1][j]) / (2 * dz);
					tau_vis[i][j].f[6] = (taurtheta[i][j + 1] - taurtheta[i][j - 1]) / (2 * dr) + taurtheta[i][j] * pnir[i][j] / tau_ni[i][j] + taurtheta[i][j] / (dr * j) + tauthetaz[i][j] * pniz[i][j] / tau_ni[i][j] + (tauthetaz[i + 1][j] - tauthetaz[i - 1][j]) / (2 * dz);
					tau_vis[i][j].f[7] = (taurz[i][j + 1] - taurz[i][j - 1]) / (2 * dr) + taurz[i][j] * pnir[i][j] / tau_ni[i][j] + taurz[i][j] / (dr * j) + tauzz[i][j] * pniz[i][j] / tau_ni[i][j] + (tauzz[i + 1][j] - tauzz[i - 1][j]) / (2 * dz);
					tau_vis[i][j].f[12] = tau_vis[i][j].f[5] * U_bar[i][j].u[5] + tau_vis[i][j].f[6] * U_bar[i][j].u[6] + tau_vis[i][j].f[7] * U_bar[i][j].u[7];
				}

			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				if (tau_ni[i][j] != 0)
				{
					tau_vis[i][j].f[5] = taurr[i][j] * pnir[i][j] / tau_ni[i][j] + taurz[i][j] * pniz[i][j] / tau_ni[i][j] + (taurz[i + 1][j] - taurz[i - 1][j]) / (2 * dz);
					tau_vis[i][j].f[6] = taurtheta[i][j] * pnir[i][j] / tau_ni[i][j] + tauthetaz[i][j] * pniz[i][j] / tau_ni[i][j] + (tauthetaz[i + 1][j] - tauthetaz[i - 1][j]) / (2 * dz);
					tau_vis[i][j].f[7] = taurz[i][j] * pnir[i][j] / tau_ni[i][j] + tauzz[i][j] * pniz[i][j] / tau_ni[i][j] + (tauzz[i + 1][j] - tauzz[i - 1][j]) / (2 * dz);
					tau_vis[i][j].f[12] = tau_vis[i][j].f[5] * U_bar[i][j].u[5] + tau_vis[i][j].f[6] * U_bar[i][j].u[6] + tau_vis[i][j].f[7] * U_bar[i][j].u[7];

				}
			}

		}

	}
#pragma omp	parallel for
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			for (k = 0; k < 13; k++)
			{
				s_bar[i][j].f[k] += tau_vis[i][j].f[k];
			}
		}
	}
}

