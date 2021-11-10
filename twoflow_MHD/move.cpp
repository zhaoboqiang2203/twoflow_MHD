#include "MHD.h"

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


			//µç×ÓÀë×Ó¿âÂ×Åö×²


		}
	}
}