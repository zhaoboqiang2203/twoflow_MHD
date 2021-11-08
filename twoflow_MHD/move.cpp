#include "MHD.h"

void move()
{
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nz; j++)
		{
			MPDT[i][j].ver = -QE * Er[i][j] / ME * dt;
			MPDT[i][j].vez = -QE * Ez[i][j] / ME * dt;

			MPDT[i][j].ee = 0.5 * MPDT[i][j].ne * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
			//MPDT[i][j].pe = MPDT[i][j].ne * MPDT[i][j].ee;

			MPDT[i][j].vir = QE * Er[i][j] / MI * dt;
			MPDT[i][j].viz = QE * Ez[i][j] / MI * dt;

			MPDT[i][j].ei = 0.5 * MPDT[i][j].ni * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);
			//MPDT[i][j].pi = MPDT[i][j].ni * MPDT[i][j].ei;


		}
	}
}