#include "atomFluid.h"



void atom_flow()
{
	struct _AF Qr;
	struct _AF Qz;
	register int i, j, k;

	
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{

			Uq[i][j] = cal_atom_u(i, j);
			Fqr[i][j] = cal_atom_Fqr(Uq[i][j]);
			Fqz[i][j] = cal_atom_Fqz(Uq[i][j]);
			sq[i][j] = cal_atom_s(Uq[i][j]);

			Uq_bar[i][j].z = i;
			Uq_bar[i][j].r = j;
			Uq_bar2[i][j].z = i;
			Uq_bar2[i][j].r = j;

		}
	}



	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{

			if (btype[i][j] == 1)
			{
				for (k = 0; k < 5; k++)
				{
					Uq_bar[i][j].f[k] = Uq[i][j].f[k] - dat / dr * (Fqr[i][j + 1].f[k] - Fqr[i][j].f[k]) - dat / dz * (Fqz[i + 1][j].f[k] - Fqz[i][j].f[k]) + dat * sq[i][j].f[k];
				}

			}
			else if (btype[i][j] == 0)
			{
				for (k = 0; k < 5; k++)
				{
					Uq_bar[i][j].f[k] = 0;
				}
			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				for (k = 0; k < 5; k++)
				{
					Uq_bar[i][j].f[k] = Uq[i][j].f[k] - dat / dz * (Fqz[i + 1][j].f[k] - Fqz[i][j].f[k]) + dat * sq[i][j].f[k];
				}
			}
			else
			{
				for (k = 0; k < 5; k++)
				{
					Uq_bar[i][j].f[k] = Uq[i][j].f[k];
				}
			}

			Fqr_bar[i][j] = cal_atom_Fqr(Uq_bar[i][j]);
			Fqz_bar[i][j] = cal_atom_Fqz(Uq_bar[i][j]);
			sq_bar[i][j] = cal_atom_s(Uq_bar[i][j]);
		}
	}



	//矫正步
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				for (k = 0; k < 5; k++)
				{
					Uq_bar2[i][j].f[k] = 0.5 * (Uq[i][j].f[k] + Uq_bar[i][j].f[k] - dat / dr * (Fqr_bar[i][j].f[k] - Fqr_bar[i][j - 1].f[k]) - dat / dz * (Fqz_bar[i][j].f[k] - Fqz_bar[i - 1][j].f[k]) + dat * sq_bar[i][j].f[k]);
				}
			}
			else if (btype[i][j] == 0)
			{
				for (k = 0; k < 5; k++)
				{
					Uq_bar2[i][j].f[k] = 0;
				}
			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				for (k = 0; k < 5; k++)
				{
					Uq_bar2[i][j].f[k] = 0.5 * (Uq[i][j].f[k] + Uq_bar[i][j].f[k] - dat / dz * (Fqz_bar[i][j].f[k] - Fqz_bar[i - 1][j].f[k]) + dat * sq_bar[i][j].f[k]);
				}
			}
			else
			{
				for (k = 0; k < 5; k++)
				{
					Uq_bar2[i][j].f[k] = Uq_bar[i][j].f[k];
				}
			}

		}
	}

	//计算				

	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				Qz = arti_atom_vis(atom[i + 1][j], atom[i][j], atom[i - 1][j]);
				Qr = arti_atom_vis(atom[i][j + 1], atom[i][j], atom[i][j - 1]);

				Uq[i][j].f[0] = Uq_bar2[i][j].f[0] + Qr.f[0] / 4 * (atom[i][j + 1].den - 2 * atom[i][j].den + atom[i][j - 1].den) + Qz.f[0] / 4 * (atom[i + 1][j].den - 2 * atom[i][j].den + atom[i - 1][j].den);
				Uq[i][j].f[1] = Uq_bar2[i][j].f[1] + Qr.f[1] / 4 * (atom[i][j + 1].vr - 2 * atom[i][j].vr + atom[i][j - 1].vr) + Qz.f[1] / 4 * (atom[i + 1][j].vr - 2 * atom[i][j].vr + atom[i - 1][j].vr);
				Uq[i][j].f[2] = Uq_bar2[i][j].f[2] + Qr.f[2] / 4 * (atom[i][j + 1].vtheta - 2 * atom[i][j].vtheta + atom[i][j - 1].vtheta) + Qz.f[2] / 4 * (atom[i + 1][j].vtheta - 2 * atom[i][j].vtheta + atom[i - 1][j].vtheta);
				Uq[i][j].f[3] = Uq_bar2[i][j].f[3] + Qr.f[3] / 4 * (atom[i][j + 1].vz - 2 * atom[i][j].vz + atom[i][j - 1].vz) + Qz.f[3] / 4 * (atom[i + 1][j].vz - 2 * atom[i][j].vz + atom[i - 1][j].vz);
				Uq[i][j].f[4] = Uq_bar2[i][j].f[4] + Qr.f[4] / 4 * (atom[i][j + 1].eng - 2 * atom[i][j].eng + atom[i][j - 1].eng) + Qz.f[4] / 4 * (atom[i + 1][j].eng - 2 * atom[i][j].eng + atom[i - 1][j].eng);



			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				Qz = arti_atom_vis(atom[i + 1][j], atom[i][j], atom[i - 1][j]);
				Qr = arti_atom_vis(atom[i][j + 2], atom[i][j + 1], atom[i][j]);

				Uq[i][j].f[0] = Uq_bar2[i][j].f[0] + Qr.f[0] / 2 * (atom[i][j + 1].den - 2 * atom[i][j].den + atom[i][j - 1].den);
				Uq[i][j].f[1] = Uq_bar2[i][j].f[1] + Qr.f[1] / 2 * (atom[i][j + 1].vr - 2 * atom[i][j].vr + atom[i][j - 1].vr);
				Uq[i][j].f[2] = Uq_bar2[i][j].f[2] + Qr.f[2] / 2 * (atom[i][j + 1].vtheta - 2 * atom[i][j].vtheta + atom[i][j - 1].vtheta);
				Uq[i][j].f[3] = Uq_bar2[i][j].f[3] + Qr.f[3] / 2 * (atom[i][j + 1].vz - 2 * atom[i][j].vz + atom[i][j - 1].vz) ;
				Uq[i][j].f[4] = Uq_bar2[i][j].f[4] + Qr.f[4] / 2 * (atom[i][j + 1].eng - 2 * atom[i][j].eng + atom[i][j - 1].eng) ;

			}
			else 
			{
				for (k = 0; k < 5; k++)
				{
					Uq[i][j].f[k] = Uq_bar2[i][j].f[k];
				}
			}
		}
	}

	//更新边界信息

	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			atom[i][j].den = Uq[i][j].f[0];
			atom[i][j].vr = Uq[i][j].f[1];
			atom[i][j].vtheta = Uq[i][j].f[2];
			atom[i][j].vz = Uq[i][j].f[3];
			atom[i][j].eng = Uq[i][j].f[4];
		}
	}


	return;
}

struct _AF cal_atom_u(int i, int j)
{
	struct _AF uij;
	uij.r = j;
	uij.z = i;
	uij.f[0] = atom[i][j].den;
	uij.f[1] = atom[i][j].vr;
	uij.f[2] = atom[i][j].vtheta;
	uij.f[3] = atom[i][j].vz;
	uij.f[4] = atom[i][j].eng;

	return uij;
}

struct _AF cal_atom_Fqr(struct _AF uij)
{
	struct _AF fij;

	double den = uij.f[0];
	double vr = uij.f[1];
	double vtheta = uij.f[2];
	double vz = uij.f[3];
	double eng = uij.f[4];


	fij.r = uij.r;
	fij.z = uij.z;
	fij.f[0] = den * vr;

	fij.f[1] = -(gamma - 1) * (0.5 * (sqr(vr) + sqr(vtheta) + sqr(vz)) - eng);
	fij.f[2] = 0;
	fij.f[3] = 0;


	fij.f[4] = -vr * (gamma - 1) * (0.5 * (sqr(vr) + sqr(vtheta) + sqr(vz)) - eng);

	return fij;
}

struct _AF cal_atom_Fqz(struct _AF uij)
{
	struct _AF fij;

	double den = uij.f[0];
	double vr = uij.f[1];
	double vtheta = uij.f[2];
	double vz = uij.f[3];
	double eng = uij.f[4];

	

	fij.r = uij.r;
	fij.z = uij.z;

	fij.f[0] = den * vz;

	fij.f[1] = 0;
	fij.f[2] = 0;
	fij.f[3] = -(gamma - 1) * (0.5 * (vr * vr + vtheta * vtheta + vz * vz) - eng);


	fij.f[4] = -vz * (gamma - 1) * (0.5 * (vr * vr + vtheta * vtheta + vz * vz) - eng);
	return fij;
}

struct _AF cal_atom_s(struct _AF uij)
{

	struct _AF fij;

	double den = uij.f[0];
	double vr = uij.f[1];
	double vtheta = uij.f[2];
	double vz = uij.f[3];
	double eng = uij.f[4];

	double r1 = 0;
	r1 = -uij.r * dr;

	fij.r = uij.r;
	fij.z = uij.z;
	fij.f[0] = r1 * den * vr;

	fij.f[1] = -r1 * (gamma - 1) * (0.5 * (sqr(vr) + sqr(vtheta) + sqr(vz)) - eng);
	fij.f[2] = 0;
	fij.f[3] = 0;


	fij.f[4] = -r1 * vr * (gamma - 1) * (0.5 * (sqr(vr) + sqr(vtheta) + sqr(vz)) - eng);

	return fij;
}

struct _AF arti_atom_vis(struct _particle np, struct  _particle n, struct  _particle nn)
{
	struct _AF Q;
	double eta = 0.03;

	for (int k = 0; k < 5; k++)
	{
		Q.f[k] = 0;
	}

	if (abs(np.den - n.den) + abs(n.den - nn.den) != 0)
	{
		Q.f[0] = eta * abs(abs(np.den - n.den) - abs(n.den - nn.den)) / abs(abs(np.den - n.den) + abs(n.den - nn.den));
	}


	Q.f[1] = Q.f[0];
	Q.f[2] = Q.f[0];
	Q.f[3] = Q.f[0];
	Q.f[4] = Q.f[0];


	return Q;
}

void atom_boundary_condition()
{
	register int i, j, k;

	double pin = 400;
	double flow_rate = 0.1;
	int ncell = 0;
	double scell = 0;
	double atomden = 0;
	_for(k, 0, bnd_size)
	{
		if (boundary_array[k].physics_type == INLET)
		{
			double lr = boundary_array[k].end.r - boundary_array[k].start.r;
			double lz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (lr > lz)
			{
				double ins_z = lz / lr;
				i = 0;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					ncell++;
					scell += PI * sqr(dr) * PI * (2 * j - 1);
					i++;
				}
			}
			else
			{
				double ins_r = lr / lz;
				j = 0;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					ncell++;
					scell += PI * sqr(dr) * PI * (2 * j - 1);
					j++;
				}
			}
		}
	}

	atomden = flow_rate / REL_MASS * NA / scell;

	_for(k, 0, bnd_size)
	{
		if (boundary_array[k].physics_type == INLET)
		{
			double lr = boundary_array[k].end.r - boundary_array[k].start.r;
			double lz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (lr > lz)
			{
				double ins_z = lz / lr;
				i = 0;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					atom[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].den = atomden;
					atom[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vz = (gamma - 1) * pin / (MI * atomden);
					atom[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].eng = 0.6 * sqr(atom[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vz);
					i++;
				}
			}
			else
			{
				double ins_r = lr / lz;
				j = 0;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					atom[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].den = atomden;
					atom[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vr = (gamma - 1) * pin / (MI * atomden);
					atom[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].eng = 0.6 * sqr((gamma - 1) * pin / (MI * atomden));
					j++;
				}
			}
		}
	}
}