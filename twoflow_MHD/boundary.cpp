#include "MHD.h"

using namespace std;

int world[ZMAX][RMAX];
int btype[ZMAX][RMAX];
int ptype[ZMAX][RMAX];
Boundary boundary_array[BND_NUM];
int bnd_size;
double bg_den = 2.41e8 * 1e9;
double inter_e_den = 1.1905e8 * 1e9;
double inter_pla_den = 6.68e7 * 1e9;


double cathod_cell;
double anode_cell;
double out_e_den;

double inlet_cell;
double out_cell;


int cath_cell_r[10 * ZMAX];
int cath_cell_z[10 * ZMAX];
int cath_cell_dir[10 * ZMAX];
int cath_num = 0;

int anode_cell_r[10 * ZMAX];
int anode_cell_z[10 * ZMAX];
int anode_cell_dir[10 * ZMAX];
int anode_num = 0;

int cera_cell_r[10 * ZMAX];
int cera_cell_z[10 * ZMAX];
int cera_cell_dir[10 * ZMAX];
int cera_num = 0;


int initial()
{
	register int i, j, k;
	
	cath_num = 0;
	anode_num = 0;
	cera_num = 0;

	memset(world, 0, sizeof(world));
	memset(btype, 0, sizeof(world));
	memset(ptype, 0, sizeof(world));

	_for(k, 0, bnd_size)
	{
		double dr = boundary_array[k].end.r - boundary_array[k].start.r;
		double dz = boundary_array[k].end.z - boundary_array[k].start.z;

		if (dr > dz)
		{
			double ins_z = dz / dr;
			i = 0;
			_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
			{
				world[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j] += 2;
				btype[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j] |= boundary_array[k].boundary_type;
				ptype[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j] |= boundary_array[k].physics_type;
				i++;
			}
		}
		else
		{
			double ins_r = dr / dz;
			j = 0;
			_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
			{
				world[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))] += 2;
				btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))] |= boundary_array[k].boundary_type;
				ptype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))] |= boundary_array[k].physics_type;
				j++;
			}
		}

	}

	fill_plasma(30 * scale, 30 * scale, 1);
	fill_plasma(0, 0, 110);
	//fill_plasma(120 * scale, 60 * scale, 210);
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (world[i][j] == 1)
			{
				btype[i][j] = 1;
			}
			else if (world[i][j] == 110)
			{
				btype[i][j] = 110;
			}
			else if (world[i][j] == 210)
			{
				btype[i][j] = 210;
			}
			//else if (btype[i][j] == (LEFT + UP) && (i != 0 && j != (RMAX - 1)))
			//{
			//	btype[i][j] = 1;
			//}
			//else if (btype[i][j] == (LEFT + DOWN) && (i != 0 && j != 0))
			//{
			//	btype[i][j] = 1;
			//}
			//else if (btype[i][j] == (RIGHT + DOWN) && (i != (ZMAX - 1) && j != 0))
			//{
			//	btype[i][j] = 1;
			//}
			//else if (btype[i][j] == (RIGHT + UP) && (i != (ZMAX - 1) && j != (RMAX - 1)))
			//{
			//	btype[i][j] = 1;
			//}
			//else if (btype[i][j] == (LEFT + LEFT))
			//{
			//	btype[i][j] = LEFT;
			//}
			//else if (btype[i][j] == (DOWN + DOWN))
			//{
			//	btype[i][j] = DOWN;
			//}
			//else if (btype[i][j] == (RIGHT + RIGHT))
			//{
			//	btype[i][j] = RIGHT;
			//}
			//else if (btype[i][j] == (UP + UP))
			//{
			//	btype[i][j] = UP;
			//}

		}
	}


	_for(k, 0, bnd_size)
	{
		double lr = boundary_array[k].end.r - boundary_array[k].start.r;
		double lz = boundary_array[k].end.z - boundary_array[k].start.z;

		if (lr == 0 || lz == 0) continue;
		if (lr > dz)
		{
			double ins_z = lz / lr;
			i = 0;
			_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
			{
				btype[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j] = judge_conner((int)(boundary_array[k].start.z + ceil(i * ins_z)), j);
				i++;
			}
		}
		else
		{
			double ins_r = lr / lz;
			j = 0;
			_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
			{
				btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))] = judge_conner(i,(int)(boundary_array[k].start.r + ceil(j * ins_r)));
				j++;
			}
		}

	}
	electrode_side_cell();


#ifdef BOUNDARY_DEBUG
	printf("cathod_cell = %lf\nanode_cell = %lf\n", cathod_cell, anode_cell);
#endif	
	matrix_int_to_csv((int**)world, ZMAX, RMAX, RMAX, (char*)(".\\output\\world.csv"));
	matrix_int_to_csv((int**)btype, ZMAX, RMAX, RMAX, (char*)(".\\output\\btype.csv"));
	matrix_int_to_csv((int**)ptype, ZMAX, RMAX, RMAX, (char*)(".\\output\\ptype.csv"));

	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] != 0 && btype[i][j] != 110 && btype[i][j] != 210)
			{
				MPDT[i][j].ne = bg_den * 0.01;
				MPDT[i][j].ni = bg_den * 0.01;
			}
			else
			{
				MPDT[i][j].ne = 0;
				MPDT[i][j].ni = 0;
			}

			MPDT[i][j].ver = 0;
			MPDT[i][j].vez = 0;
			MPDT[i][j].vetheta = 0;
			MPDT[i][j].vir = 0;
			MPDT[i][j].viz = 0;
			MPDT[i][j].vitheta = 0;
			MPDT[i][j].br = 0;
			MPDT[i][j].bz = 0;
			MPDT[i][j].btheta = 0;
			MPDT[i][j].pe = 0;
			MPDT[i][j].pi = 0;


			if (btype[i][j] != 0 && btype[i][j] != 110)
			{
				MPDT[i][j].ee = MPDT[i][j].pe / (MPDT[i][j].ne * ME) / (gamma - 1) + 0.5 * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
				MPDT[i][j].ei = MPDT[i][j].pi / (MPDT[i][j].ni * MI) / (gamma - 1) + 0.5 * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);

			}
			else
			{
				MPDT[i][j].ee = 0;
				MPDT[i][j].ei = 0;
			}


			if (btype[i][j] == 1)
			{
				atom[i][j].den = bg_den;
				//atom[i][j].eng = 1e-3 / (MI * bg_den);
			}
			else
			{
				atom[i][j].den = 0;
				atom[i][j].eng = 0;
			}
			atom[i][j].vr = 0;
			atom[i][j].vtheta = 0;
			atom[i][j].vz = 0;

			btheta[i][j] = 0;

			}
	}

	return 0;
}

int fill_plasma(int tz, int tr, int fill_n)
{
	queue<int> quer, quez;
	int tnr, tnz;
	quer.push(tr);
	quez.push(tz);
	while (!quer.empty() && !quez.empty())
	{
		tnr = quer.front();
		tnz = quez.front();
		quer.pop();
		quez.pop();
		if (tnr >= 0 && tnr < nr && tnz >= 0 && tnz < nz && world[tnz][tnr] == 0)
		{
			world[tnz][tnr] = fill_n;

			quer.push(tnr - 1);
			quez.push(tnz);

			quer.push(tnr + 1);
			quez.push(tnz);

			quer.push(tnr);
			quez.push(tnz - 1);

			quer.push(tnr);
			quez.push(tnz + 1);
		}
	}
	return 1;
}

//边界条件处理
void  boundary_condition()
{

	register int i, j, k;

	double eletron_emi = 0;
	double eletron_assi = 0;
	//固体边界

	_for(k, 0, bnd_size)
	{
		//if (boundary_array[k].physics_type == CATHODE_BOUNDARY)
		//{
		//	double dr = boundary_array[k].end.r - boundary_array[k].start.r;
		//	double dz = boundary_array[k].end.z - boundary_array[k].start.z;

		//	if (dr > dz)
		//	{
		//		double ins_z = dz / dr;
		//		i = 0;
		//		_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
		//		{
		//			MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ver = 0;
		//			MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vir = 0;
		//			MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vetheta = 0;
		//			MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vitheta = 0;
		//			MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vez = 0;
		//			MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].viz = 0;
		//			i++;
		//		}
		//	}
		//	else
		//	{
		//		double ins_r = dr / dz;
		//		j = 0;
		//		_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
		//		{
		//			MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ver = 0;
		//			MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vir = 0;
		//			MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vetheta = 0;
		//			MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vitheta = 0;
		//			MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vez = 0;
		//			MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].viz = 0;
		//			j++;
		//		}
		//	}
		//}

		//陶瓷边界二次电子发射
		if (boundary_array[k].physics_type == DIELECTRIC_SURFACE_BOUNDARY)
		{
			double lr = boundary_array[k].end.r - boundary_array[k].start.r;
			double lz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (lr > lz)
			{
				double ins_z = lz / lr;
				i = 0;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					if (MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ne < bg_den / 10)
					{
						MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ne = bg_den / 10;
					}
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ver = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vir = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vetheta = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vitheta = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vez = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].viz = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].pe = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].pi = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ee = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ei = 0;
					i++;
				}
			}
			else
			{
				double ins_r = lr / lz;
				j = 0;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					if (MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne < bg_den / 10)
					{
						MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne = bg_den / 10;
					}
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ver = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vir = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vetheta = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vitheta = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vez = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].viz = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].pe = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].pi = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ee = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ei = 0;
					j++;
				}
			}
		}

		if (boundary_array[k].physics_type == ANODE_BOUNDARY)
		{
			double lr = boundary_array[k].end.r - boundary_array[k].start.r;
			double lz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (lr > lz)
			{
				double ins_z = lz / lr;
				i = 0;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ver = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vir = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vetheta = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vitheta = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vez = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].viz = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].pe = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].pi = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ee = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ei = 0;
					i++;
				}
			}
			else
			{
				double ins_r = lr / lz;
				j = 0;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ver = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vir = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vetheta = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vitheta = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vez = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].viz = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].pe = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].pi = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ee = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ei = 0;
					j++;
				}
			}
		}

	}

	updata_ion_edge();
	updata_electron_edge();

	_for(k, 0, bnd_size)
	{
		if (boundary_array[k].physics_type == CATHODE_BOUNDARY)
		{
			double lr = boundary_array[k].end.r - boundary_array[k].start.r;
			double lz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (lr > lz)
			{
				double ins_z = lz / lr;
				i = 0;
				_for(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					//if (btype[(int)(boundary_array[k].start.z + ceil(i * ins_z)) + 1][j] != 1) continue;
					//MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z)) + 1][j].neq = inter_e_den / scale;
					//MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].neq = inter_e_den / scale;
					//MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vnqz = 0.2 * 30000;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ne = inter_e_den / scale;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vez = 0.2 * 30000;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ee = 0.5 * sqr(0.21 * 30000);
					i++;
				}
			}
			else
			{
				double ins_r = lr / lz;
				j = 0;
				_for(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					if (btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))] == DOWN)
					{
						//if (btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) + 1] != 1) continue;
						//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) + 1].neq = inter_e_den / scale;
						//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].neq = inter_e_den / scale;
						//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vnqr = 0.2 * 30000;
						MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne = inter_e_den / scale;
						MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ver = 0.2 * 30000;
						MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ee = 0.5 * sqr(0.21 * 30000);
					}
					else
					{
						//if (btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1] != 1) continue;
						//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1].neq = inter_e_den / scale;
						//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].neq = inter_e_den / scale;
						//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vnqr = -0.2 * 30000;
						MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne = inter_e_den / scale;
						MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ver = -0.2 * 30000;
						MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ee = 0.5 * sqr(0.21 * 30000);
					}
					j++;
				}
			}
		}


		if (boundary_array[k].physics_type == VACCUM_BOUNDARY)
		{
			double lr = boundary_array[k].end.r - boundary_array[k].start.r;
			double lz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (lr > lz)
			{
				double ins_z = lz / lr;
				i = 0;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					//if (MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ne < 1e4) continue;
					//if (btype[(int)(boundary_array[k].start.z + ceil(i * ins_z)) - 1][j] != 1) continue;
					//MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z)) - 1][j].ne /= 2;
					//MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z)) - 1][j].ni /= 2;
					//MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ne /= 2;
					//MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ni /= 2;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ne = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ver = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vetheta = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vez = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ni = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vir = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vitheta = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].viz = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ee = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ei = 0;
					i++;
				}
			}
			else
			{
				double ins_r = lr / lz;
				j = 0;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					//if (MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne < 1e4) continue;
					//if (btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1] != 1) continue;
					//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1].ne /= 2;
					//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1].ni /= 2;
					//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne /= 2;
					//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ni /= 2;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne = bg_den / scale;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ver = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vetheta = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vez = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ni = bg_den / scale;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vir = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vitheta = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].viz = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ee = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ei = 0;
					j++;
				}
			}
		}






		/*if (boundary_array[k].physics_type == INLET)
		{
			double lr = boundary_array[k].end.r - boundary_array[k].start.r;
			double lz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (lr > lz)
			{
				double ins_z = lz / lr;
				i = 0;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					
					i++;
				}
			}
			else
			{
				double ins_r = lr / lz;
				j = 0;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{

					j++;
				}
			}
		}*/

	}
}

int judge_conner(int i,int j)
{
	uint8_t state;
	uint8_t state_left = 8;
	uint8_t state_right = 4;
	uint8_t state_down = 2;
	uint8_t state_up = 1;
	state = 0;
	if (i == 0)
	{
		state += state_left;
	}
	else
	{
		if (btype[i - 1][j] == 0 || btype[i - 1][j] == 110)
		{
			state += state_left;
		}
	}
	
	if(i == ZMAX - 1)
	{
		state += state_right;
	}
	else
	{
		if (btype[i + 1][j] == 0 || btype[i + 1][j] == 110)
		{
			state += state_right;
		}
	}

	if (j == 0)
	{
		state += state_down;
	}
	else
	{
		if (btype[i][j - 1] == 0 || btype[i][j - 1] == 110)
		{
			state += state_down;
		}
	}

	if (j == RMAX - 1)
	{
		state += state_up;
	}
	else
	{
		if (btype[i][j + 1] == 0 || btype[i][j + 1] == 110)
		{
			state += state_up;
		}
	}

	if (state == state_left)
	{
		return LEFT;
	}
	else if (state == state_right)
	{
		return RIGHT;
	}
	else if (state == state_up)
	{
		return UP;
	}
	else if (state == state_down)
	{
		return DOWN;
	}
	else if (state == state_left + state_up)
	{
		return LEFT + UP;
	}
	else if (state == state_left + state_down)
	{
		return LEFT + DOWN;
	}
	else if (state == state_right + state_up)
	{
		return RIGHT + UP;
	}
	else if (state == state_right + state_down)
	{
		return RIGHT + DOWN;
	}
	else if (state == state_left + state_down + state_up)
	{
		return LEFT + DOWN;
	}
	else if (state == state_right + state_down + +state_up)
	{
		return RIGHT + DOWN;
	}
	else if (state == 0)
	{
		return btype[i][j];
	}
	else
	{
		printf("boundary condition error %d %d %d\n",state,i,j);
		return 0;
	}

	return 0;
}

void electrode_side_cell()
{
	register int i, j, k;
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				if (is_cathode_ngh(i, j))
				{
					break;
				}

				if (is_anode_ngh(i, j))
				{
					break;
				}

				if (is_cera_ngh(i, j))
				{
					break;
				}
			}
		}
	}
}

int is_cathode_ngh(int i, int j)
{
	if (i - 1 > 0)
	{
		if (ptype[i - 1][j] == CATHODE_BOUNDARY)
		{
			cath_cell_z[cath_num] =i;
			cath_cell_r[cath_num] =j;
			cath_cell_dir[cath_num] = LEFT;
			cath_num++;

			return 1;
		}
	}

	if (i + 1 < ZMAX)
	{
		if (ptype[i + 1][j] == CATHODE_BOUNDARY )
		{
			cath_cell_z[cath_num] = i;
			cath_cell_r[cath_num] = j;
			cath_cell_dir[cath_num] = RIGHT;
			cath_num++;

			return 1;
		}
	}

	if (j - 1 > 0)
	{
		if (ptype[i][j - 1] == CATHODE_BOUNDARY)
		{
			cath_cell_z[cath_num] = i;
			cath_cell_r[cath_num] = j;
			cath_cell_dir[cath_num] = DOWN;
			cath_num++;

			return 1;
		}
	}

	if (i + 1 < ZMAX)
	{
		if (ptype[i + 1][j] == CATHODE_BOUNDARY )
		{
			cath_cell_z[cath_num] = i;
			cath_cell_r[cath_num] = j;
			cath_cell_dir[cath_num] = UP;
			cath_num++;

			return 1;
		}
	}
	
	return 0;
}


int is_anode_ngh(int i, int j)
{
	if (i - 1 > 0)
	{
		if (ptype[i - 1][j] == ANODE_BOUNDARY )
		{
			anode_cell_z[anode_num] = i;
			anode_cell_r[anode_num] = j;
			anode_cell_dir[anode_num] = LEFT;
			anode_num++;

			return 1;
		}
	}

	if (i + 1 < ZMAX)
	{
		if (ptype[i + 1][j] == ANODE_BOUNDARY )
		{
			anode_cell_z[anode_num] = i;
			anode_cell_r[anode_num] = j;
			anode_cell_dir[anode_num] = RIGHT;
			anode_num++;

			return 1;
		}
	}

	if (j - 1 > 0)
	{
		if (ptype[i][j - 1] == ANODE_BOUNDARY)
		{
			anode_cell_z[anode_num] = i;
			anode_cell_r[anode_num] = j;
			anode_cell_dir[anode_num] = DOWN;
			anode_num++;

			return 1;
		}
	}

	if (i + 1 < ZMAX)
	{
		if (ptype[i + 1][j] == ANODE_BOUNDARY)
		{
			anode_cell_z[anode_num] = i;
			anode_cell_r[anode_num] = j;
			anode_cell_dir[anode_num] = UP;
			anode_num++;

			return 1;
		}
	}

	return 0;
}



int is_cera_ngh(int i, int j)
{
	if (i - 1 > 0)
	{
		if (ptype[i - 1][j] == DIELECTRIC_SURFACE_BOUNDARY)
		{
			cera_cell_z[cera_num] = i;
			cera_cell_r[cera_num] = j;
			cera_cell_dir[cera_num] = LEFT;
			cera_num++;

			return 1;
		}
	}

	if (i + 1 < ZMAX)
	{
		if (ptype[i + 1][j] == DIELECTRIC_SURFACE_BOUNDARY)
		{
			cera_cell_z[cera_num] = i;
			cera_cell_r[cera_num] = j;
			cera_cell_dir[cera_num] = RIGHT;
			cera_num++;

			return 1;
		}
	}

	if (j - 1 > 0)
	{
		if (ptype[i][j - 1] == DIELECTRIC_SURFACE_BOUNDARY)
		{
			cera_cell_z[cera_num] = i;
			cera_cell_r[cera_num] = j;
			cera_cell_dir[cera_num] = DOWN;
			cera_num++;

			return 1;
		}
	}

	if (i + 1 < ZMAX)
	{
		if (ptype[i + 1][j] == DIELECTRIC_SURFACE_BOUNDARY)
		{
			cera_cell_z[cath_num] = i;
			cera_cell_r[cath_num] = j;
			cera_cell_dir[cath_num] = UP;
			cera_num++;

			return 1;
		}
	}

	return 0;
}

void updata_atom_edge()
{
	int i, j;
	for (int k = 0; k < cath_num; k++)
	{
		i = cath_cell_z[k];
		j = cath_cell_r[k];
		if (atom[i][j].vz < 0 && cath_cell_dir[k] == LEFT)
		{
			atom[i][j].vz = -atom[i][j].vz;
		}

		if (atom[i][j].vz > 0 && cath_cell_dir[k] == RIGHT)
		{
			atom[i][j].vz = -atom[i][j].vz;
		}
		
		if (atom[i][j].vr < 0 && cath_cell_dir[k] == DOWN)
		{
			atom[i][j].vr = -atom[i][j].vr;
		}

		if (atom[i][j].vr > 0 && cath_cell_dir[k] == UP)
		{
			atom[i][j].vr = -atom[i][j].vr;
		}
	}

	for (int k = 0; k < anode_num; k++)
	{
		i = anode_cell_z[k];
		j = anode_cell_r[k];
		if (atom[i][j].vz < 0 && anode_cell_dir[k] == LEFT)
		{
			atom[i][j].vz = -atom[i][j].vz;
		}

		if (atom[i][j].vz > 0 && anode_cell_dir[k] == RIGHT)
		{
			atom[i][j].vz = -atom[i][j].vz;
		}

		if (atom[i][j].vr < 0 && anode_cell_dir[k] == DOWN)
		{
			atom[i][j].vr = -atom[i][j].vr;
		}

		if (atom[i][j].vr > 0 && anode_cell_dir[k] == UP)
		{
			atom[i][j].vr = -atom[i][j].vr;
		}
	}

	for (int k = 0; k < cera_num; k++)
	{
		i = cera_cell_z[k];
		j = cera_cell_r[k];
		if (atom[i][j].vz < 0 && cera_cell_dir[k] == LEFT)
		{
			atom[i][j].vz = -atom[i][j].vz;
		}

		if (atom[i][j].vz > 0 && cera_cell_dir[k] == RIGHT)
		{
			atom[i][j].vz = -atom[i][j].vz;
		}

		if (atom[i][j].vr < 0 && cera_cell_dir[k] == DOWN)
		{
			atom[i][j].vr = -atom[i][j].vr;
		}

		if (atom[i][j].vr > 0 && cera_cell_dir[k] == UP)
		{
			atom[i][j].vr = -atom[i][j].vr;
		}
	}
}

void updata_ion_edge()
{
	int i, j;
	for (int k = 0; k < cath_num; k++)
	{
		i = cath_cell_z[k];
		j = cath_cell_r[k];
		if (MPDT[i][j].viz < 0 && cath_cell_dir[k] == LEFT)
		{
			MPDT[i][j].viz = -MPDT[i][j].viz;
		}

		if (MPDT[i][j].viz > 0 && cath_cell_dir[k] == RIGHT)
		{
			MPDT[i][j].viz = -MPDT[i][j].viz;
		}

		if (MPDT[i][j].vir < 0 && cath_cell_dir[k] == DOWN)
		{
			MPDT[i][j].vir = -MPDT[i][j].vir;
		}

		if (MPDT[i][j].vir > 0 && cath_cell_dir[k] == UP)
		{
			MPDT[i][j].vir = -MPDT[i][j].vir;
		}
	}

	for (int k = 0; k < anode_num; k++)
	{
		i = anode_cell_z[k];
		j = anode_cell_r[k];
		if (MPDT[i][j].viz < 0 && anode_cell_dir[k] == LEFT)
		{
			MPDT[i][j].viz = -MPDT[i][j].viz;
		}

		if (MPDT[i][j].viz > 0 && anode_cell_dir[k] == RIGHT)
		{
			MPDT[i][j].viz = -MPDT[i][j].viz;
		}

		if (MPDT[i][j].vir < 0 && anode_cell_dir[k] == DOWN)
		{
			MPDT[i][j].vir = -MPDT[i][j].vir;
		}

		if (MPDT[i][j].vir > 0 && anode_cell_dir[k] == UP)
		{
			MPDT[i][j].vir = -MPDT[i][j].vir;
		}
	}

	for (int k = 0; k < cera_num; k++)
	{
		i = cera_cell_z[k];
		j = cera_cell_r[k];
		if (MPDT[i][j].viz < 0 && cera_cell_dir[k] == LEFT)
		{
			MPDT[i][j].viz = -MPDT[i][j].viz;
		}

		if (MPDT[i][j].viz > 0 && cera_cell_dir[k] == RIGHT)
		{
			MPDT[i][j].viz = -MPDT[i][j].viz;
		}

		if (MPDT[i][j].vir < 0 && cera_cell_dir[k] == DOWN)
		{
			MPDT[i][j].vir = -MPDT[i][j].vir;
		}

		if (MPDT[i][j].vir > 0 && cera_cell_dir[k] == UP)
		{
			MPDT[i][j].vir = -MPDT[i][j].vir;
		}
	}
}

void updata_electron_edge()
{
	int i, j;
	for (int k = 0; k < cath_num; k++)
	{
		i = cath_cell_z[k];
		j = cath_cell_r[k];
		if (MPDT[i][j].vez < 0 && cath_cell_dir[k] == LEFT)
		{
			MPDT[i][j].vez = -MPDT[i][j].vez;
		}

		if (MPDT[i][j].vez > 0 && cath_cell_dir[k] == RIGHT)
		{
			MPDT[i][j].vez = -MPDT[i][j].vez;
		}

		if (MPDT[i][j].ver < 0 && cath_cell_dir[k] == DOWN)
		{
			MPDT[i][j].ver = -MPDT[i][j].ver;
		}

		if (MPDT[i][j].ver > 0 && cath_cell_dir[k] == UP)
		{
			MPDT[i][j].ver = -MPDT[i][j].ver;
		}
	}

	/*for (int k = 0; k < anode_num; k++)
	{
		i = anode_cell_z[k];
		j = anode_cell_r[k];
		if (MPDT[i][j].vez < 0 && anode_cell_dir[k] == LEFT)
		{
			MPDT[i][j].vez = -MPDT[i][j].vez;
		}

		if (MPDT[i][j].vez > 0 && anode_cell_dir[k] == RIGHT)
		{
			MPDT[i][j].vez = -MPDT[i][j].vez;
		}

		if (MPDT[i][j].ver < 0 && anode_cell_dir[k] == DOWN)
		{
			MPDT[i][j].ver = -MPDT[i][j].ver;
		}

		if (MPDT[i][j].ver > 0 && anode_cell_dir[k] == UP)
		{
			MPDT[i][j].ver = -MPDT[i][j].ver;
		}
	}*/

	for (int k = 0; k < cera_num; k++)
	{
		i = cera_cell_z[k];
		j = cera_cell_r[k];
		if (MPDT[i][j].vez < 0 && cera_cell_dir[k] == LEFT)
		{
			MPDT[i][j].vez = -MPDT[i][j].vez;
		}

		if (MPDT[i][j].vez > 0 && cera_cell_dir[k] == RIGHT)
		{
			MPDT[i][j].vez = -MPDT[i][j].vez;
		}

		if (MPDT[i][j].ver < 0 && cera_cell_dir[k] == DOWN)
		{
			MPDT[i][j].ver = -MPDT[i][j].ver;
		}

		if (MPDT[i][j].ver > 0 && cera_cell_dir[k] == UP)
		{
			MPDT[i][j].ver = -MPDT[i][j].ver;
		}
	}
}