#include "MHD.h"

using namespace std;

double world[ZMAX][RMAX];
double btype[ZMAX][RMAX];
double ptype[ZMAX][RMAX];
Boundary boundary_array[BND_NUM];
//double bg_den = 2.41e8 * 1e9;
//double inter_e_den = 1.1905e8 * 1e9;
//double inter_pla_den = 6.68e7 * 1e9;
double bg_den = 2.41e8;
double inter_e_den = 1.1905e8;
double inter_pla_den = 6.68e7;
int initial()
{
	register int i, j, k;
	//		 ___________________________________________
	//		 |  40  | 50      |            100			|
	//		 |				  |							|
	//		 |50			  |							|
	//		 |				  |							|
	//_______|				  |							|
	//|  10					  |							|
	//|10					  |							|
	//|___30______		      |							|
	//		     |		      |							|
	//           |20	      |							|
	//		     |___70_______|___________100___________|

	boundary_array[0].physics_type = CYLINDRICAL_AXIS;
	boundary_array[0].start.r = 0 * scale;
	boundary_array[0].start.z = 30 * scale;
	boundary_array[0].end.r = 0 * scale;
	boundary_array[0].end.z = 200 * scale;
	boundary_array[0].bnd_dir = Z_DIR;
	boundary_array[0].boundary_type = DOWN;

	boundary_array[1].physics_type = VACCUM_BOUNDARY;
	boundary_array[1].start.r = 0 * scale;
	boundary_array[1].start.z = 200 * scale;
	boundary_array[1].end.r = 80 * scale;
	boundary_array[1].end.z = 200 * scale;
	boundary_array[1].bnd_dir = R_DIR;
	boundary_array[1].boundary_type = RIGHT;

	boundary_array[2].physics_type = VACCUM_BOUNDARY;
	boundary_array[2].start.r = 80 * scale;
	boundary_array[2].start.z = 100 * scale;
	boundary_array[2].end.r = 80 * scale;
	boundary_array[2].end.z = 200 * scale;
	boundary_array[2].bnd_dir = Z_DIR;
	boundary_array[2].boundary_type = UP;

	boundary_array[3].physics_type = ANODE_BOUNDARY;
	boundary_array[3].start.r = 80 * scale;
	boundary_array[3].start.z = 50 * scale;
	boundary_array[3].end.r = 80 * scale;
	boundary_array[3].end.z = 100 * scale;
	boundary_array[3].bnd_dir = Z_DIR;
	boundary_array[3].boundary_type = UP;

	boundary_array[4].physics_type = DIELECTRIC_SURFACE_BOUNDARY;
	boundary_array[4].start.r = 80 * scale;
	boundary_array[4].start.z = 10 * scale;
	boundary_array[4].end.r = 80 * scale;
	boundary_array[4].end.z = 50 * scale;
	boundary_array[4].bnd_dir = Z_DIR;
	boundary_array[4].boundary_type = UP;

	boundary_array[5].physics_type = DIELECTRIC_SURFACE_BOUNDARY;
	boundary_array[5].start.r = 30 * scale;
	boundary_array[5].start.z = 10 * scale;
	boundary_array[5].end.r = 80 * scale;
	boundary_array[5].end.z = 10 * scale;
	boundary_array[5].bnd_dir = R_DIR;
	boundary_array[5].boundary_type = LEFT;

	boundary_array[6].physics_type = DIELECTRIC_SURFACE_BOUNDARY;
	boundary_array[6].start.r = 30 * scale;
	boundary_array[6].start.z = 0 * scale;
	boundary_array[6].end.r = 30 * scale;
	boundary_array[6].end.z = 10 * scale;
	boundary_array[6].bnd_dir = Z_DIR;
	boundary_array[6].boundary_type = UP;

	boundary_array[7].physics_type = DIELECTRIC_SURFACE_BOUNDARY;
	boundary_array[7].start.r = 20 * scale;
	boundary_array[7].start.z = 0 * scale;
	boundary_array[7].end.r = 30 * scale;
	boundary_array[7].end.z = 0 * scale;
	boundary_array[7].bnd_dir = R_DIR;
	boundary_array[7].boundary_type = LEFT;

	boundary_array[8].physics_type = CATHODE_BOUNDARY;
	boundary_array[8].start.r = 20 * scale;
	boundary_array[8].start.z = 0 * scale;
	boundary_array[8].end.r = 20 * scale;
	boundary_array[8].end.z = 30 * scale;
	boundary_array[8].bnd_dir = Z_DIR;
	boundary_array[8].boundary_type = DOWN;

	boundary_array[9].physics_type = INLET;
	boundary_array[9].start.r = 0 * scale;
	boundary_array[9].start.z = 30 * scale;
	boundary_array[9].end.r = 20 * scale;
	boundary_array[9].end.z = 30 * scale;
	boundary_array[9].bnd_dir = R_DIR;
	boundary_array[9].boundary_type = LEFT;



	memset(world, 0, sizeof(world));
	memset(btype, 0, sizeof(world));
	memset(ptype, 0, sizeof(world));

	_for(k, 0, BND_NUM)
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
				btype[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j] += boundary_array[k].boundary_type;
				ptype[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j] += boundary_array[k].physics_type;
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
				btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))] += boundary_array[k].boundary_type;
				ptype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))] += boundary_array[k].physics_type;
				j++;
			}
		}

	}

	fill_plasma(30 * scale, 30 * scale, 1);
	fill_plasma(0, 0, 110);
	for (int i = 0; i < nz; i++)
	{
		for (int j = 0; j < nr; j++)
		{
			if (world[i][j] == 1)
			{
				btype[i][j] = 1;
			}
			else if (world[i][j] == 110)
			{
				btype[i][j] = 110;
			}
			else if (btype[i][j] == (LEFT + UP) && (i != 0 && j != (RMAX - 1)))
			{
				btype[i][j] = 1;
			}
			else if (btype[i][j] == (LEFT + DOWN) && (i != 0 && j != 0))
			{
				btype[i][j] = 1;
			}
			else if (btype[i][j] == (RIGHT + DOWN) && (i != (ZMAX - 1) && j != 0))
			{
				btype[i][j] = 1;
			}
			else if (btype[i][j] == (RIGHT + UP) && (i != (ZMAX - 1) && j != (RMAX - 1)))
			{
				btype[i][j] = 1;
			}
			else if (btype[i][j] == (LEFT + LEFT))
			{
				btype[i][j] = LEFT;
			}
			else if (btype[i][j] == (DOWN + DOWN))
			{
				btype[i][j] = DOWN;
			}
			else if (btype[i][j] == (RIGHT + RIGHT))
			{
				btype[i][j] = RIGHT;
			}
			else if (btype[i][j] == (UP + UP))
			{
				btype[i][j] = UP;
			}

		}
	}


	_for(k, 0, BND_NUM)
	{
		double dr = boundary_array[k].end.r - boundary_array[k].start.r;
		double dz = boundary_array[k].end.z - boundary_array[k].start.z;

		if (dr == 0 || dz == 0) continue;
		if (dr > dz)
		{
			double ins_z = dz / dr;
			i = 0;
			_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
			{
				
				btype[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j] = judge_conner((int)(boundary_array[k].start.z + ceil(i * ins_z)), j);

				i++;
			}
		}
		else
		{
			double ins_r = dr / dz;
			j = 0;
			_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
			{
				btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))] = judge_conner(i,(int)(boundary_array[k].start.r + ceil(j * ins_r)));
				j++;
			}
		}

	}
	//for (int i = 0; i < nz; i++)
	//{
	//	for (int j = 0; j < nr; j++)
	//	{
	//		btype[i][j] = judge_conner(i, j);
	//	}
	//}

#ifdef BOUNDARY_DEBUG
	matrix_to_csv((double**)world, ZMAX, RMAX, RMAX, (char*)(".\\output\\world.csv"));
	matrix_to_csv((double**)btype, ZMAX, RMAX, RMAX, (char*)(".\\output\\btype.csv"));
	matrix_to_csv((double**)ptype, ZMAX, RMAX, RMAX, (char*)(".\\output\\ptype.csv"));
#endif
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] != 0 && btype[i][j] != 110)
			{
				MPDT[i][j].ne = bg_den / scale;
				MPDT[i][j].ni = bg_den / scale;
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
			//MPDT[i][j].ee = 0;
			//MPDT[i][j].ei = 0;

			MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) + 0.5 * MPDT[i][j].ne * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
			MPDT[i][j].ei = MPDT[i][j].pi / (gamma - 1) + 0.5 * MPDT[i][j].ni * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);

			btheta[i][j] = 0;
		}
	}

	return 0;
}

int fill_plasma(int tr, int tz, int fill_n)
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

//�߽���������
void  boundary_condition()
{

	register int i, j, k;

	//����߽�

	_for(k, 0, BND_NUM)
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

		//�մɱ߽���ε��ӷ���
		if (boundary_array[k].physics_type == DIELECTRIC_SURFACE_BOUNDARY)
		{
			double dr = boundary_array[k].end.r - boundary_array[k].start.r;
			double dz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (dr > dz)
			{
				double ins_z = dz / dr;
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
					i++;
				}
			}
			else
			{
				double ins_r = dr / dz;
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
					j++;
				}
			}
		}

		if (boundary_array[k].physics_type == ANODE_BOUNDARY)
		{
			double dr = boundary_array[k].end.r - boundary_array[k].start.r;
			double dz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (dr > dz)
			{
				double ins_z = dz / dr;
				i = 0;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ver = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vir = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vetheta = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vitheta = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].vez = 0;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].viz = 0;
					i++;
				}
			}
			else
			{
				double ins_r = dr / dz;
				j = 0;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ver = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vir = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vetheta = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vitheta = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vez = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].viz = 0;
					j++;
				}
			}
		}

	}

	//��������͵��ӽ�����뿪��������
	_for(k, 0, BND_NUM)
	{

		if (boundary_array[k].physics_type == CATHODE_BOUNDARY)
		{
			double dr = boundary_array[k].end.r - boundary_array[k].start.r;
			double dz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (dr > dz)
			{
				double ins_z = dz / dr;
				i = 0;
				_for(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					//if (btype[(int)(boundary_array[k].start.z + ceil(i * ins_z)) + 1][j] != 1) continue;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z)) + 1][j].ne += inter_e_den / scale;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ne += inter_e_den / scale;
					i++;
				}
			}
			else
			{
				double ins_r = dr / dz;

				int num = dz / 2;
				double den_step = (inter_e_den - 6.02e5) / num;
				j = 0;

				_for(i, boundary_array[k].start.z, (boundary_array[k].start.z + boundary_array[k].end.z) / 2)
				{
					if (btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))] == DOWN)
					{
						//if (btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) + 1] != 1) continue;
						MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) + 1].ne += bg_den / scale;
						MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne += bg_den / scale;
					}
					else
					{
						//if (btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1] != 1) continue;
						MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1].ne += bg_den / scale;
						MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne += bg_den / scale;
					}
					j++;
				}
			}
		}


		if (boundary_array[k].physics_type == VACCUM_BOUNDARY)
		{
			double dr = boundary_array[k].end.r - boundary_array[k].start.r;
			double dz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (dr > dz)
			{
				double ins_z = dz / dr;
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
				double ins_r = dr / dz;
				j = 0;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					//if (MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne < 1e4) continue;
					//if (btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1] != 1) continue;
					//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1].ne /= 2;
					//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1].ni /= 2;
					//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne /= 2;
					//MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ni /= 2;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ver = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vetheta = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vez = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ni = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vir = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].vitheta = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].viz = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ee = 0;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ei = 0;
					j++;
				}
			}
		}




		//��֤������������������ĵ����ܶ����
		if (boundary_array[k].physics_type == ANODE_BOUNDARY)
		{
			double dr = boundary_array[k].end.r - boundary_array[k].start.r;
			double dz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (dr > dz)
			{
				double ins_z = dz / dr;
				i = 0;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					if (MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ne < 1e4) continue;
					if (btype[(int)(boundary_array[k].start.z + ceil(i * ins_z)) + 1][j] != 1) continue;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z)) + 1][j].ne /= 1.2;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ne /= 1.2;
					i++;
				}
			}
			else
			{
				double ins_r = dr / dz;
				j = 0;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					if (MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne < 1e4) continue;
					if (btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1] != 1) continue;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) - 1].ne /= 1.2;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne /= 1.2;
					j++;
				}
			}
		}

		if (boundary_array[k].physics_type == INLET)
		{
			double dr = boundary_array[k].end.r - boundary_array[k].start.r;
			double dz = boundary_array[k].end.z - boundary_array[k].start.z;

			if (dr > dz)
			{
				double ins_z = dz / dr;
				i = 0;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r - scale)
				{
					//if (btype[(int)(boundary_array[k].start.z + ceil(i * ins_z)) + 1][j] != 1) continue;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z)) + 1][j].ne += inter_pla_den / scale;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z)) + 1][j].ni += inter_pla_den / scale;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ne += inter_pla_den / scale;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ni += inter_pla_den / scale;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z)) + 1][j].ne += inter_e_den / scale;
					MPDT[(int)(boundary_array[k].start.z + ceil(i * ins_z))][j].ne += inter_e_den / scale;
					i++;
				}
			}
			else
			{
				double ins_r = dr / dz;
				j = 0;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z - scale)
				{
					//if (btype[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) + 1] != 1) continue;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) + 1].ne += inter_pla_den / scale;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) + 1].ni += inter_pla_den / scale;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne += inter_pla_den / scale;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ni += inter_pla_den / scale;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r)) + 1].ne += inter_e_den / scale;
					MPDT[i][(int)(boundary_array[k].start.r + ceil(j * ins_r))].ne += inter_e_den / scale;
					j++;
				}
			}
		}

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
		return 1;
	}
	else
	{
		printf("boundary condition error %d %d %d\n",state,i,j);
		return 0;
	}

	return 0;
}