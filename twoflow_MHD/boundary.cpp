#include "MHD.h"

using namespace std;

double world[ZMAX][RMAX];
double btype[ZMAX][RMAX];
double ptype[ZMAX][RMAX];
Boundary boundary_array[BND_NUM];

int initial()
{
	int i, j, k;
	//			    __________________
	//				|        50       |
	//				|				  |
	//				|50				  |
	//				|				  |
	//______________|				  |
	//|       50					  |
	//|28							  |
	//|___20___						  |
	//		  |2______________80______|

	boundary_array[0].physics_type = CYLINDRICAL_AXIS;
	boundary_array[0].start.r = 0 * scale;
	boundary_array[0].start.z = 20 * scale;
	boundary_array[0].end.r = 0 * scale;
	boundary_array[0].end.z = 100 * scale;
	boundary_array[0].bnd_dir = Z_DIR;
	boundary_array[0].boundary_type = DOWN;

	boundary_array[1].physics_type = VACCUM_BOUNDARY;
	boundary_array[1].start.r = 0 * scale;
	boundary_array[1].start.z = 100 * scale;
	boundary_array[1].end.r = 80 * scale;
	boundary_array[1].end.z = 100 * scale;
	boundary_array[1].bnd_dir = R_DIR;
	boundary_array[1].boundary_type = RIGHT;

	boundary_array[2].physics_type = ANODE_BOUNDARY;
	boundary_array[2].start.r = 80 * scale;
	boundary_array[2].start.z = 50 * scale;
	boundary_array[2].end.r = 80 * scale;
	boundary_array[2].end.z = 100 * scale;
	boundary_array[2].bnd_dir = Z_DIR;
	boundary_array[2].boundary_type = UP;

	boundary_array[3].physics_type = DIELECTRIC_SURFACE_BOUNDARY;
	boundary_array[3].start.r = 30 * scale;
	boundary_array[3].start.z = 50 * scale;
	boundary_array[3].end.r = 80 * scale;
	boundary_array[3].end.z = 50 * scale;
	boundary_array[3].bnd_dir = R_DIR;
	boundary_array[3].boundary_type = LEFT;

	boundary_array[4].physics_type = DIELECTRIC_SURFACE_BOUNDARY;
	boundary_array[4].start.r = 30 * scale;
	boundary_array[4].start.z = 0 * scale;
	boundary_array[4].end.r = 30 * scale;
	boundary_array[4].end.z = 50 * scale;
	boundary_array[4].bnd_dir = Z_DIR;
	boundary_array[4].boundary_type = UP;

	boundary_array[5].physics_type = DIELECTRIC_SURFACE_BOUNDARY;
	boundary_array[5].start.r = 2 * scale;
	boundary_array[5].start.z = 0 * scale;
	boundary_array[5].end.r = 30 * scale;
	boundary_array[5].end.z = 0 * scale;
	boundary_array[5].bnd_dir = R_DIR;
	boundary_array[5].boundary_type = LEFT;

	boundary_array[6].physics_type = CATHODE_BOUNDARY;
	boundary_array[6].start.r = 2 * scale;
	boundary_array[6].start.z = 0 * scale;
	boundary_array[6].end.r = 2 * scale;
	boundary_array[6].end.z = 20 * scale;
	boundary_array[6].bnd_dir = Z_DIR;
	boundary_array[6].boundary_type = DOWN;

	boundary_array[7].physics_type = INLET;
	boundary_array[7].start.r = 0 * scale;
	boundary_array[7].start.z = 20 * scale;
	boundary_array[7].end.r = 2 * scale;
	boundary_array[7].end.z = 20 * scale;
	boundary_array[7].bnd_dir = R_DIR;
	boundary_array[7].boundary_type = LEFT;


	//初始化 ro_boundary ，

	memset(world, 0, sizeof(world));
	memset(btype, 0, sizeof(world));
	memset(ptype, 0, sizeof(world));

	_for(k, 0, BND_NUM)
	{
		if (boundary_array[k].bnd_dir == R_DIR)
		{
			i = boundary_array[k].start.z;
			_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
			{
				world[i][j] += 2;
				btype[i][j] += boundary_array[k].boundary_type;
				ptype[i][j] += boundary_array[k].physics_type;
			}
		}
		else if (boundary_array[k].bnd_dir == Z_DIR)
		{
			j = boundary_array[k].start.r;
			_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
			{
				world[i][j] += 2;
				btype[i][j] += boundary_array[k].boundary_type;
				ptype[i][j] += boundary_array[k].physics_type;
			}
		}
	}

	fill_plasma(4 * scale, 2 * scale, 1);
	fill_plasma(0,0, 110);
	for (int i = 0; i < nz; i++)
	{
		for (int j = 0; j < nr; j++)
		{
			if (world[i][j] == 1)
			{
				btype[i][j] = 1;
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
			else if (world[i][j] == 110)
			{
				btype[i][j] = 110;
			}
		}
	}
#ifdef DADI_DEBUG
	matrix_to_csv((double**)world, ZMAX, RMAX, RMAX, (char*)(".\\output\\world.csv"));
	matrix_to_csv((double**)btype, ZMAX, RMAX, RMAX, (char*)(".\\output\\btype.csv"));
	matrix_to_csv((double**)ptype, ZMAX, RMAX, RMAX, (char*)(".\\output\\ptype.csv"));
#endif
	for (int i = 0; i < nz; i++)
	{
		for (int j = 0; j < nr; j++)
		{
			if (btype[i][j] != 0)
			{
				MPDT[i][j].ne = 6.02e5 / scale;
				MPDT[i][j].ni = 6.02e5 / scale;
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

//边界条件处理
void  boundary_condition()
{

	int i, j, k;
	double inter_e_den = 3.500e9;
	double inter_pla_den = 1.11e9;
	//固体边界

	_for(k, 0, BND_NUM)
	{
		if (boundary_array[k].physics_type == CATHODE_BOUNDARY)
		{
			if (boundary_array[k].bnd_dir == R_DIR)
			{
				i = boundary_array[k].start.z;
				_feq(j, boundary_array[k].start.r + 1, boundary_array[k].end.r - 1)
				{
					//if (btype[i + 1][j] != 1) continue;

					//double pd = MPDT[i][j].vez * MPDT[i][j].ne + MPDT[i + 1][j].ne * MPDT[i + 1][j].vez;
					
					//MPDT[i + 1][j].ne += MPDT[i][j].ne;
					//MPDT[i + 1][j].vez = -MPDT[i + 1][j].vez;

					//MPDT[i][j].ne = 0;
					//MPDT[i][j].ver = 0;
					//MPDT[i][j].vez = -MPDT[i][j].vez;

				}
			}
			else if (boundary_array[k].bnd_dir == Z_DIR)
			{
				j = boundary_array[k].start.r;
				_feq(i, boundary_array[k].start.z + 1, boundary_array[k].end.z - 1)
				{
					//if (btype[i][j + 1] != 1) continue;
					//double pd = MPDT[i][j].ver * MPDT[i][j].ne + MPDT[i + 1][j].ne * MPDT[i][j + 1].ver;
					//MPDT[i][j + 1].ne += MPDT[i][j].ne;
					//MPDT[i][j + 1].ver = -pd / MPDT[i][j + 1].ne;

					//MPDT[i][j + 1].ver = -MPDT[i][j + 1].ver;
					//MPDT[i][j].ne = 0;
					//MPDT[i][j].ver = -MPDT[i][j].ver;
					//MPDT[i][j].vez = 0;
				}
			}
		}

		//if (boundary_array[k].physics_type == DIELECTRIC_SURFACE_BOUNDARY)
		//{
		//	if (boundary_array[k].bnd_dir == R_DIR)
		//	{
		//		i = boundary_array[k].start.z;
		//		_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
		//		{
		//			//if (btype[i + 1][j] != 1) continue;

		//			//double pd = MPDT[i][j].vez * MPDT[i][j].ne + MPDT[i + 1][j].ne * MPDT[i + 1][j].vez;

		//			//MPDT[i + 1][j].ne += MPDT[i][j].ne;
		//			//MPDT[i + 1][j].vez = -pd / MPDT[i + 1][j].ne;
		//			MPDT[i + 1][j].vez =  -MPDT[i + 1][j].vez;
		//			//MPDT[i][j].ne = 0;
		//			//MPDT[i][j].ver = 0;
		//			MPDT[i][j].vez = -MPDT[i][j].vez;
		//		}
		//	}
		//	else if (boundary_array[k].bnd_dir == Z_DIR)
		//	{
		//		j = boundary_array[k].start.r;
		//		_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
		//		{
		//			//if (btype[i][j - 1] != 1) continue;
		//			//double pd = MPDT[i][j].ver * MPDT[i][j].ne + MPDT[i - 1][j].ne * MPDT[i][j - 1].ver;
		//			//MPDT[i][j - 1].ne += MPDT[i][j].ne;
		//			//MPDT[i][j - 1].ver = -pd / MPDT[i][j - 1].ne;
		//			MPDT[i][j - 1].ver = -MPDT[i][j - 1].ver;

		//			//MPDT[i][j].ne = 0;
		//			MPDT[i][j].ver = -MPDT[i][j].ver;
		//			//MPDT[i][j].vez = 0;
		//		}
		//	}
		//}

		if (boundary_array[k].physics_type == ANODE_BOUNDARY)
		{
			if (boundary_array[k].bnd_dir == R_DIR)
			{
				i = boundary_array[k].start.z;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					//if (btype[i + 1][j] != 1) continue;

					//double pd = MPDT[i][j].vez * MPDT[i][j].ne + MPDT[i + 1][j].ne * MPDT[i + 1][j].vez;

					//MPDT[i + 1][j].ne += MPDT[i][j].ne;
					//MPDT[i + 1][j].vez = -pd / MPDT[i + 1][j].ne;
					MPDT[i + 1][j].vez = -MPDT[i + 1][j].vez;

					//MPDT[i][j].ne = 0;
					//MPDT[i][j].ver = 0;
					MPDT[i][j].vez = -MPDT[i][j].vez;
				
				}
			}
			else if (boundary_array[k].bnd_dir == Z_DIR)
			{
				j = boundary_array[k].start.r;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					//阳极只有r方向
					//if (btype[i][j - 1] != 1) continue;
					//double pd = MPDT[i][j].ver * MPDT[i][j].ne + MPDT[i - 1][j].ne * MPDT[i][j - 1].ver;
					//MPDT[i][j - 1].ne += MPDT[i][j].ne;
					//MPDT[i][j - 1].ver = -pd / MPDT[i][j - 1].ne;

					//MPDT[i][j].ne = 0;
					//MPDT[i][j].ver = 0;
					//MPDT[i][j].vez = 0;
				}
			}
		}
	}

	//等离子体和电子进入和离开仿真区域
	_for(k, 0, BND_NUM)
	{
		if (boundary_array[k].physics_type == CATHODE_BOUNDARY)
		{
			if (boundary_array[k].bnd_dir == R_DIR)
			{
				i = boundary_array[k].start.z;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r - scale)
				{
					//MPDT[i][j].ver = 100;
					//MPDT[i][j].vez = 100;
					if (btype[i + 1][j] != 1) continue;
					MPDT[i + 1][j].ne += inter_e_den / scale;
					MPDT[i][j].ne += inter_e_den / scale;
				}
			}
			else if (boundary_array[k].bnd_dir == Z_DIR)
			{
				j = boundary_array[k].start.r;
				int con = 1;
				int num = (boundary_array[k].end.z - boundary_array[k].start.z) / 2;
				double den_step = (inter_e_den - 6.02e5) / num;

				_feq(i, (boundary_array[k].start.z + boundary_array[k].end.z) / 2, boundary_array[k].end.z - scale)
				{
					//MPDT[i][j].ver = 100;
					//MPDT[i][j].vez = 100;
					if (btype[i][j + 1] != 1) continue;
					MPDT[i][j + 1].ne += (6.02e5 +  con * den_step) / scale;
					MPDT[i][j].ne += (6.02e5 + con * den_step) / scale;
				}
			}
		}

		if (boundary_array[k].physics_type == VACCUM_BOUNDARY)
		{
			if (boundary_array[k].bnd_dir == R_DIR)
			{
				i = boundary_array[k].start.z;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					if (btype[i - 1][j] != 1) continue;
					MPDT[i - 1][j].ne /= 2;
					MPDT[i - 1][j].ni /= 2;
					MPDT[i][j].ne /= 2;
					MPDT[i][j].ni /= 2;
				}
			}
			else if (boundary_array[k].bnd_dir == Z_DIR)
			{
				j = boundary_array[k].start.r;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					if (btype[i][j - 1] != 1) continue;
					MPDT[i][j - 1].ne /= 2;
					MPDT[i][j - 1].ni /= 2;
					MPDT[i][j].ne /= 2;
					MPDT[i][j].ni /= 2;
				}
			}
		}



		//保证阴极进入和阳极流出的电子密度相等
		if (boundary_array[k].physics_type == ANODE_BOUNDARY)
		{
			if (boundary_array[k].bnd_dir == R_DIR)
			{
				i = boundary_array[k].start.z;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					if (btype[i + 1][j] != 1) continue;
					MPDT[i + 1][j].ne /= 1000;
					MPDT[i][j].ne /= 1000;
				}
			}
			else if (boundary_array[k].bnd_dir == Z_DIR)
			{
				j = boundary_array[k].start.r;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					if (btype[i][j - 1] != 1) continue;
					MPDT[i][j - 1].ne /= 1000;
					MPDT[i][j].ne /= 1000;
				}
			}
		}

		if (boundary_array[k].physics_type == INLET)
		{
			if (boundary_array[k].bnd_dir == R_DIR)
			{
				i = boundary_array[k].start.z;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r - scale)
				{
					//MPDT[i][j].vez = 0;
					//MPDT[i][j].viz = 0;
					if (btype[i + 1][j] != 1) continue;
					MPDT[i + 1][j].ne += inter_pla_den / scale;
					MPDT[i + 1][j].ni += inter_pla_den / scale;
					MPDT[i][j].ne += inter_pla_den / scale;
					MPDT[i][j].ni += inter_pla_den / scale;
					MPDT[i + 1][j].ne += inter_e_den / scale;
					MPDT[i][j].ne += inter_e_den / scale;
				}
			}
			else if (boundary_array[k].bnd_dir == Z_DIR)
			{
				j = boundary_array[k].start.r;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z - scale)
				{
					//MPDT[i][j].vez = 0;
					//MPDT[i][j].viz = 0;
					if (btype[i][j + 1] != 1) continue;
					MPDT[i][j + 1].ne += inter_pla_den / scale;
					MPDT[i][j + 1].ni += inter_pla_den / scale;
					MPDT[i][j].ne += inter_pla_den / scale;
					MPDT[i][j].ni += inter_pla_den / scale;
					MPDT[i][j + 1].ne += inter_e_den / scale;
					MPDT[i][j].ne += inter_e_den / scale;
				}
			}
		}
	}
}