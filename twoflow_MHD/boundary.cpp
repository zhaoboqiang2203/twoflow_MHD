#include "MHD.h"

using namespace std;

float world[RMAX][ZMAX];
float btype[RMAX][ZMAX];
float tb[RMAX][ZMAX];
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

	boundary_array[2].physics_type = VACCUM_BOUNDARY;
	boundary_array[2].start.r = 80 * scale;
	boundary_array[2].start.z = 50 * scale;
	boundary_array[2].end.r = 80 * scale;
	boundary_array[2].end.z = 100 * scale;
	boundary_array[2].bnd_dir = Z_DIR;
	boundary_array[2].boundary_type = UP;

	boundary_array[3].physics_type = ANODE_BOUNDARY;
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
	//memset(world, 0, sizeof(world));

	_for(k, 0, BND_NUM)
	{
		if (boundary_array[k].bnd_dir == R_DIR)
		{
			i = boundary_array[k].start.z;
			_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
			{
				world[i][j] += 2;
				btype[i][j] += boundary_array[k].boundary_type;
			}
		}
		else if (boundary_array[k].bnd_dir == Z_DIR)
		{
			j = boundary_array[k].start.r;
			_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
			{
				world[i][j] += 2;
				btype[i][j] += boundary_array[k].boundary_type;
			}
		}
	}

	fill_plasma(4 * scale, 2 * scale, 1);

	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nz; j++)
		{
			if (world[i][j] == 1)
			{
				btype[i][j] = 1;
			}
		}
	}
	matrix_to_csv((float**)world, ZMAX, RMAX, RMAX, (char*)(".\\output\\world.csv"));
	matrix_to_csv((float**)btype, ZMAX, RMAX, RMAX, (char*)(".\\output\\btype.csv"));
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nz; j++)
		{
			MPDT[i][j].ne = 6.02e10;
			MPDT[i][j].ni = 6.02e10;
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
			MPDT[i][j].ee = 0;
			MPDT[i][j].ei = 0;
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
	_for(k, 0, BND_NUM)
	{
		if (boundary_array[k].physics_type == CATHODE_BOUNDARY)
		{
			if (boundary_array[k].bnd_dir == R_DIR)
			{
				i = boundary_array[k].start.z;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					MPDT[i + 1][j].ne += 7.4e17;
				}
			}
			else if (boundary_array[k].bnd_dir == Z_DIR)
			{
				j = boundary_array[k].start.r;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					MPDT[i][j + 1].ne += 7.4e17;
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
					MPDT[i - 1][j].ne /= 2;
					MPDT[i - 1][j].ni /= 2;
				}
			}
			else if (boundary_array[k].bnd_dir == Z_DIR)
			{
				j = boundary_array[k].start.r;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					MPDT[i][j - 1].ne /= 2;
					MPDT[i][j - 1].ni /= 2;
				}
			}
		}

		//if(boundary_array[k].physics_type == VACCUM_BOUNDARY)
		//{
		//	if (boundary_array[k].bnd_dir == R_DIR)
		//	{
		//		i = boundary_array[k].start.z;
		//		_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
		//		{
		//			MPDT[i][j].ne /= 2;
		//			MPDT[i][j].ni /= 2;
		//		}
		//	}
		//	else if (boundary_array[k].bnd_dir == Z_DIR)
		//	{
		//		j = boundary_array[k].start.r;
		//		_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
		//		{
		//			MPDT[i][j].ne /= 2;
		//			MPDT[i][j].ni /= 2;
		//		}
		//	}
		//}

		//保证阴极进入和阳极流出的电子密度相等
		if (boundary_array[k].physics_type == ANODE_BOUNDARY)
		{
			if (boundary_array[k].bnd_dir == R_DIR)
			{
				i = boundary_array[k].start.z;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					MPDT[i + 1][j].ne /= 1000;

				}
			}
			else if (boundary_array[k].bnd_dir == Z_DIR)
			{
				j = boundary_array[k].start.r;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					MPDT[i][j - 1].ne /= 1000;
				}
			}
		}

		if (boundary_array[k].physics_type == INLET)
		{
			if (boundary_array[k].bnd_dir == R_DIR)
			{
				i = boundary_array[k].start.z;
				_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
				{
					MPDT[i + 1][j].ne += 6.02e12;
					MPDT[i + 1][j].ni += 6.02e12;
				}
			}
			else if (boundary_array[k].bnd_dir == Z_DIR)
			{
				j = boundary_array[k].start.r;
				_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
				{
					MPDT[i][j + 1].ne += 6.02e12;
					MPDT[i][j + 1].ni += 6.02e12;
				}
			}
		}
	}
}