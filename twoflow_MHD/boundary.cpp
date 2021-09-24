#include "MHD.h"


int world[RMAX][ZMAX];
Boundary boundary_array[BND_NUM];

int initial()
{
	int i, j, k;
	//				__________________
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
	boundary_array[0].start.r = 0;
	boundary_array[0].start.z = 20;
	boundary_array[0].end.r = 0;
	boundary_array[0].end.z = 100;
	boundary_array[0].bnd_dir = Z_DIR;
	boundary_array[0].boundary_type = DOWN;

	boundary_array[1].physics_type = VACCUM_BOUNDARY;
	boundary_array[1].start.r = 0;
	boundary_array[1].start.z = 100;
	boundary_array[1].end.r = 80;
	boundary_array[1].end.z = 100;
	boundary_array[1].bnd_dir = R_DIR;
	boundary_array[1].boundary_type = RIGHT;

	boundary_array[2].physics_type = VACCUM_BOUNDARY;
	boundary_array[2].start.r = 80;
	boundary_array[2].start.z = 50;
	boundary_array[2].end.r = 80;
	boundary_array[2].end.z = 100;
	boundary_array[2].bnd_dir = Z_DIR;
	boundary_array[2].boundary_type = UP;

	boundary_array[3].physics_type = ANODE_BOUNDARY;
	boundary_array[3].start.r = 30;
	boundary_array[3].start.z = 50;
	boundary_array[3].end.r = 80;
	boundary_array[3].end.z = 50;
	boundary_array[3].bnd_dir = R_DIR;
	boundary_array[3].boundary_type = LEFT;

	boundary_array[4].physics_type = DIELECTRIC_SURFACE_BOUNDARY;
	boundary_array[4].start.r = 30;
	boundary_array[4].start.z = 0;
	boundary_array[4].end.r = 30;
	boundary_array[4].end.z = 50;
	boundary_array[4].bnd_dir = Z_DIR;
	boundary_array[4].boundary_type = UP;

	boundary_array[5].physics_type = DIELECTRIC_SURFACE_BOUNDARY;
	boundary_array[5].start.r = 2;
	boundary_array[5].start.z = 0;
	boundary_array[5].end.r = 30;
	boundary_array[5].end.z = 0;
	boundary_array[5].bnd_dir = R_DIR;
	boundary_array[5].boundary_type = LEFT;

	boundary_array[6].physics_type = CATHODE_BOUNDARY;
	boundary_array[6].start.r = 2;
	boundary_array[6].start.z = 0;
	boundary_array[6].end.r = 2;
	boundary_array[6].end.z = 20;
	boundary_array[6].bnd_dir = Z_DIR;
	boundary_array[6].boundary_type = DOWN;

	boundary_array[7].physics_type = CATHODE_BOUNDARY;
	boundary_array[7].start.r = 0;
	boundary_array[7].start.z = 20;
	boundary_array[7].end.r = 2;
	boundary_array[7].end.z = 20;
	boundary_array[7].bnd_dir = R_DIR;
	boundary_array[7].boundary_type = RIGHT;


	//≥ı ºªØ ro_boundary £¨

	memset(world, 0, sizeof(world));


	_for(k, 0, BND_NUM)
	{
		if (boundary_array[k].bnd_dir == R_DIR)
		{
			i = boundary_array[k].start.z;
			_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
			{
				world[i][j] += 2;
				//cmprg_gradient[i][j] += boundary_array[k].boundary_type;
			}
		}
		else if (boundary_array[k].bnd_dir == Z_DIR)
		{
			j = boundary_array[k].start.r;
			_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
			{
				world[i][j] += 2;
				//cmprg_gradient[i][j] += boundary_array[k].boundary_type;
			}
		}
	}

	matrix_to_csv((float**)world, ZMAX + 1, RMAX + 1, RMAX + 1, (char*)(".\\output\\world.csv"));
	return 0;
}