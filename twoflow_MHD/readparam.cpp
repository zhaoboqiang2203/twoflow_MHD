#include "json/json.h"
#include <iostream>
#include <memory>
#include "MHD.h"

using namespace std;

int parameter_read()
{
	//Json::Value root;
	//Json::Value magnet;
	//Json::StreamWriterBuilder builder;
	//const std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());

	//root["Magnetic Field"] = 11;
	//root["Boundary"] = 20;
	////magnet["length"] = 10;
	//writer->write(root, &cout);

	Json::Value root;
	Json::Value bdyC;
	Json::Value mgC;
	double cur_para;
	double mag_field;

	std::ifstream ifs;
	ifs.open("mpdtParam.json");

	Json::CharReaderBuilder builder;
	builder["collectComments"] = true;
	JSONCPP_STRING errs;
	if (!parseFromStream(builder, ifs, &root, &errs)) {
		std::cout << errs << std::endl;
		return EXIT_FAILURE;
	}



	index = root["index"].asInt();
	cathode_I = root["cathod current"].asInt();
	REL_MASS = root["particle"][0]["REL_MASS"].asDouble();             //相对原子质量
	Ionization_Energy = root["particle"][0]["IONIZATION_ENERGY"].asDouble();           //原子电离能
	EPS_PLA = root["EPS_PLA"].asDouble();
	MI = REL_MASS * AMU;		        // kg, electron mass

	is_atom_sim = root["ATOM_SIM"].asInt();
	max_speed = root["max speed"].asDouble();

#ifdef BOUNDARY_DEBUG
	std::cout << root << std::endl;
#endif
	
	printf("index = %d\n", index);
	printf("cathode_I = %lf\n", cathode_I);
	printf("EPS_PLA = %e\n", EPS_PLA);
	printf("REL_MASS = %lf\n", REL_MASS);

	bnd_size = root["boundary"].size();
	for (int i = 0; i < root["boundary"].size(); i++)
	{
		bdyC = root["boundary"][i];

		boundary_array[i].physics_type = ptype_trans(bdyC["physics_type"].asString());
		boundary_array[i].start.r = bdyC["start.r"].asDouble() * scale;
		boundary_array[i].start.z = bdyC["start.z"].asDouble() *scale;
		boundary_array[i].end.r = bdyC["end.r"].asDouble() * scale;
		boundary_array[i].end.z = bdyC["end.z"].asDouble() * scale;
		boundary_array[i].bnd_dir = dir_trans(bdyC["bnd_dir"].asString());
		boundary_array[i].boundary_type = btype_trans(bdyC["boundary_type"].asString());
		//cout << "start.r = " << boundary_array[i].start.r << endl;
		//cout << "start.z = " << boundary_array[i].start.z << endl;
		//cout << "end.r = " << boundary_array[i].end.r << endl;
		//cout << "end.z = " << boundary_array[i].end.z << endl;
	}

	mgC = root["magnetic coil"][0];

	//coil_I = mgC["current"].asDouble();//线圈电流
	coil_R = mgC["radius"].asDouble();//线圈内径
	coil_W = mgC["weith"].asDouble();//线圈宽度（内外径差值）
	coil_L = mgC["length"].asDouble();//线圈厚度

	cur_para = mgC["currenr para"].asDouble();
	mag_field = mgC["magnetic field"].asDouble();
	
	coil_I = mag_field / cur_para;  //线圈电流

	printf("coil current = %lf\n", coil_I);
	printf("length = %lf\n", coil_L);
	printf("weith = %lf\n", coil_W);
	printf("radius = %lf\n", coil_R);

	pid.set_current = cathode_I;
	pid.Kp = root["Kp"].asDouble();
	pid.Ki = root["Ki"].asDouble();
	pid.Kd = root["Kd"].asDouble();
	pid.integral = 0;


	ion_area_rmin = root["ionization area"][0]["rmin"].asInt();
	ion_area_rmax = root["ionization area"][0]["rmax"].asInt();
	ion_area_zmin = root["ionization area"][0]["zmin"].asInt();
	ion_area_zmax = root["ionization area"][0]["zmax"].asInt();

	cur_pos_z = root["current surface"][0]["pos_z"].asInt();
	cur_pos_r = root["current surface"][0]["pos_r"].asInt();

	printf("cur_pos_r = %d\n", cur_pos_r);
	printf("cur_pos_z = %d\n", cur_pos_z);
	printf("ion_area_rmin = %d\n", ion_area_rmin);
	printf("ion_area_rmax = %d\n", ion_area_rmax);
	printf("ion_area_zmin = %d\n", ion_area_zmin);
	printf("ion_area_zmax = %d\n", ion_area_zmax);

	ifs.close();
	return EXIT_SUCCESS;

}


PTypes ptype_trans(string ss)
{
		if (ss == "VACCUM_BOUNDARY")
		{
			return VACCUM_BOUNDARY;
		}
		else if(ss =="DIELECTRIC_SURFACE_BOUNDARY")
		{
			return DIELECTRIC_SURFACE_BOUNDARY;
		}
		else if (ss=="EXTERN_INTERIOR_BOUNDARY")
		{
			return EXTERN_INTERIOR_BOUNDARY;

		}
		else if (ss == "PERIODIC_BOUNDARY")
		{
			return PERIODIC_BOUNDARY;
		}
		else if (ss == "INLET")
		{
			return INLET;
		}
		else if (ss == "MIRROR_REFLECTION_BOUNDARY")
		{
			return MIRROR_REFLECTION_BOUNDARY;
		}
		else if (ss == "ANODE_BOUNDARY")
		{
			return ANODE_BOUNDARY;
		}
		else if (ss == "CATHODE_BOUNDARY")
		{
			return CATHODE_BOUNDARY;
		}
		else if (ss == "CYLINDRICAL_AXIS")
		{
			return CYLINDRICAL_AXIS;
		}
		else if (ss == "CONDUCTING_BOUNDARY")
		{
			return CONDUCTING_BOUNDARY;
		}

}

int dir_trans(string ss)
{
	if (ss == "R_DIR")
	{
		return R_DIR;
	}
	else if (ss == "Z_DIR")
	{
		return Z_DIR;
	}
	return -1;
}

int btype_trans(string ss)
{
	if (ss == "LEFT")
	{
		return LEFT;
	}
	else if (ss == "UP")
	{
		return UP;
	}
	else if (ss == "RIGHT")
	{
		return RIGHT;
	}
	else if (ss == "DOWN")
	{
		return DOWN;
	}
	return -1;
}
