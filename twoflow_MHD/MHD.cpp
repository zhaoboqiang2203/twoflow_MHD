#include "MHD.h"
#include "fluid.h"
#include "dadi.h"
#include "testcode.h"
using namespace std;

int nr;
int nz;
double dr, dz, dtheta;
double dt;

int index;

struct node MPDT[ZMAX][RMAX];
struct node MPDT_in[ZMAX][RMAX];
struct _U U[ZMAX][RMAX], U_bar[ZMAX][RMAX], U_bar2[ZMAX][RMAX];
struct _F Fr[ZMAX][RMAX], Fr_bar[ZMAX][RMAX], Fr_bar2[ZMAX][RMAX];
struct _F Fz[ZMAX][RMAX], Fz_bar[ZMAX][RMAX], Fz_bar2[ZMAX][RMAX];
struct _F s[ZMAX][RMAX], s_bar[ZMAX][RMAX], s_bar2[ZMAX][RMAX];

struct _U Uq[ZMAX][RMAX], Uq_bar[ZMAX][RMAX], Uq_bar2[ZMAX][RMAX];
struct _F Fqr[ZMAX][RMAX], Fqr_bar[ZMAX][RMAX], Fqr_bar2[ZMAX][RMAX];
struct _F Fqz[ZMAX][RMAX], Fqz_bar[ZMAX][RMAX], Fqz_bar2[ZMAX][RMAX];
struct _F sq[ZMAX][RMAX], sq_bar[ZMAX][RMAX], sq_bar2[ZMAX][RMAX];

struct _F tau_vis[ZMAX][RMAX];

double taurr[ZMAX][RMAX];
double taurtheta[ZMAX][RMAX];
double taurz[ZMAX][RMAX];
double tautheta2[ZMAX][RMAX];
double tauthetaz[ZMAX][RMAX];
double tauzz[ZMAX][RMAX];

double pniz[ZMAX][RMAX];
double pnir[ZMAX][RMAX];
double tau_ni[ZMAX][RMAX];
struct _pid pid;

double phi[ZMAX][RMAX];
double rho[ZMAX][RMAX];

double phi1[ZMAX][RMAX];
double rou[ZMAX][RMAX];

int scale;

double vez[ZMAX][RMAX];
double ver[ZMAX][RMAX];
double vethera[ZMAX][RMAX];

double viz[ZMAX][RMAX];
double vir[ZMAX][RMAX];
double vithera[ZMAX][RMAX];

double app_Bz[ZMAX][RMAX];
double app_Br[ZMAX][RMAX];

double Jz[ZMAX][RMAX];
double Jr[ZMAX][RMAX];

double res_out[ZMAX][RMAX];
double res_out1[ZMAX][RMAX];

double btheta[ZMAX][RMAX];


double max_phi;
double set_phi;

double phi_sigma;
double cathode_I;
double current_I;
double last_I;
double REL_MASS;             //相对原子质量
double MI;		            // kg, electron mass
double EPS_PLA;				//相对电导率
double para_p, para_i, para_d;

int main()
{
	int nq;
	double last_current = 0;
	nz = ZMAX;
	nr = RMAX;


	//仿真参数定义区

	scale = ZMAX / 500;
	parameter_read();
	max_phi = 0;
	//清空文件
	ofstream myfile(".\\output\\PHI_I_data.csv");
	myfile.close();

	exit(0);
	dr = 0.001 / scale;
	dz = 0.001 / scale;
	dtheta = PI / 180;
	dt = dr / 1e5 / 3;
	nq = 100;
	//根据背景压强，通气流量，电流密度计算
	bg_den = 1e-3 / (K * 300);
	inter_e_den = cathode_I * 1.1e19;
	inter_pla_den = 0.04  / 40 * NA / 360 / 20 * 1e9;



	//dt = 0.05 * ((dr * dr) + (dz * dz));
	printf("dt = %e\n", dt);
	printf("inter_e_den = %e\n", inter_e_den);
	printf("inter_pla_den = %e\n", inter_pla_den);

	initial();
	magnetic_field_initial();

	//dadi initial
	init_solve();

	//判断是否存在已经仿真好的数据
	if (is_read_datfile())
	{
		read_datfile();
		potential_solve();
	}
	
	while (index--)
	{
		//printf("index %d\n", index);

		//clock_t start = clock();
		boundary_condition();
		electron_flow_v2();

		if (index % nq == 0)
		{
			printf("index %d\n", index);
			max_write();
			potential_solve();
			current_caulate();
			if (abs(current_I - last_current) < 0.3 * abs(current_I))
			{
				current_control();
			}
			last_current = current_I;
			//current_control();
		}


		move();
		//current_caulate();
		//current_control();
		//mag_phi();

		out_judge();

		if (index % 1000 == 0)
		{
			wirte_datfile();
		}

		if (index % 1000 == 0)
			//if(index < 36000)
		{
			output();

		}

		//clock_t ends = clock();
		//cout << "Running Time : " << (double)(ends - start) / CLOCKS_PER_SEC << endl;
	}

	return 0;
}

void matrix_int_to_csv(int** a, int N, int M, int array_size, char* filename)
{
	fstream myfile(filename, ios::out);
	if (!myfile.is_open())
	{
		cout << "未成功打开文件" << endl;
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			//cout << *((double*)a + (array_size * i + j)) << endl;
			myfile << *((int*)a + (array_size * i + j)) << ",";
		}
		myfile << endl;
	}
	myfile.close();
}


void matrix_to_csv(double** a, int N, int M, int array_size, char* filename)
{
	fstream myfile(filename, ios::out);
	if (!myfile.is_open())
	{
		cout << "未成功打开文件" << endl;
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			//cout << *((double*)a + (array_size * i + j)) << endl;
			myfile << fixed << setprecision(10) << *((double*)a + (array_size * i + j)) << ",";
		}
		myfile << endl;
	}
	myfile.close();
}

void matrix_to_binary(char* a, unsigned long long array_size, char* filename)
{
	ofstream fout;
	fout.open(filename, ios::out | ios::binary);

	if (!fout.is_open())
	{
		cout << "未成功打开文件" << endl;
	}

	fout.write((char*)a, array_size);
	fout.close();
}

void output()
{
	// 输出电子密度
	char fname[100];
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ne;
		}
	}
	sprintf_s(fname, (".\\output\\electron density\\electron_density_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//输出电子速度

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ver;
		}
	}
	sprintf_s(fname, (".\\output\\electron ver\\electron_ver_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vez;
		}
	}
	sprintf_s(fname, (".\\output\\electron vez\\electron_vez_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vetheta;
		}
	}

	sprintf_s(fname, (".\\output\\electron vetheta\\electron_vetheta_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


	//输出离子密度
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ni;
		}
	}
	sprintf_s(fname, (".\\output\\ion density\\ion_density_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	// 输出离子速度
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vir;
		}
	}

	sprintf_s(fname, (".\\output\\ion vir\\ion_vir_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].viz;
		}
	}

	sprintf_s(fname, (".\\output\\ion viz\\ion_viz_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].vitheta;
		}
	}

	sprintf_s(fname, (".\\output\\ion vitheta\\ion_vitheta_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	// 输出电势分布

	sprintf_s(fname, (".\\output\\phi\\phi_%d.csv"), index);
	matrix_to_csv((double**)phi, ZMAX, RMAX, RMAX, fname);

	//sprintf_s(fname, (".\\output\\phi\\phi1_%d.csv"), index);
	//matrix_to_csv((double**)phi1, ZMAX, RMAX, RMAX, fname);

	sprintf_s(fname, (".\\output\\rho\\rho_%d.csv"), index);
	matrix_to_csv((double**)rou, ZMAX, RMAX, RMAX, fname);
	// 输出磁场分布

	//输出电子压强分布
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].pe;
		}
	}

	sprintf_s(fname, (".\\output\\electron pe\\electron_pe_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//输出离子压强分布
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].pi;
		}
	}

	sprintf_s(fname, (".\\output\\ion pi\\ion_pi_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	// 输出电子能量分布
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ee * ME/ QE;
		}
	}

	sprintf_s(fname, (".\\output\\electron ee\\electron_ee_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	// 输出离子能量分布
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = MPDT[i][j].ei * MI / QE;
		}
	}

	sprintf_s(fname, (".\\output\\ion ei\\ion_ei_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


	//Ez分布
	sprintf_s(fname, (".\\output\\Ez\\Ez_%d.csv"), index);
	matrix_to_csv((double**)Ez, ZMAX, RMAX, RMAX, fname);

	//Er分布
	sprintf_s(fname, (".\\output\\Er\\Er_%d.csv"), index);
	matrix_to_csv((double**)Er, ZMAX, RMAX, RMAX, fname);

	// 输出感生磁场分布
	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].br;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\br\\br_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].btheta;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\btheta\\btheta_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].bz;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\bz\\bz_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	 
	//sprintf_s(fname, (".\\output\\btheta\\btheta_%d.csv"), index);
	//matrix_to_csv((double**)btheta, ZMAX, RMAX, RMAX, fname);

	////电子离子碰撞频率
	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].mu_ie;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\mu_ie\\mu_ie_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
	
	//电流密度r方向
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = (MPDT[i][j].ni * MPDT[i][j].vir - MPDT[i][j].ne * MPDT[i][j].ver) * QE;
		}
	}

	sprintf_s(fname, (".\\output\\Jr\\Jr_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//电流密度z方向
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			res_out[i][j] = (MPDT[i][j].ni * MPDT[i][j].viz - MPDT[i][j].ne * MPDT[i][j].vez) * QE;
		}
	}

	sprintf_s(fname, (".\\output\\Jz\\Jz_%d.csv"), index);
	matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	////电子向离子转移能量
	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].delta_ei;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\delta_ei\\delta_ei_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	////电子离子碰撞截面
	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].sigma_Q;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\sigma\\sigma_Q_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	////多余电荷密度
	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].peq;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\peq\\peq_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].neq;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\neq\\neq_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].vnqr;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\vnqr\\vnqr_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].vnqz;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\vnqz\\vnqz_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].vnqtheta;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\vnqtheta\\vnqtheta_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);



	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].vpqr;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\vpqr\\vpqr_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].vpqz;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\vpqz\\vpqz_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].vpqtheta;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\vpqtheta\\vpqtheta_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);



	//for (int i = 0; i < ZMAX; i++)
	//{
	//	for (int j = 0; j < RMAX; j++)
	//	{
	//		res_out[i][j] = MPDT[i][j].angle_b_vi;
	//	}
	//}

	//sprintf_s(fname, (".\\output\\angle\\angle_%d.csv"), index);
	//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


}

void out_judge()
{
	int i, j;
	int err_flag = 0;
	for (i = 0; i < ZMAX; i++)
	{
		for (j = 0; j < RMAX; j++)
		{
			if (btype[i][j] != 1) continue;
			if (MPDT[i][j].ne != MPDT[i][j].ne || MPDT[i][j].ne < 0)
			{
				err_flag = 1;
				break;
			}

			if (MPDT[i][j].ni != MPDT[i][j].ni || MPDT[i][j].ni < 0)
			{
				err_flag = 1;
				break;
			}

			if (MPDT[i][j].ver != MPDT[i][j].ver)
			{
				err_flag = 1;
				break;
			}

			if (MPDT[i][j].vetheta != MPDT[i][j].vetheta)
			{
				err_flag = 1;
				break;
			}
			if (MPDT[i][j].vez != MPDT[i][j].vez)
			{
				err_flag = 1;
				break;
			}

			if (MPDT[i][j].vir != MPDT[i][j].vir)
			{
				err_flag = 1;
				break;
			}

			if (MPDT[i][j].vitheta != MPDT[i][j].vitheta)
			{
				err_flag = 1;
				break;
			}
			if (MPDT[i][j].viz != MPDT[i][j].viz)
			{
				err_flag = 1;
				break;
			}

			if (MPDT[i][j].ee != MPDT[i][j].ee)
			{
				err_flag = 1;
				break;
			}

			if (MPDT[i][j].ei != MPDT[i][j].ei)
			{
				err_flag = 1;
				break;
			}
		}

		if (err_flag == 1)
		{
			break;
		}
	}

	if (err_flag == 1)
	{
		output();
		printf("i = %d j = %d\n", i, j);
	}
}


void read_csv()
{
	char fname[100];
	string ss;
	register int i, j;
	fstream infile;

	// 读取电子密度
	sprintf_s(fname, (".\\output\\electron density\\electron_density_%d.csv"), index);

	infile.open(fname, ios::in);
	if (!infile.is_open())
	{
		cout << "未成功打开文件" << endl;
	}

	i = 0;
	while (getline(infile, ss))
	{
		string number;
		istringstream readstr(ss); //string数据流化
		//将一行数据按'，'分割
		for (j = 0; j < RMAX; j++)
		{
			getline(readstr, number, ','); //循环读取数据
			res_out[i][j] = atof(number.c_str());
		}
		i++;
	}
	//printf("ne i = %d\n", i);
	infile.close();
	for (i = 0; i < ZMAX; i++)
	{
		for (j = 0; j < RMAX; j++)
		{
			if (MPDT[i][j].ne != res_out[i][j])
			{
				printf("MPDT[%d][%d].ne = %lf, read = %lf\n", i, j, MPDT[i][j].ne, res_out[i][j]);
			}
		}
	}

	// 读取电子速度
	sprintf_s(fname, (".\\output\\electron ver\\electron_ver_%d.csv"), index);

	infile.open(fname, ios::in);
	if (!infile.is_open())
	{
		cout << "未成功打开文件" << endl;
	}

	i = 0;
	while (getline(infile, ss))
	{
		string number;
		istringstream readstr(ss); //string数据流化
		//将一行数据按'，'分割
		for (j = 0; j < RMAX; j++)
		{
			getline(readstr, number, ','); //循环读取数据
			res_out[i][j] = atof(number.c_str());
		}
		i++;
	}
	infile.close();
	for (i = 0; i < ZMAX; i++)
	{
		for (j = 0; j < RMAX; j++)
		{
			if (abs(MPDT[i][j].ver - res_out[i][j]) > 1e6)
			{
				printf("MPDT[%d][%d].ver = %lf, read = %lf\n", i, j, MPDT[i][j].ver, res_out[i][j]);
			}
		}
	}

		////输出电子速度

		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].ver;
		//	}
		//}
		//sprintf_s(fname, (".\\output\\electron ver\\electron_ver_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].vez;
		//	}
		//}
		//sprintf_s(fname, (".\\output\\electron vez\\electron_vez_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].vetheta;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\electron vetheta\\electron_vetheta_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


		////输出离子密度
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].ni;
		//	}
		//}
		//sprintf_s(fname, (".\\output\\ion density\\ion_density_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//// 输出离子速度
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].vir;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\ion vir\\ion_vir_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].viz;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\ion viz\\ion_viz_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].vitheta;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\ion vitheta\\ion_vitheta_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//// 输出电势分布

		//sprintf_s(fname, (".\\output\\phi\\phi_%d.csv"), index);
		//matrix_to_csv((double**)phi, ZMAX, RMAX, RMAX, fname);

		////sprintf_s(fname, (".\\output\\phi\\phi1_%d.csv"), index);
		////matrix_to_csv((double**)phi1, ZMAX, RMAX, RMAX, fname);

		//sprintf_s(fname, (".\\output\\rho\\rho_%d.csv"), index);
		//matrix_to_csv((double**)rou, ZMAX, RMAX, RMAX, fname);
		//// 输出磁场分布

		////输出电子压强分布
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].pe;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\electron pe\\electron_pe_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		////输出离子压强分布
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].pi;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\ion pi\\ion_pi_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
		//// 输出电子能量分布
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].ee / QE;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\electron ee\\electron_ee_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);
		//// 输出离子能量分布
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].ei / QE;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\ion ei\\ion_ei_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


		////Ez分布
		//sprintf_s(fname, (".\\output\\Ez\\Ez_%d.csv"), index);
		//matrix_to_csv((double**)Ez, ZMAX, RMAX, RMAX, fname);

		////Er分布
		//sprintf_s(fname, (".\\output\\Er\\Er_%d.csv"), index);
		//matrix_to_csv((double**)Er, ZMAX, RMAX, RMAX, fname);

		//// 输出感生磁场分布
		////for (int i = 0; i < ZMAX; i++)
		////{
		////	for (int j = 0; j < RMAX; j++)
		////	{
		////		res_out[i][j] = MPDT[i][j].br;
		////	}
		////}

		////sprintf_s(fname, (".\\output\\br\\br_%d.csv"), index);
		////matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		////for (int i = 0; i < ZMAX; i++)
		////{
		////	for (int j = 0; j < RMAX; j++)
		////	{
		////		res_out[i][j] = MPDT[i][j].btheta;
		////	}
		////}

		////sprintf_s(fname, (".\\output\\btheta\\btheta_%d.csv"), index);
		////matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		////for (int i = 0; i < ZMAX; i++)
		////{
		////	for (int j = 0; j < RMAX; j++)
		////	{
		////		res_out[i][j] = MPDT[i][j].bz;
		////	}
		////}

		////sprintf_s(fname, (".\\output\\bz\\bz_%d.csv"), index);
		////matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//sprintf_s(fname, (".\\output\\btheta\\btheta_%d.csv"), index);
		//matrix_to_csv((double**)btheta, ZMAX, RMAX, RMAX, fname);

		////电子离子碰撞频率
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].mu_ie;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\mu_ie\\mu_ie_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		////电流密度r方向
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = (MPDT[i][j].ni * MPDT[i][j].vir - MPDT[i][j].ne * MPDT[i][j].ver) * QE / dt;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\Jr\\Jr_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		////电流密度z方向
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = (MPDT[i][j].ni * MPDT[i][j].viz - MPDT[i][j].ne * MPDT[i][j].vez) * QE / dt;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\Jz\\Jz_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		////电子向离子转移能量
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].delta_ei;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\delta_ei\\delta_ei_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		////电子离子碰撞截面
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].sigma_Q;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\sigma\\sigma_Q_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		////多余电荷密度
		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].peq;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\peq\\peq_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].neq;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\neq\\neq_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].vnqr;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\vnqr\\vnqr_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].vnqz;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\vnqz\\vnqz_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].vnqtheta;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\vnqtheta\\vnqtheta_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);



		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].vpqr;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\vpqr\\vpqr_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].vpqz;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\vpqz\\vpqz_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);

		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].vpqtheta;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\vpqtheta\\vpqtheta_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);



		//for (int i = 0; i < ZMAX; i++)
		//{
		//	for (int j = 0; j < RMAX; j++)
		//	{
		//		res_out[i][j] = MPDT[i][j].angle_b_vi;
		//	}
		//}

		//sprintf_s(fname, (".\\output\\angle\\angle_%d.csv"), index);
		//matrix_to_csv((double**)res_out, ZMAX, RMAX, RMAX, fname);


	}

void wirte_datfile()
{
	char fname[100];
	ofstream fout;
	sprintf_s(fname, (".\\output\\MPDT\\MPDT_%d.dat"), index);
	fout.open(fname, ofstream::binary);
	if (!fout.is_open())
	{
		cout << fname << "文件未能打开" << endl;
	}
	fout.write((char*)&MPDT, sizeof(MPDT));
	fout.close();
}

void read_datfile()
{
	char fname[100];
	ifstream fin;
	sprintf_s(fname, (".\\output\\MPDT\\MPDT_%d.dat"), index);

	fin.open(fname, ifstream::binary);
	if (!fin.is_open())
	{
		cout << fname << "文件未能打开" << endl;
		return;
	}
	fin.read((char*)&MPDT, sizeof(MPDT));
	fin.close();

}

bool is_read_datfile()
{
	char fname[100];
	ifstream fin;
	sprintf_s(fname, (".\\output\\MPDT\\MPDT_%d.dat"), index);

	fin.open(fname, ifstream::binary);
	if (!fin.is_open())
	{
		cout << fname << "文件未能打开" << endl;
		return false;
	}
	fin.close();

	return true;
}

void judge_bit()
{
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			if (MPDT_in[i][j].ne != MPDT[i][j].ne)
			{
				printf("MPDT[%d][%d].ne = %lf\n", i, j, MPDT[i][j].ne);
			}

			if (MPDT_in[i][j].ver != MPDT[i][j].ver)
			{
				printf("MPDT[%d][%d].ver = %lf\n", i, j, MPDT[i][j].ver);
			}

			if (MPDT_in[i][j].vez != MPDT[i][j].vez)
			{
				printf("MPDT[%d][%d].vez = %lf\n", i, j, MPDT[i][j].vez);
			}

		}
	}
}

void max_write()
{
	fstream myfile(".\\output\\PHI_I_data.csv", ios::app);
	if (!myfile.is_open())
	{
		cout << "未成功打开文件" << endl;
	}

	myfile << index << "," << max_phi << "," << current_I << endl;

	myfile.close();
}
