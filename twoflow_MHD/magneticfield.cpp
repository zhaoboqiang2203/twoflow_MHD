#include "MHD.h"

double c0 = 2e-7;// mu0 / 2pi = 2e-7
double I = 18000;// current, unit: A
double a = 0.5;// radius of current coil, unit: m
double eps = 1e-4;

double coil_I;//线圈电流
double coil_R;//线圈内径
double coil_W;//线圈宽度（内外径差值）
double coil_L;//线圈厚度

// Set field for single coil
double rho1(double x, double y, double z)
{
    return sqrt(x * x + y * y);
}
double theta(double x, double y, double z)
{
    return atan2(y, x);
}
double k2_0(double x, double y, double z)
{
   return 4 * a * rho1(x, y, z)/ (sqr(a + rho1(x, y, z)) + z * z);
}
double k2(double x, double y, double z)
{
    if (k2_0(x, y, z) == 1)
    {
        return 1 - eps;
    }
    else
    {
        return k2_0(x, y, z);
    }
}

double K1(double x, double y, double z)
{
    return PI / 2 * (1 + 0.25 * k2(x, y, z) + sqr(0.375) * sqr(k2(x, y, z)) + sqr(0.3125) * cube(k2(x, y, z)));
}
double E(double x, double y, double z)
{
    double k2_rz = k2(x, y, z);
    return PI / 2 * (1 - 0.25 * k2_rz + sqr(0.375) * sqr(k2_rz)/ 3.0 - sqr(0.3125) * cube(k2_rz)/ 5.0);
}
double a1(double x, double y, double z)
{
    return sqrt(sqr(a + rho1(x, y, z)) + sqr(z));
}
double b1_0(double x, double y, double z)
{
    return sqr(a - rho1(x, y, z)) + sqr(z);
}
double b1(double x, double y, double z)
{
    if (b1_0(x, y, z) == 0)
    {
        return eps;
    }
    else
    {
        return b1_0(x, y, z);
    }
}
double b2(double x, double y, double z)
{
    return a * a + sqr(rho1(x, y, z)) + z * z;
}
double b3(double x, double y, double z)
{
    return a * a  - sqr(rho1(x, y, z)) - z * z;
}
double Brho(double x, double y, double z)
{
    if (rho1(x, y, z) == 0)
    {
        return c0 * I * z / (eps * a1(x, y, z)) * (b2(x, y, z) / b1(x, y, z) * E(x, y, z) - K1(x, y, z));
    }
    else
    {
        return c0 * I * z / (rho1(x, y, z) * a1(x, y, z)) * (b2(x, y, z) / b1(x, y, z) * E(x, y, z) - K1(x, y, z));
    }
    
}

double Bscz(double x, double y, double z)
{
    return c0 * I / a1(x, y, z) * (b3(x, y, z) / b1(x, y, z) * E(x, y, z) + K1(x, y, z));
}

double Btheta(double x, double y, double z)
{
    return 0;
}

//void magnetic_field_initial()
//{
//    double bias;
//    a = 0.1;
//    bias = 0.01;
//    for (int i = 0; i < nz; i++)
//    {
//        for (int j = 0; j < nr; j++)
//        {
//            app_Br[i][j] = Brho(j * dr, 0, i * dz - bias);
//            app_Bz[i][j] = Bscz(j * dr, 0, i * dz - bias);
//        }
//    }
//
//    matrix_to_csv((double**)app_Br, ZMAX, RMAX, RMAX, (char*)(".\\output\\app_Br.csv"));
//    matrix_to_csv((double**)app_Bz, ZMAX, RMAX, RMAX, (char*)(".\\output\\app_Bz.csv"));
//
//   
//}

void magnetic_field_initial()
{
    double length;
    double weith;

    //length = 0.076;
    //weith = 0.001;
    length = coil_L;
    weith = coil_W;

    double nd = 0.001;
    double wth = 0;
    double lth = 0;
    double t0;
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nr; j++)
        {
            app_Br[i][j] = 0;
            app_Bz[i][j] = 0;
        }
    }

    I = coil_I / (length / nd * weith / nd);

    while (wth < weith)
    {
        a = coil_R + wth;
        while (lth < length)
        {
            for (int i = 0; i < nz; i++)
            {
                for (int j = 0; j < nr; j++)
                {
                    app_Br[i][j] += Brho(j * dr, 0, i * dz - lth);
                    app_Bz[i][j] += Bscz(j * dr, 0, i * dz - lth);

                    t0 = 2 * PI * ME / (QE * sqrt(sqr(app_Br[i][j]) + sqr(app_Bz[i][j])));


                    
                    e_half[i][i] = -fmod(dt, t0) * QE / ME / 2;
                    e_hrho[i][i] = e_half[i][i] * app_Br[i][j];
                    e_htheta[i][i] = 0;
                    e_hz[i][i] = e_half[i][i] * app_Bz[i][j];
                    e_h2[i][i] = sqr(e_hrho[i][i]) + sqr(e_hz[i][i]);
                    
                    e_srho[i][i] = 2 * e_hrho[i][i] / (1 + e_h2[i][i]);
                    e_stheta[i][i] = 2 * e_htheta[i][i] / (1 + e_h2[i][i]);
                    e_sz[i][i] = 2 * e_hz[i][i] / (1 + e_h2[i][i]);

                    t0 = 2 * PI * MI / (QE * sqrt(sqr(app_Br[i][j]) + sqr(app_Bz[i][j])));

                    i_half[i][j] = fmod(dt, t0) * QE / MI / 2;
                    i_hrho[i][j] = i_half[i][j] * app_Br[i][j];
                    i_htheta[i][j] = 0;
                    i_hz[i][j] = i_half[i][j] * app_Bz[i][j];
                    i_h2[i][j] = sqr(i_hrho[i][j]) + sqr(i_hz[i][j]);

                    i_srho[i][j] = 2 * i_hrho[i][j] / (1 + i_h2[i][j]);
                    i_stheta[i][j] = 2 * i_htheta[i][j] / (1 + i_h2[i][j]);
                    i_sz[i][j] = 2 * i_hz[i][j] / (1 + i_h2[i][j]);
                }
            }

            lth += nd;
        }

        wth += nd;
    }

    matrix_to_csv((double**)app_Br, ZMAX, RMAX, RMAX, (char*)(".\\output\\app_Br.csv"));
    matrix_to_csv((double**)app_Bz, ZMAX, RMAX, RMAX, (char*)(".\\output\\app_Bz.csv"));

}

void magnetic_display()
{

}