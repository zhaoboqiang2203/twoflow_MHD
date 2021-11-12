#include "MHD.h"

double c0 = 2e-7;// mu0 / 2pi = 2e-7
double I = 18000;// current, unit: A
double a = 0.5;// radius of current coil, unit: m
double eps = 1e-4;

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
    return k2_0(x, y, z) - eps * (k2_0(x, y, z) == 1);
}

double K1(double x, double y, double z)
{
    return PI / 2 * (1 + sqr(1 / 2) * k2(x, y, z) + sqr(1 / 2 * 3 / 4) * sqr(k2(x, y, z)) + sqr(1 / 2 * 3 / 4 * 5 / 6) * cube(k2(x, y, z)));
}
double E(double x, double y, double z)
{
    return PI / 2. * (1 - sqr(1 / 2) * k2(x, y, z) + sqr(1 / 2 * 3 / 4) * pow(k2(x, y, z), 2 / 3) - sqr(1 / 2 * 3 / 4 * 5 / 6) * pow(k2(x, y, z), 3 / 5));
}
double a1(double x, double y, double z)
{
    return sqrt(sqr(a + rho1(x, y, z)) + sqr(z));
}
double b1_0(double x, double y, double z)
{
    return sqr(a + rho1(x, y, z)) + sqr(z);
}
double b1(double x, double y, double z)
{
    return b1_0(x, y, z) + eps * (b1_0(x, y, z) == 0);
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
    return c0 * I * z  / ((rho1(x, y, z) + eps * (rho1(x, y, z) == 0)) * a1(x, y, z)) * (b2(x, y, z)  / b1(x, y, z) * E(x, y, z) - K1(x, y, z));
}

double Bscz(double x, double y, double z)
{
    return c0 * I / a1(x, y, z) * (b3(x, y, z) / b1(x, y, z) * E(x, y, z) + K1(x, y, z));
}

double Btheta(double x, double y, double z)
{
    return 0;
}

void magnetic_field_initial()
{
    double bias;
    a = 0.1;

    bias = 0;
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            app_Br[i][j] = Brho(i * dr - bias, 0, j * dz);
            app_Bz[i][j] = Bscz(i * dr - bias, 0, j * dz);
        }
    }

    matrix_to_csv((double**)app_Br, ZMAX, RMAX, RMAX, (char*)(".\\output\\app_Br.csv"));
    matrix_to_csv((double**)app_Bz, ZMAX, RMAX, RMAX, (char*)(".\\output\\app_Bz.csv"));

}