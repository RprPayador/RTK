#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

const double GM_BDS = 3.986004418e14;
const double Omega_BDS = 7.2921150e-5;

struct Eph {
    double SqrtA, e, i0, iDot, M0, DeltaN, OMEGA, OMEGADot, omega;
    double Cuc, Cus, Crc, Crs, Cic, Cis;
    double toe;
};

void calculate(Eph eph, double t_bds_trans) {
    double tk = t_bds_trans - eph.toe;
    double n0 = sqrt(GM_BDS) / pow(eph.SqrtA, 3);
    double n = n0 + eph.DeltaN;
    double M = eph.M0 + n * tk;
    double E = M;
    for (int i=0; i<10; i++) E = M + eph.e * sin(E);
    
    double v = atan2(sqrt(1 - eph.e*eph.e) * sin(E), cos(E) - eph.e);
    double phi = eph.omega + v;
    double u = phi + eph.Cuc * cos(2*phi) + eph.Cus * sin(2*phi);
    double r = eph.SqrtA * eph.SqrtA * (1 - eph.e * cos(E)) + eph.Crc * cos(2*phi) + eph.Crs * sin(2*phi);
    double i_inc = eph.i0 + eph.iDot * tk + eph.Cic * cos(2*phi) + eph.Cis * sin(2*phi);
    
    double xp = r * cos(u);
    double yp = r * sin(u);
    double L = eph.OMEGA + (eph.OMEGADot - Omega_BDS) * tk - Omega_BDS * eph.toe;
    
    double x = xp * cos(L) - yp * cos(i_inc) * sin(L);
    double y = xp * sin(L) + yp * cos(i_inc) * cos(L);
    double z = yp * sin(i_inc);
    
    cout << fixed << setprecision(4);
    cout << "Time (BDS Trans): " << t_bds_trans << endl;
    cout << "Position: " << x << " " << y << " " << z << endl;
    
    // Rates
    double dE = n / (1 - eph.e * cos(E));
    double dv = dE * sqrt(1 - eph.e*eph.e) / (1 - eph.e * cos(E));
    double du = dv + 2 * dv * (eph.Cus * cos(2*phi) - eph.Cuc * sin(2*phi));
    double dr = eph.e * eph.SqrtA * eph.SqrtA * dE + 2 * dv * (eph.Crs * cos(2*phi) - eph.Crc * sin(2*phi));
    double di = eph.iDot + 2 * dv * (eph.Cis * cos(2*phi) - eph.Cic * sin(2*phi));
    double dL = eph.OMEGADot - Omega_BDS;
    
    double dxp = dr * cos(u) - r * du * sin(u);
    double dyp = dr * sin(u) + r * du * cos(u);
    
    double vx = dxp * cos(L) - xp * dL * sin(L) - dyp * cos(i_inc) * sin(L) + yp * di * sin(i_inc) * sin(L) - yp * cos(i_inc) * dL * cos(L);
    double vy = dxp * sin(L) + xp * dL * cos(L) + dyp * cos(i_inc) * cos(L) - yp * di * sin(i_inc) * cos(L) - yp * cos(i_inc) * dL * sin(L);
    double vz = dyp * sin(i_inc) + yp * di * cos(i_inc);
    
    cout << "Velocity: " << vx << " " << vy << " " << vz << endl;
    cout << "Norm: " << sqrt(vx*vx + vy*vy + vz*vz) << endl;
}

int main() {
    Eph eph;
    eph.SqrtA = 6.493744514465332E+03;
    eph.e = 4.696009680628777E-03;
    eph.i0 = 8.322711880299253E-01;
    eph.iDot = -2.507247294056509E-10;
    eph.M0 = -2.611632933657392E+00;
    eph.DeltaN = 1.587923286235789E-09;
    eph.OMEGA = 1.873874061141915E+00;
    eph.OMEGADot = -2.111516524567248E-09;
    eph.omega = -2.581484528323869E+00;
    eph.Cuc = 1.340871676802635E-05;
    eph.Cus = 3.649387508630753E-06;
    eph.Crc = 9.704687500000000E+01;
    eph.Crs = 4.088593750000000E+02;
    eph.Cic = -1.448206603527069E-07;
    eph.Cis = 2.468004822731018E-08;
    eph.toe = 273600.000;
    
    // Search around 273995 - 14.12
    for (double dt = -14.2; dt < -14.0; dt += 0.001) {
        double t = 273995.0 + dt;
        double tk = t - eph.toe;
        double n = sqrt(GM_BDS) / pow(eph.SqrtA, 3) + eph.DeltaN;
        double M = eph.M0 + n * tk;
        double E = M; for (int k=0; k<10; k++) E = M + eph.e * sin(E);
        double v = atan2(sqrt(1 - eph.e*eph.e)*sin(E), cos(E)-eph.e);
        double phi = eph.omega + v;
        double u = phi + eph.Cuc*cos(2*phi) + eph.Cus*sin(2*phi);
        double r = eph.SqrtA*eph.SqrtA*(1-eph.e*cos(E)) + eph.Crc*cos(2*phi) + eph.Crs*sin(2*phi);
        double i = eph.i0 + eph.iDot*tk + eph.Cic*cos(2*phi) + eph.Cis*sin(2*phi);
        double L = eph.OMEGA + (eph.OMEGADot - Omega_BDS) * tk - Omega_BDS * eph.toe;
        double xp = r * cos(u); double yp = r * sin(u);
        double x = xp * cos(L) - yp * cos(i) * sin(L);
        if (fabs(x - (-3556091.913)) < 1.0) {
            calculate(eph, t);
            break;
        }
    }
    return 0;
}
