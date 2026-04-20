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

void get_pos(Eph eph, double t, double* pos) {
    double tk = t - eph.toe;
    if (tk > 302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;

    double n0 = sqrt(GM_BDS) / pow(eph.SqrtA, 3);
    double n = n0 + eph.DeltaN;
    double M = eph.M0 + n * tk;
    double E = M;
    for (int k=0; k<10; k++) E = M + eph.e * sin(E);
    
    double v = atan2(sqrt(1 - eph.e*eph.e) * sin(E), cos(E) - eph.e);
    double phi = eph.omega + v;
    double u = phi + eph.Cuc * cos(2*phi) + eph.Cus * sin(2*phi);
    double r = eph.SqrtA * eph.SqrtA * (1 - eph.e * cos(E)) + eph.Crc * cos(2*phi) + eph.Crs * sin(2*phi);
    double i = eph.i0 + eph.iDot * tk + eph.Cic * cos(2*phi) + eph.Cis * sin(2*phi);
    
    double xp = r * cos(u);
    double yp = r * sin(u);
    double L = eph.OMEGA + (eph.OMEGADot - Omega_BDS) * tk - Omega_BDS * eph.toe;
    
    pos[0] = xp * cos(L) - yp * cos(i) * sin(L);
    pos[1] = xp * sin(L) + yp * cos(i) * cos(L);
    pos[2] = yp * sin(i);
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

    double target_x = -3556091.913;
    
    cout << fixed << setprecision(7);
    for (double dt = -14.15; dt < -14.10; dt += 0.001) {
        double t = 273995.0 + dt;
        double p[3];
        get_pos(eph, t, p);
        if (fabs(p[0] - target_x) < 5.0) {
            double p1[3], p2[3];
            double delta = 0.01;
            get_pos(eph, t - delta, p1);
            get_pos(eph, t + delta, p2);
            cout << "dt: " << dt << " Pos: " << p[0] << " " << p[1] << " " << p[2] << endl;
            cout << "Vel: " << (p2[0]-p1[0])/(2*delta) << " " << (p2[1]-p1[1])/(2*delta) << " " << (p2[2]-p1[2])/(2*delta) << endl;
        }
    }
    return 0;
}
