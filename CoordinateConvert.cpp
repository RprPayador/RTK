#include "RTK_Structs.h"
#include <cmath>

// 空间直角坐标转大地坐标
void XYZToBLH(const double XYZ[3], double BLH[3], const double R, const double F)
{
    double X = XYZ[0];
    double Y = XYZ[1];
    double Z = XYZ[2];
    
    double e2 = 2.0 * F - F * F;
    double p = sqrt(X * X + Y * Y);
    
    double L = atan2(Y, X);
    double B = atan2(Z, p * (1.0 - e2));
    double B0 = 0.0;
    double N = 0.0;
    
    // 迭代求解纬度 B
    while (fabs(B - B0) > 1e-10) {
        B0 = B;
        N = R / sqrt(1.0 - e2 * sin(B0) * sin(B0));
        B = atan2(Z + N * e2 * sin(B0), p);
    }
    
    double H = p / cos(B) - N;
    
    BLH[0] = B;
    BLH[1] = L;
    BLH[2] = H;
}

// 大地坐标转空间直角坐标
void BLHToXYZ(const double BLH[3], double XYZ[3], const double R, const double F)
{
    double B = BLH[0];
    double L = BLH[1];
    double H = BLH[2];
    
    double e2 = 2.0 * F - F * F;
    double N = R / sqrt(1.0 - e2 * sin(B) * sin(B));
    
    XYZ[0] = (N + H) * cos(B) * cos(L);
    XYZ[1] = (N + H) * cos(B) * sin(L);
    XYZ[2] = (N * (1.0 - e2) + H) * sin(B);
}

// 计算从地心地固坐标系(ECEF)到站心坐标系(ENU)的旋转矩阵
void BLHToNEUMat(const double BLH[], double Mat[])
{
    double B = BLH[0];
    double L = BLH[1];
    
    double sinB = sin(B);
    double cosB = cos(B);
    double sinL = sin(L);
    double cosL = cos(L);
    
    Mat[0] = -sinL;
    Mat[1] = cosL;
    Mat[2] = 0.0;
    
    Mat[3] = -sinB * cosL;
    Mat[4] = -sinB * sinL;
    Mat[5] = cosB;
    
    Mat[6] = cosB * cosL;
    Mat[7] = cosB * sinL;
    Mat[8] = sinB;
}

// 卫星高度角方位角计算函数
void CompSatElAz(const double Xr[], const double Xs[], double *Elev, double *Azim)
{
    double dx = Xs[0] - Xr[0];
    double dy = Xs[1] - Xr[1];
    double dz = Xs[2] - Xr[2];
    
    double blh[3];

    XYZToBLH(Xr, blh, R_WGS84, F_WGS84);
    
    double mat[9];
    BLHToNEUMat(blh, mat);
    

    double de = mat[0] * dx + mat[1] * dy + mat[2] * dz;
    double dn = mat[3] * dx + mat[4] * dy + mat[5] * dz;
    double dh = mat[6] * dx + mat[7] * dy + mat[8] * dz;
    
    if (Elev) {
        *Elev = atan2(dh, sqrt(dn * dn + de * de));
    }
    if (Azim) {
        double az = atan2(de, dn); 
        if (az < 0.0) az += 2.0 * PAI;
        *Azim = az;
    }
}

// 定位误差计算函数
void Comp_dEnu(const double X0[], const double Xr[], double dNeu[])
{
    double dx = Xr[0] - X0[0];
    double dy = Xr[1] - X0[1];
    double dz = Xr[2] - X0[2];
    
    double blh[3];

    XYZToBLH(X0, blh, R_WGS84, F_WGS84);
    
    double mat[9];
    BLHToNEUMat(blh, mat);
    
    // 投影到站心坐标系计算东北天误差 (N, E, U)
    dNeu[0] = mat[0] * dx + mat[1] * dy + mat[2] * dz; // N
    dNeu[1] = mat[3] * dx + mat[4] * dy + mat[5] * dz; // E
    dNeu[2] = mat[6] * dx + mat[7] * dy + mat[8] * dz; // U
}
