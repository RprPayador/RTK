#include "RTK_Structs.h"
#include <math.h>

double GetDiffTime(const GPSTIME* t2, const GPSTIME* t1) {
    return (double)(t2->Week - t1->Week) * 604800.0 + (t2->SecOfWeek - t1->SecOfWeek);
}

bool CompSatClkOff(const int Prn, const GNSSSys Sys, const GPSTIME* t, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, SATMIDRES* Mid) {
    GPSEPHREC* eph = NULL;

    if (Sys == GPS) {
        eph = &GPSEph[Prn - 1];
    }
    else if (Sys == BDS) {
        eph = &BDSEph[Prn - 1];
    }

    if (eph == NULL || eph->PRN == 0) return false;

    GPSTIME t_sys = *t;
    if (Sys == BDS) { //如果是BDS，从GPST转换到BDST，再计算到TOC的间隔时间
        t_sys.SecOfWeek -= 14.0;
        t_sys.Week -= 1356;
        if (t_sys.SecOfWeek < 0) { t_sys.SecOfWeek += 604800.0; t_sys.Week--; }
    }

    double dt = GetDiffTime(&t_sys, &eph->TOC);//确保相同时间系统间求差

    Mid->SatClkOft = eph->ClkBias + eph->ClkDrift * dt + eph->ClkDriftRate * dt * dt;//钟差
    Mid->SatClkSft = eph->ClkDrift + 2 * eph->ClkDriftRate * dt;//钟速

    Mid->Tgd1 = eph->TGD1;
    if (Sys == GPS) {
        Mid->Tgd2 = (77.0 / 60.0) * (77.0 / 60.0) * eph->TGD1;
    } else {
        Mid->Tgd2 = eph->TGD2;
    }

    return true;
}


bool CompGPSSatPVT(const int Prn, const GPSTIME* t, const GPSEPHREC* eph, SATMIDRES* Mid) {
    if (eph == NULL || eph->PRN == 0) return false;

    double tk = GetDiffTime(t, &eph->TOE);

    double n0 = sqrt(GM_Earth) / pow(eph->SqrtA, 3);//computed mean motion
    double n = n0 + eph->DeltaN;//corrected mean motion
    double M = eph->M0 + n * tk;//mean anomaly
    double E = M, E_old;//eccentric anomaly by iteration
    for (int i = 0; i < 20; i++) {
        E_old = E;
        E = M + eph->e * sin(E);
        if (fabs(E - E_old) < 1.0e-12) break;
    }
    
    double F = -2.0 * sqrt(GM_Earth) / (C_Light * C_Light);
    Mid->SatClkOft += F * eph->e * eph->SqrtA * sin(E);//钟差的相对论改正

    double v = atan2(sqrt(1 - eph->e * eph->e) * sin(E), cos(E) - eph->e);//true anomaly
    double phi = eph->omega + v;//argument of latitude
    double du = eph->Cuc * cos(2 * phi) + eph->Cus * sin(2 * phi);//second harmonic
    double dr = eph->Crc * cos(2 * phi) + eph->Crs * sin(2 * phi);
    double di = eph->Cic * cos(2 * phi) + eph->Cis * sin(2 * phi);
    double u = phi + du;//corrected latitude
    double r = eph->SqrtA * eph->SqrtA * (1 - eph->e * cos(E)) + dr;//corrected radius
    double i = eph->i0 + eph->iDot * tk + di;//corrected inclination

    double x_plane = r * cos(u);//x position in orbital
    double y_plane = r * sin(u);//y position in orbital

    double L = eph->OMEGA + (eph->OMEGADot - Omega_WGS) * tk - Omega_WGS * eph->TOE.SecOfWeek;//Corrected longitude of ascending node(Omeaga_k)

    Mid->SatPos[0] = x_plane * cos(L) - y_plane * cos(i) * sin(L);//earth fixed
    Mid->SatPos[1] = x_plane * sin(L) + y_plane * cos(i) * cos(L);
    Mid->SatPos[2] = y_plane * sin(i);

    double dE=n/(1-eph->e*cos(E));//eccentric anomaly rate
    
    Mid->SatClkSft += F * eph->e * eph->SqrtA * cos(E) * dE;//钟速的相对论改正

    double dv=dE*sqrt(1-eph->e*eph->e)/(1-eph->e*cos(E));//true anomaly rate
    di=eph->iDot+2*dv*(eph->Cis*cos(2*phi)-eph->Cic*sin(2*phi));//Corrected Inclination Angle Rate
    du=dv+2*dv*(eph->Cus*cos(2*phi)-eph->Cuc*sin(2*phi));//Corrected Argument of Latitude Rate
    dr=eph->e*eph->SqrtA*eph->SqrtA*dE*sin(E)+2*dv*(eph->Crs*cos(2*phi)-eph->Crc*sin(2*phi));//Corrected Radius Rate
    double dL=eph->OMEGADot-Omega_WGS;//Longitude of Ascending Node Rate
    
    double dx_plane=dr*cos(u)-r*du*sin(u);//In-plane xvelocity
    double dy_plane=dr*sin(u)+r*du*cos(u);//In-plane xvelocity
  
    Mid->SatVel[0]=-x_plane*dL*sin(L)+dx_plane*cos(L)-dy_plane*sin(L)*cos(i)-y_plane*(dL*cos(L)*cos(i)-di*sin(L)*sin(i));//Earth-Fixed x velocity 
    Mid->SatVel[1]=x_plane*dL*cos(L)+dx_plane*sin(L)+dy_plane*cos(L)*cos(i)-y_plane*(dL*sin(L)*cos(i)+di*sin(L)*sin(i));//Earth-Fixed y velocity 
    Mid->SatVel[2]=dy_plane*sin(i)+y_plane*di*cos(i);//Earth-Fixed z velocity 

    Mid->Valid = true;
    return true;
}


bool CompBDSSatPVT(const int Prn, const GPSTIME* t, const GPSEPHREC* eph, SATMIDRES* Mid) {
    if (eph == NULL || eph->PRN == 0) return false;

    /*GPST转换到BDST*/
    GPSTIME t_bds = *t;
    t_bds.SecOfWeek -= 14.0;
    t_bds.Week -= 1356;
    if (t_bds.SecOfWeek < 0) {
        t_bds.SecOfWeek += 604800.0;
        t_bds.Week -= 1;
    }
    // 确保 tk 是相对于 TOE 的合理数值（应为正负几小时内）
    double tk = GetDiffTime(&t_bds, &eph->TOE);
    /*
    if (tk > 302400.0) tk -= 604800.0;
    else if (tk < -302400.0) tk += 604800.0;
    */


    double n0 = sqrt(GM_BDS) / pow(eph->SqrtA, 3);//卫星平均角速度
    double n = n0 + eph->DeltaN;//改正平均角速度
    double M = eph->M0 + n * tk;//平近点角
    double E = M, E_old;//偏近点角
    for (int k = 0; k < 20; k++) {
        E_old = E;
        E = M + eph->e * sin(E);
        if (fabs(E - E_old) < 1.0e-12) break;
    }

    double F = -2.0 * sqrt(GM_BDS) / (C_Light * C_Light);
    Mid->SatClkOft += F * eph->e * eph->SqrtA * sin(E);//钟差的相对论改正

    double v = atan2(sqrt(1 - eph->e * eph->e) * sin(E), cos(E) - eph->e);//真近点角
    double phi = eph->omega + v;//计算纬度幅角
    double du = eph->Cuc * cos(2 * phi) + eph->Cus * sin(2 * phi);//纬度幅角改正项
    double dr = eph->Crc * cos(2 * phi) + eph->Crs * sin(2 * phi);//径向改正项
    double di = eph->Cic * cos(2 * phi) + eph->Cis * sin(2 * phi);//轨道倾角改正项
    double u = phi + du;//改正后的纬度幅角
    double r = eph->SqrtA * eph->SqrtA * (1 - eph->e * cos(E)) + dr;//改正后的径向
    double i = eph->i0 + eph->iDot * tk + di;//改正后的轨道倾角

    double x_plane = r * cos(u);//轨道平面内的x坐标
    double y_plane = r * sin(u);//轨道平面内的y坐标

    bool is_geo = (Prn >= 1 && Prn <= 5) || (Prn >= 59 && Prn <= 63);//判断是否为GEO卫星

    if (is_geo) {
        double L = eph->OMEGA + eph->OMEGADot * tk - Omega_BDS * eph->TOE.SecOfWeek;//历元升交点经度（惯性系）(Omega_k)
        double x_gk = x_plane * cos(L) - y_plane * cos(i) * sin(L);//自定义坐标系中的x坐标
        double y_gk = x_plane * sin(L) + y_plane * cos(i) * cos(L);//自定义坐标系中的y坐标
        double z_gk = y_plane * sin(i);//自定义坐标系中的z坐标

        double ang_phi = Omega_BDS * tk;
        double Rx = -5.0 * PAI / 180.0;

        double xt = x_gk;
        double yt = y_gk * cos(Rx) + z_gk * sin(Rx);
        double zt = -y_gk * sin(Rx) + z_gk * cos(Rx);

        Mid->SatPos[0] = xt * cos(ang_phi) + yt * sin(ang_phi);//BDCS坐标系中的坐标
        Mid->SatPos[1] = -xt * sin(ang_phi) + yt * cos(ang_phi);
        Mid->SatPos[2] = zt;

        double dE=n/(1-eph->e*cos(E));//eccentric anomaly rate
        Mid->SatClkSft += F * eph->e * eph->SqrtA * cos(E) * dE;//钟速的相对论改正
        double dv=dE*sqrt(1-eph->e*eph->e)/(1-eph->e*cos(E));//true anomaly rate
        di=eph->iDot+2*dv*(eph->Cis*cos(2*phi)-eph->Cic*sin(2*phi));//Corrected Inclination Angle Rate
        du=dv+2*dv*(eph->Cus*cos(2*phi)-eph->Cuc*sin(2*phi));//Corrected Argument of Latitude Rate
        dr=eph->e*eph->SqrtA*eph->SqrtA*dE*sin(E)+2*dv*(eph->Crs*cos(2*phi)-eph->Crc*sin(2*phi));//Corrected Radius Rate
        double dL=eph->OMEGADot;//Longitude of Ascending Node Rate

        double dx_plane=dr*cos(u)-r*du*sin(u);//In-plane xvelocity
        double dy_plane=dr*sin(u)+r*du*cos(u);//In-plane xvelocity

        double dx_gk=-x_plane*dL*sin(L)+dx_plane*cos(L)-dy_plane*sin(L)*cos(i)-y_plane*(dL*cos(L)*cos(i)-di*sin(L)*sin(i));//自定义坐标系下x速度
        double dy_gk=x_plane*dL*cos(L)+dx_plane*sin(L)+dy_plane*cos(L)*cos(i)-y_plane*(dL*sin(L)*cos(i)+di*cos(L)*sin(i));//自定义坐标系下y速度
        double dz_gk=dy_plane*sin(i)+y_plane*di*cos(i);//自定义坐标系下z速度

        double c=cos(-5.0/180.0*PAI);
        double s=sin(-5.0/180.0*PAI);
        double theta=Omega_BDS*tk;
        Mid->SatVel[0] = cos(theta)*dx_gk + c*sin(theta)*dy_gk + s*sin(theta)*dz_gk+ Omega_BDS*(-sin(theta)*x_gk + c*cos(theta)*y_gk + s*cos(theta)*z_gk);
        Mid->SatVel[1] = -sin(theta)*dx_gk + c*cos(theta)*dy_gk + s*cos(theta)*dz_gk+ Omega_BDS*(-cos(theta)*x_gk - c*sin(theta)*y_gk - s*sin(theta)*z_gk);
        Mid->SatVel[2] = -s*dy_gk + c*dz_gk; 




    } else {
        double L = eph->OMEGA + (eph->OMEGADot - Omega_BDS) * tk - Omega_BDS * eph->TOE.SecOfWeek;//历元升交点经度（地固系）
        Mid->SatPos[0] = x_plane * cos(L) - y_plane * cos(i) * sin(L);//BDCS坐标系中的x坐标
        Mid->SatPos[1] = x_plane * sin(L) + y_plane * cos(i) * cos(L);//BDCS坐标系中的y坐标
        Mid->SatPos[2] = y_plane * sin(i);//BDCS坐标系中的z坐标

        /*速度计算同GPS*/
        double dE=n/(1-eph->e*cos(E));//eccentric anomaly rate
        Mid->SatClkSft += F * eph->e * eph->SqrtA * cos(E) * dE;//钟速的相对论改正
        double dv=dE*sqrt(1-eph->e*eph->e)/(1-eph->e*cos(E));//true anomaly rate
        di=eph->iDot+2*dv*(eph->Cis*cos(2*phi)-eph->Cic*sin(2*phi));//Corrected Inclination Angle Rate
        du=dv+2*dv*(eph->Cus*cos(2*phi)-eph->Cuc*sin(2*phi));//Corrected Argument of Latitude Rate
        dr=eph->e*eph->SqrtA*eph->SqrtA*dE*sin(E)+2*dv*(eph->Crs*cos(2*phi)-eph->Crc*sin(2*phi));//Corrected Radius Rate
        double dL=eph->OMEGADot-Omega_BDS;//Longitude of Ascending Node Rate

        double dx_plane=dr*cos(u)-r*du*sin(u);//In-plane xvelocity
        double dy_plane=dr*sin(u)+r*du*cos(u);//In-plane xvelocity

        Mid->SatVel[0]=-x_plane*dL*sin(L)+dx_plane*cos(L)-dy_plane*sin(L)*cos(i)-y_plane*(dL*cos(L)*cos(i)-di*sin(L)*sin(i));//Earth-Fixed x velocity 
        Mid->SatVel[1]=x_plane*dL*cos(L)+dx_plane*sin(L)+dy_plane*cos(L)*cos(i)-y_plane*(dL*sin(L)*cos(i)+di*cos(L)*sin(i));//Earth-Fixed y velocity 
        Mid->SatVel[2]=dy_plane*sin(i)+y_plane*di*cos(i);//Earth-Fixed z velocity 

    }

    Mid->Valid = true;
    return true;
}
