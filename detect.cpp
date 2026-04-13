#include "RTK_Structs.h"
#include <cmath>

double Hopfield(double H,double Elev){
    double RH=RH0*(exp(-0.0006396*(H-HIGHT0)));
    double p=PRESSURE0*pow((1-0.0000226*(H-HIGHT0)),5.225);
    double T=TEMPERATURE0-0.0065*(H-HIGHT0);
    double e=RH*exp(-37.2465+0.213166*T-0.000256908*T*T);
    double h_w=11000;
    double h_d=40136+148.72*(TEMPERATURE0-273.16);
    double K_w=155.2e-7*4810/T/T*e*(h_w-H);
    double K_d=155.2e-7*p/T*(h_d-H);
    return K_d/sin(sqrt(Elev*Elev+6.25))+K_w/sin(sqrt(Elev*Elev+2.25));
}


void DetectOutlier(EPOCHOBS* Obs){
    int validcount=0;
    MWGF CurComObs[MAXCHANNUM];
    for(int i=0;i<Obs->SatNum;i++){
        SATOBS& sat=Obs->SatObs[i];

        CurComObs[i].Prn=sat.Prn;
        CurComObs[i].Sys=sat.System;

        if(sat.c1==0.0||sat.l1==0.0||sat.p2==0.0||sat.l2==0.0) {
            sat.Valid=false;
            continue;
        }
        double GF=sat.l1-sat.l2;
        double MW,PIF;
        if(sat.System==GPS){
            MW=(FG1_GPS*sat.l1-FG2_GPS*sat.l2)/(FG1_GPS-FG2_GPS)-(FG1_GPS*sat.c1+FG2_GPS*sat.p2)/(FG1_GPS+FG2_GPS);
            PIF=(FG1_GPS*FG1_GPS*sat.c1-FG2_GPS-FG2_GPS*sat.p2)/(FG1_GPS*FG1_GPS-FG2_GPS*FG2_GPS);
        }

        else{
            MW=(FG1_BDS*sat.l1-FG2_BDS*sat.l2)/(FG1_BDS-FG2_BDS)-(FG1_BDS*sat.c1+FG2_BDS*sat.p2)/(FG1_BDS+FG2_BDS);
            PIF=(FG1_BDS*FG1_BDS*sat.c1-FG2_BDS-FG2_BDS*sat.p2)/(FG1_BDS*FG1_BDS-FG2_BDS*FG2_BDS);
        } 
        double dGF=fabs(GF-Obs->ComObs[i].GF);
        double dMW=fabs(MW-Obs->ComObs[i].MW);
        if(dGF>0.05||dMW>3){
            sat.Valid=false;
        }
        else{
            sat.Valid=true;
            validcount++;
            CurComObs[i].GF=GF;
            CurComObs[i].MW=(Obs->ComObs[i].n*Obs->ComObs[i].MW+MW)/(Obs->ComObs[i].n+1);
            CurComObs[i].n=Obs->ComObs[i].n+1;
            CurComObs[i].PIF=PIF;  
        }
    }
    memcpy(CurComObs,Obs->ComObs,MAXCHANNUM*sizeof(MWGF));
}