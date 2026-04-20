#include "RTK_Structs.h"
#include <cmath>

double Hopfield(double H,double Elev){
    // 增加稳健性判断：若高度异常（如迭代初期在地心），返回0以防止 exp 溢出
    if (H < -1000.0 || H > 50000.0) return 0.0;

    double RH=RH0*(exp(-0.0006396*(H-HIGHT0)));
    double p=PRESSURE0*pow((1-0.0000226*(H-HIGHT0)),5.225);
    double T=TEMPERATURE0-0.0065*(H-HIGHT0);
    double e=RH*exp(-37.2465+0.213166*T-0.000256908*T*T);
    double h_w=11000;
    double h_d=40136+148.72*(TEMPERATURE0-273.16);
    double K_w=155.2e-7*4810/T/T*e*(h_w-H);
    double K_d=155.2e-7*p/T*(h_d-H);
    double degEl = Elev * Deg;
    return K_d/sin(sqrt(degEl*degEl+6.25)*Rad)+K_w/sin(sqrt(degEl*degEl+2.25)*Rad);
}


void DetectOutlier(EPOCHOBS* Obs){
    int validcount=0;
    MWGF CurComObs[MAXCHANNUM];
    
    // 1. 备份并清除状态，准备记录当前历元有效信息
    MWGF PreComObs[MAXCHANNUM];
    memcpy(PreComObs, Obs->ComObs, MAXCHANNUM * sizeof(MWGF));
    memset(Obs->ComObs, 0, MAXCHANNUM * sizeof(MWGF));

    for(int i=0; i<Obs->SatNum; i++){
        SATOBS& sat = Obs->SatObs[i];
        
        // 初始化当前卫星的记录
        CurComObs[i].Prn = sat.Prn;
        CurComObs[i].Sys = sat.System;

        if(sat.c1==0.0 || sat.l1==0.0 || sat.p2==0.0 || sat.l2==0.0) {
            sat.Valid = false;
            continue;
        }

        // 2. 计算当前组合值
        double GF = sat.l1 - sat.l2;
        double MW, PIF;
        if(sat.System == GPS){
            MW = (FG1_GPS*sat.l1 - FG2_GPS*sat.l2)/(FG1_GPS - FG2_GPS) - (FG1_GPS*sat.c1 + FG2_GPS*sat.p2)/(FG1_GPS + FG2_GPS);
            PIF = (FG1_GPS*FG1_GPS*sat.c1 - FG2_GPS*FG2_GPS*sat.p2)/(FG1_GPS*FG1_GPS - FG2_GPS*FG2_GPS);
        }
        else {
            // BDS B1/B3 组合 (注意：之前代码误用了 FG2_BDS)
            MW = (FG1_BDS*sat.l1 - FG3_BDS*sat.l2)/(FG1_BDS - FG3_BDS) - (FG1_BDS*sat.c1 + FG3_BDS*sat.p2)/(FG1_BDS + FG3_BDS);
            PIF = (FG1_BDS*FG1_BDS*sat.c1 - FG3_BDS*FG3_BDS*sat.p2)/(FG1_BDS*FG1_BDS - FG3_BDS*FG3_BDS);
        }

        // 3. 搜索上一历元匹配项
        int preIdx = -1;
        for(int j=0; j<MAXCHANNUM; j++){
            if(PreComObs[j].Prn == sat.Prn && PreComObs[j].Sys == sat.System){
                preIdx = j;
                break;
            }
        }

        // 4. 周跳检测与状态平滑
        if(preIdx != -1 && PreComObs[preIdx].n > 0){
            double dGF = fabs(GF - PreComObs[preIdx].GF);
            double dMW = fabs(MW - PreComObs[preIdx].MW);
            
            if(dGF > 0.05 || dMW > 3.0){
                sat.Valid = false;
                CurComObs[i].MW = MW;
                CurComObs[i].n = 1;
            }
            else{
                sat.Valid = true;
                validcount++;
                CurComObs[i].MW = (PreComObs[preIdx].n * PreComObs[preIdx].MW + MW) / (PreComObs[preIdx].n + 1);
                CurComObs[i].n = PreComObs[preIdx].n + 1;
            }
        }
        else {
            // 新卫星接入
            sat.Valid = true;
            CurComObs[i].MW = MW;
            CurComObs[i].n = 1;
            validcount++;
        }
        
        CurComObs[i].GF = GF;
        CurComObs[i].PIF = PIF;
    }

    // 5. 将当前计算结果更新回持久化缓冲区
    memcpy(Obs->ComObs, CurComObs, MAXCHANNUM * sizeof(MWGF));
}