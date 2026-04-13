#include <iostream>
#include <iomanip>
#include "RTK_Structs.h"

int main()
{
    FILE* fp, * fout;
    unsigned char Buff[MAXRAWLEN];
    int Len, LenRead, val;
    EPOCHOBS Obs;
    GPSEPHREC GpsEph[MAXGPSNUM], BdsEph[MAXBDSNUM];

    if ((fp = fopen("D:\\GNSS Algorithm\\RTK\\RTK\\oem719-202603111200.bin", "rb")) == NULL) { // C4996
        printf("The file 'oem719-202603111200.bin' was not opened\n");
        return 0;
    }

    Len = 0;
    do {
        LenRead = fread(Buff + Len, sizeof(unsigned char), MAXRAWLEN - Len, fp);
        Len = Len + LenRead;

        val = DecodeNovOem7Dat(Buff, Len, &Obs, GpsEph, BdsEph);

        if (val == 1) { // 成功解码出一个历元的观测值
            for (int i = 0; i < Obs.SatNum; i++) {
                int prn = Obs.SatObs[i].Prn;
                GNSSSys sys = Obs.SatObs[i].System;

                // 计算卫星钟差和位置，计算出的结果保存在 Obs.SatPVT[i] 中
                CompSatClkOff(prn, sys, &Obs.Time, GpsEph, BdsEph, &Obs.SatPVT[i]);
                if (sys == GPS) {
                    CompGPSSatPVT(prn, &Obs.Time, &GpsEph[prn - 1], &Obs.SatPVT[i]);
                }
                else if (sys == BDS) {
                    CompBDSSatPVT(prn, &Obs.Time, &BdsEph[prn - 1], &Obs.SatPVT[i]);
                }
            }
        }

    } while (!feof(fp));
    return 1;
}

