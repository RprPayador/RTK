#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include "RTK_Structs.h"
#include "sockets.cpp"

int main()
{
    FILE* fp;
    std::ofstream fout("D:\\GNSS Algorithm\\RTK_Raw\\result.txt");
    unsigned char Buff[MAXRAWLEN];
    int Len = 0, LenRead, val;
    SOCKET NetGps;

    if(FILEMODE == 0){
        if(OpenSocket(NetGps, "47.114.134.129", 7190)==false){
            printf( "This ip& port was not opened.\n");
            return 0;
        }
    }
    
    EPOCHOBS Obs; 
    memset(&Obs, 0, sizeof(Obs));
    GPSEPHREC GpsEph[MAXGPSNUM], BdsEph[MAXBDSNUM];

    memset(GpsEph, 0, sizeof(GpsEph));
    memset(BdsEph, 0, sizeof(BdsEph));

    if ((fp = fopen("D:\\GNSS Algorithm\\RTK\\RTK\\oem719-202603111200.bin", "rb")) == NULL) { 
        printf("The file 'oem719-202603111200.bin' was not opened\n");
        return 0;
    }

    if (!fout.is_open()) {
        printf("Failed to open result.txt for writing\n");
    }

    printf("Starting SPP Debugging Flow...\n");

    int total_epochs = 0;
    int success_epochs = 0;

    do {
        if(FILEMODE == 1){
            LenRead = fread(Buff + Len, sizeof(unsigned char), MAXRAWLEN - Len, fp);
            if (LenRead <= 0 && Len <= 0) break;
            Len = Len + LenRead;
        }
        if(FILEMODE == 0){
            LenRead=recv(NetGps, (char*)Buff, MAXRAWLEN, 0);
            if (LenRead <= 0 && Len <= 0) break;
            Len = Len + LenRead;
        }


        int initial_len = Len;
        val = DecodeNovOem7Dat(Buff, Len, &Obs, GpsEph, BdsEph);

        if (val == 2) {
            // 如果 SPP 还没成功过，用 PSRPOS 作为初值
            PPRESULT psr_res;
            decode_psrpos(Buff, &psr_res);
            if (Obs.Pos[0] == 0 && Obs.Pos[1] == 0 && Obs.Pos[2] == 0) {
                Obs.Pos[0] = psr_res.Position[0];
                Obs.Pos[1] = psr_res.Position[1];
                Obs.Pos[2] = psr_res.Position[2];
            }
        }

        if (val == 1) { 
            total_epochs++;
            // 粗差探测
            DetectOutlier(&Obs);

            // 统计通过探测的卫星
            int detected_sats = 0;
            for(int i=0; i<Obs.SatNum; i++) if(Obs.SatObs[i].Valid) detected_sats++;

            // 单点定位解算
            PPRESULT Res;
            RAWDAT Raw;
            for (int i = 0; i < MAXGPSNUM; i++) Raw.GpsEph[i] = GpsEph[i];
            for (int i = 0; i < MAXBDSNUM; i++) Raw.BdsEph[i] = BdsEph[i];

            if (SPP(&Obs, &Raw, &Res)) {
                // 执行单点测速
                SPV(&Obs, &Res);

                success_epochs++;
                
                // 详细卫星信息输出
                fout << std::fixed << std::uppercase;
                for (int i = 0; i < Obs.SatNum; i++) {
                    if (!Obs.SatPVT[i].Valid) continue;
                    SATOBS& sat = Obs.SatObs[i];
                    SATMIDRES& pvt = Obs.SatPVT[i];
                    
                    fout << (sat.System == GPS ? 'G' : 'C') << std::setfill('0') << std::setw(2) << (int)sat.Prn << std::setfill(' ')
                         << " X=" << std::setw(14) << std::setprecision(3) << pvt.SatPos[0]
                         << " Y=" << std::setw(14) << std::setprecision(3) << pvt.SatPos[1]
                         << " Z=" << std::setw(14) << std::setprecision(3) << pvt.SatPos[2]
                         << " Clk=" << std::scientific << std::setw(15) << std::setprecision(6) << pvt.SatClkOft
                         << std::fixed << " Vx=" << std::setw(11) << std::setprecision(4) << pvt.SatVel[0]
                         << " Vy=" << std::setw(11) << std::setprecision(4) << pvt.SatVel[1]
                         << " Vz=" << std::setw(11) << std::setprecision(4) << pvt.SatVel[2]
                         << " Clkd=" << std::scientific << std::setw(15) << std::setprecision(5) << pvt.SatClkSft
                         << std::fixed << " PIF=" << std::setw(14) << std::setprecision(4) << Obs.ComObs[i].PIF
                         << " Trop=" << std::setw(8) << std::setprecision(3) << pvt.TropCorr
                         << " E=" << std::setw(7) << std::setprecision(3) << pvt.Elevation * Deg << "deg" << std::endl;
                }

                // SPP 汇总输出
                double blh[3];
                XYZToBLH(Res.Position, blh, R_WGS84, F_WGS84);
                
                fout << "SPP: " << (int)Res.Time.Week << " " << std::fixed << std::setprecision(3) << Res.Time.SecOfWeek
                     << " X:" << std::setprecision(4) << Res.Position[0]
                     << " Y:" << std::setprecision(4) << Res.Position[1]
                     << " Z:" << std::setprecision(4) << Res.Position[2]
                     << " B:" << std::setw(12) << std::setprecision(8) << blh[0] * Deg
                     << " L:" << std::setw(12) << std::setprecision(8) << blh[1] * Deg
                     << " H:" << std::setw(8) << std::setprecision(3) << blh[2]
                     << " Vx:" << std::setw(9) << std::setprecision(4) << Res.Velocity[0]
                     << " Vy:" << std::setw(9) << std::setprecision(4) << Res.Velocity[1]
                     << " Vz:" << std::setw(9) << std::setprecision(4) << Res.Velocity[2]
                     << " GPS Clk:" << std::setw(12) << std::setprecision(3) << Res.RcvClkOft[0] * C_Light
                     << " BDS Clk:" << std::setw(12) << std::setprecision(3) << Res.RcvClkOft[1] * C_Light
                     << " PDOP:" << std::setw(8) << std::setprecision(3) << Res.PDOP
                     << " Sigma:" << std::setw(8) << std::setprecision(3) << Res.SigmaPos
                     << " GPSSats:" << std::setw(3) << (int)Res.GPSSatNum
                     << " BDSSats:" << std::setw(3) << (int)Res.BDSSatNum
                     << " Sats:" << std::setw(3) << (int)Res.AllSatNum << std::endl;

                if (success_epochs % 20 == 0 || success_epochs < 5) {
                    printf("Epoch %d (TOW %.1f): Success. Sats: %d | PDOP: %.2f | H: %.1f | V: [%.3f, %.3f, %.3f]\n", 
                           total_epochs, Res.Time.SecOfWeek, Res.AllSatNum, Res.PDOP, blh[2], Res.Velocity[0], Res.Velocity[1], Res.Velocity[2]);
                }
                // 热启动：将解算的坐标作为下一历元的初值
                Obs.Pos[0] = Res.Position[0];
                Obs.Pos[1] = Res.Position[1];
                Obs.Pos[2] = Res.Position[2];
            } 
            else {
                    printf("Epoch %d (TOW %.1f): SPP Failed. Detected Sats: %d\n", total_epochs, Obs.Time.SecOfWeek, detected_sats);
            }
        }
        
        if (Len == initial_len && !feof(fp) && initial_len > 0) {
            memmove(Buff, Buff + 1, --Len);
        }

    } while (!feof(fp) || Len > 0);

    fclose(fp);
    fout.close();
    printf("\nProcessing Summary:\n");
    printf("Total Epochs Decoded: %d\n", total_epochs);
    printf("Total Success SPP:   %d\n", success_epochs);
    printf("Results saved to result.txt\n");
    return 1;
}
