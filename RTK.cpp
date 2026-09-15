#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include "RTK_Structs.h"
#include "sockets.cpp"
#include "ReadConfig.cpp"

/* ============================================================
 * SPP 模式: 单路数据源 (流动站), 逐历元单点定位 + 单点测速
 * ============================================================ */
static void RunSPP(const ROVERCFGINFO& cfg, std::ofstream& fout)
{
    FILE* fp = NULL;
    SOCKET NetGps = INVALID_SOCKET;
    unsigned char Buff[MAXRAWLEN];
    int Len = 0, LenRead, val;

    if (cfg.IsFileData == 0) {
        if (OpenSocket(NetGps, cfg.RovNetIP, (unsigned short)cfg.RovNetPort) == false) {
            printf("Failed to connect to %s:%d\n", cfg.RovNetIP, cfg.RovNetPort);
            return;
        }
    } else {
        if ((fp = fopen(cfg.RovObsDatFile, "rb")) == NULL) {
            printf("The file '%s' was not opened\n", cfg.RovObsDatFile);
            return;
        }
    }

    EPOCHOBS Obs;
    memset(&Obs, 0, sizeof(Obs));
    GPSEPHREC GpsEph[MAXGPSNUM], BdsEph[MAXBDSNUM];

    memset(GpsEph, 0, sizeof(GpsEph));
    memset(BdsEph, 0, sizeof(BdsEph));

    printf("Starting SPP Debugging Flow...\n");

    int total_epochs = 0;
    int success_epochs = 0;

    do {
        if (cfg.IsFileData == 1) {
            LenRead = fread(Buff + Len, sizeof(unsigned char), MAXRAWLEN - Len, fp);
            if (LenRead == 0) break;
            Len = Len + LenRead;
        } else {
            Sleep(980);
            LenRead = recv(NetGps, (char*)(Buff + Len), MAXRAWLEN - Len, 0);
            if (LenRead <= 0) break;
            Len += LenRead;
        }


        while (true) {
            val = DecodeNovOem7Dat(Buff, Len, &Obs, GpsEph, BdsEph);
            if (val == -1) {
                break;
            }

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

                    fout << "SPP: Epoch " << total_epochs << " " << (int)Res.Time.Week << " " << std::fixed << std::setprecision(3) << Res.Time.SecOfWeek
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

                    printf("Epoch %d (TOW %.1f): Success. Sats: %d | PDOP: %.2f | XYZ: [%.3f, %.3f, %.3f] | H: %.1f | V: [%.3f, %.3f, %.3f]\n",
                           total_epochs, Res.Time.SecOfWeek, Res.AllSatNum, Res.PDOP,
                           Res.Position[0], Res.Position[1], Res.Position[2],
                           blh[2], Res.Velocity[0], Res.Velocity[1], Res.Velocity[2]);
                    // 将解算的坐标作为下一历元的初值
                    Obs.Pos[0] = Res.Position[0];
                    Obs.Pos[1] = Res.Position[1];
                    Obs.Pos[2] = Res.Position[2];
                }
                else {
                    printf("Epoch %d (TOW %.1f): SPP Failed. Detected Sats: %d\n", total_epochs, Obs.Time.SecOfWeek, detected_sats);
                }
            }
        }

    } while (true);

    if (cfg.IsFileData == 1) {
        if (fp) fclose(fp);
    } else {
        if (NetGps != INVALID_SOCKET) CloseSocket(NetGps);
    }

    printf("\nProcessing Summary:\n");
    printf("Total Epochs Decoded: %d\n", total_epochs);
    printf("Total Success SPP:   %d\n", success_epochs);
    printf("Results saved to %s\n", cfg.ResFile);
}

/* ============================================================
 * RTK 模式: 基准站/流动站双路数据, 历元同步后做单差/双差解算
 * ============================================================ */
static void RunRTK(const ROVERCFGINFO& cfg, std::ofstream& fout)
{
    const double SynThres = 0.3;   // 历元同步阈值 [s], 后续可移入 config.ini

    FILE *FBas = NULL, *FRov = NULL;
    SOCKET BasSock = INVALID_SOCKET, RovSock = INVALID_SOCKET;

    if (cfg.IsFileData == 0) {
        if (!OpenSocket(RovSock, cfg.RovNetIP, (unsigned short)cfg.RovNetPort) ||
            !OpenSocket(BasSock, cfg.BasNetIP, (unsigned short)cfg.BasNetPort)) {
            printf("Failed to connect bas/rov network streams\n");
            return;
        }
    } else {
        if ((FRov = fopen(cfg.RovObsDatFile, "rb")) == NULL) {
            printf("The file '%s' was not opened\n", cfg.RovObsDatFile);
            return;
        }
        if ((FBas = fopen(cfg.BasObsDatFile, "rb")) == NULL) {
            printf("The file '%s' was not opened\n", cfg.BasObsDatFile);
            return;
        }
    }

    RAWDAT Raw;

    printf("Starting RTK Debugging Flow...\n");

    int total_epochs = 0;

    while (true) {
        int syn = GetSynObs(FBas, FRov, BasSock, RovSock, &Raw, cfg.IsFileData, SynThres);
        if (syn == 0) break;        // 数据结束
        if (syn < 0)  continue;     // 尚未同步, 继续读下一历元

        total_epochs++;

        // TODO: 单差 -> 周跳探测 -> 基准星选取 -> 双差 -> RTK 浮点解/固定解
        // FormSDEpochObs(&Raw.BasEpk, &Raw.RovEpk, &Raw.SdObs);
        // DetectCycleSlip(&Raw.SdObs);
        // DetRefSat(&Raw.BasEpk, &Raw.RovEpk, &Raw.SdObs, &Raw.DDObs);
        // RTKFloat(&Raw, &BasRes, &RovRes);

        printf("Sync epoch %d: TOW %.3f, Bas %d sats, Rov %d sats\n",
               total_epochs, Raw.RovEpk.Time.SecOfWeek,
               Raw.BasEpk.SatNum, Raw.RovEpk.SatNum);
    }

    if (cfg.IsFileData == 1) {
        if (FBas) fclose(FBas);
        if (FRov) fclose(FRov);
    } else {
        if (BasSock != INVALID_SOCKET) CloseSocket(BasSock);
        if (RovSock != INVALID_SOCKET) CloseSocket(RovSock);
    }

    printf("\nRTK Processing Summary:\n");
    printf("Total Synced Epochs: %d\n", total_epochs);
    printf("Results saved to %s\n", cfg.ResFile);
}

int main(int argc, char* argv[])
{
    const char* cfgPath = (argc > 1) ? argv[1] : "config.ini";

    ROVERCFGINFO cfg;
    // 字符串字段缺省值（读取配置成功后会被覆盖）
    strcpy(cfg.RovObsDatFile, "./data/oem719-202603111200.bin");
    strcpy(cfg.BasObsDatFile, "./data/oem719-202603111200.bin");
    strcpy(cfg.ResFile, "./result.txt");
    strcpy(cfg.RovNetIP, "47.114.134.129");
    strcpy(cfg.BasNetIP, "47.114.134.129");
    cfg.RovNetPort = 7190;
    cfg.BasNetPort = 7190;

    if (!ReadRTKConfigInfo(cfgPath, cfg))
        printf("Warning: cannot open '%s', using default parameters.\n", cfgPath);

    std::ofstream fout(cfg.ResFile);
    if (!fout.is_open()) {
        printf("Failed to open %s for writing\n", cfg.ResFile);
    }

    // 按解算模式分流: 1 = RTK 载波相位差分, 0 = SPP 单点定位
    if (cfg.IsRTKMode)
        RunRTK(cfg, fout);
    else
        RunSPP(cfg, fout);

    fout.close();
    return 1;
}
