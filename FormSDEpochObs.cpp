#include "RTK_Structs.h"
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include "lambda.h"

using namespace Eigen;



void FormSDEpochObs(const EPOCHOBS *EpkBas, const EPOCHOBS *EpkRov, SDEPOCHOBS *SDObs){
    if (!EpkBas || !EpkRov || !SDObs) return;

    static double preLockTimeRov[6][64][2] = {0};
    static double preLockTimeBas[6][64][2] = {0};
    static bool   hasPre[6][64] = {false};

    // 清空单差卫星观测数据
    memset(SDObs->SdSatObs, 0, sizeof(SDObs->SdSatObs));
    int sdNum = 0;

    // 根据流动站卫星数进行循环
    for(int i = 0; i < EpkRov->SatNum; i++){
        const SATOBS& satRov = EpkRov->SatObs[i];

        // 1. 检查流动站卫星号和系统号是否正常，双频观测值是否完整
        if (satRov.Prn <= 0 || (satRov.System != GPS && satRov.System != BDS)) continue;
        if (satRov.c1 == 0.0 || satRov.l1 == 0.0 || satRov.p2 == 0.0 || satRov.l2 == 0.0) continue;

        // 2. 在基站观测值中查找与该流动站相同的卫星
        int idxBas = -1;
        for(int j = 0; j < EpkBas->SatNum; j++){
            if (EpkBas->SatObs[j].Prn == satRov.Prn && EpkBas->SatObs[j].System == satRov.System){
                idxBas = j;
                break;
            }
        }
        if (idxBas == -1) continue;

        const SATOBS& satBas = EpkBas->SatObs[idxBas];
        // 检查基站双频观测值是否完整
        if (satBas.c1 == 0.0 || satBas.l1 == 0.0 || satBas.p2 == 0.0 || satBas.l2 == 0.0) continue;

        // 3. 对同类型和同频率的观测值求差并保存索引号
        SDSATOBS& sdSat = SDObs->SdSatObs[sdNum];
        sdSat.Prn = satRov.Prn;
        sdSat.System = satRov.System;
        sdSat.nRov = (short)i;
        sdSat.nBas = (short)idxBas;

        // 站间单差：流动站 - 基准站 (m)
        sdSat.dP[0] = satRov.c1 - satBas.c1;
        sdSat.dL[0] = satRov.l1 - satBas.l1;
        sdSat.dP[1] = satRov.p2 - satBas.p2;
        sdSat.dL[1] = satRov.l2 - satBas.l2;

        // 4. LockTime 检查
        int sys = (int)satRov.System;
        int prn = satRov.Prn;

        if (sys >= 0 && sys < 6 && prn >= 0 && prn < 64) {
            if (hasPre[sys][prn]) {
                bool slipRov = (satRov.LockTime[0] < preLockTimeRov[sys][prn][0] ||
                                satRov.LockTime[1] < preLockTimeRov[sys][prn][1]);
                bool slipBas = (satBas.LockTime[0] < preLockTimeBas[sys][prn][0] ||
                                satBas.LockTime[1] < preLockTimeBas[sys][prn][1]);

                if (slipRov || slipBas) {
                    sdSat.Valid = 0; 
                } else {
                    sdSat.Valid = 1; 
                }
            } else {
                sdSat.Valid = 1; // 新卫星接入
            }

            // 更新静态历史
            preLockTimeRov[sys][prn][0] = satRov.LockTime[0];
            preLockTimeRov[sys][prn][1] = satRov.LockTime[1];
            preLockTimeBas[sys][prn][0] = satBas.LockTime[0];
            preLockTimeBas[sys][prn][1] = satBas.LockTime[1];
            hasPre[sys][prn] = true;
        } else {
            sdSat.Valid = 1;
        }

        sdNum++;
        if (sdNum >= MAXCHANNUM) break;
    }

    // 循环结束后，对单差观测值的观测时刻和卫星数量赋值
    SDObs->Time = EpkRov->Time;
    SDObs->SatNum = (short)sdNum;
}


// 单差周跳探测函数
void DetectCycleSlip(SDEPOCHOBS* Obs){
    if (!Obs) return;

    // 备份上一历元单差组合平滑状态
    MWGF PreSdCObs[MAXCHANNUM];
    memcpy(PreSdCObs, Obs->SdCObs, sizeof(PreSdCObs));

    // 准备当前历元记录
    MWGF CurSdCObs[MAXCHANNUM];
    memset(CurSdCObs, 0, sizeof(CurSdCObs));

    for (int i = 0; i < Obs->SatNum; i++) {
        SDSATOBS& sdSat = Obs->SdSatObs[i];

        CurSdCObs[i].Prn = sdSat.Prn;
        CurSdCObs[i].Sys = sdSat.System;

        // 1. 计算单差组合观测值
        double GF = sdSat.dL[0] - sdSat.dL[1]; // 单差几何无关组合 (m)
        double MW = 0.0;                       // 单差 MW 组合 (周)
        if (sdSat.System == GPS) {
            MW = (FG1_GPS * sdSat.dL[0] - FG2_GPS * sdSat.dL[1]) / (FG1_GPS - FG2_GPS)
               - (FG1_GPS * sdSat.dP[0] + FG2_GPS * sdSat.dP[1]) / (FG1_GPS + FG2_GPS);
        } else if (sdSat.System == BDS) {
            MW = (FG1_BDS * sdSat.dL[0] - FG3_BDS * sdSat.dL[1]) / (FG1_BDS - FG3_BDS)
               - (FG1_BDS * sdSat.dP[0] + FG3_BDS * sdSat.dP[1]) / (FG1_BDS + FG3_BDS);
        }

        // 2. 寻找上一历元单差组合记录
        int preIdx = -1;
        for (int j = 0; j < MAXCHANNUM; j++) {
            if (PreSdCObs[j].Prn == sdSat.Prn && PreSdCObs[j].Sys == sdSat.System) {
                preIdx = j;
                break;
            }
        }

        // 3. 周跳检测与状态平滑递推
        bool slip = (sdSat.Valid == 0); // 前级 LockTime 已标记失锁
        if (preIdx != -1 && PreSdCObs[preIdx].n > 0) {
            double dGF = fabs(GF - PreSdCObs[preIdx].GF);
            double dMW = fabs(MW - PreSdCObs[preIdx].MW);

            if (slip || dGF > 0.05 || dMW > 3.0) {
                // 检测到周跳或粗差：标记无效并以当前历元重置平滑
                sdSat.Valid = 0;
                CurSdCObs[i].MW = MW;
                CurSdCObs[i].n = 1;
            } else {
                // 正常连续无周跳：递推平滑 MW
                sdSat.Valid = 1;
                CurSdCObs[i].MW = (PreSdCObs[preIdx].n * PreSdCObs[preIdx].MW + MW) / (PreSdCObs[preIdx].n + 1);
                CurSdCObs[i].n = PreSdCObs[preIdx].n + 1;
            }
        } else {
            // 新接入共视卫星，初始化新弧段
            sdSat.Valid = 1;
            CurSdCObs[i].MW = MW;
            CurSdCObs[i].n = 1;
        }

        CurSdCObs[i].GF = GF;
    }

    memcpy(Obs->SdCObs, CurSdCObs, sizeof(CurSdCObs));
}

void DetRefSat(const EPOCHOBS* EpkBas, const EPOCHOBS* EpkRov, SDEPOCHOBS* SDObs, DDCOBS* DDObs){
    if (!EpkBas || !EpkRov || !SDObs || !DDObs) return;

    // 1. 初始化基准星与双差卫星数
    DDObs->RefPrn[0] = DDObs->RefPrn[1] = -1;
    DDObs->RefPos[0] = DDObs->RefPos[1] = -1;
    DDObs->DDSatNum[0] = DDObs->DDSatNum[1] = 0;
    DDObs->Sats = 0;

    double max_elev_gps = -1.0;
    double max_elev_bds = -1.0;

    // 2. 遍历单差卫星，寻找各系统高度角最高且质量优良的卫星作为基准星 (必须在两站均具有有效轨道且高度角>=10度)
    for (int i = 0; i < SDObs->SatNum; i++) {
        SDSATOBS& sdobs = SDObs->SdSatObs[i];
        int nRov = sdobs.nRov;
        int nBas = sdobs.nBas;

        if (sdobs.Valid == 0) continue; // 剔除失锁或周跳卫星

        // 检查卫星轨道与 PVT 计算是否在基站和流动站均有效
        if (!EpkRov->SatPVT[nRov].Valid || !EpkBas->SatPVT[nBas].Valid) continue;
        if (Vector3d(EpkRov->SatPVT[nRov].SatPos).norm() < 1e6 ||
            Vector3d(EpkBas->SatPVT[nBas].SatPos).norm() < 1e6) continue;

        // 载噪比质量门限 (双频 > 35.0 dB-Hz，保证基准星观测质量最高)
        if (EpkBas->SatObs[nBas].cn0[0] < 35.0 || EpkBas->SatObs[nBas].cn0[1] < 35.0 ||
            EpkRov->SatObs[nRov].cn0[0] < 35.0 || EpkRov->SatObs[nRov].cn0[1] < 35.0) continue;

        double elev = EpkRov->SatPVT[nRov].Elevation; // 卫星高度角 (rad)
        if (elev < 10.0 * Rad) continue;

        if (sdobs.System == GPS) {
            if (elev > max_elev_gps) {
                max_elev_gps = elev;
                DDObs->RefPrn[0] = sdobs.Prn;
                DDObs->RefPos[0] = i; // 记录基准星在 SDObs->SdSatObs 中的索引号
            }
        } else 
        if (sdobs.System == BDS) {
            if (elev > max_elev_bds) {
                max_elev_bds = elev;
                DDObs->RefPrn[1] = sdobs.Prn;
                DDObs->RefPos[1] = i; // 记录基准星在 SDObs->SdSatObs 中的索引号
            }
        }
    }

    // 3. 统计各系统双差卫星数 (有效单差星中剔除基准星，并要求两站轨道有效且高度角合格)
    int nGps = 0, nBds = 0;
    for (int i = 0; i < SDObs->SatNum; i++) {
        const SDSATOBS& sd = SDObs->SdSatObs[i];
        if (sd.Valid == 0) continue;
        if (!EpkRov->SatPVT[sd.nRov].Valid || !EpkBas->SatPVT[sd.nBas].Valid) continue;
        if (Vector3d(EpkRov->SatPVT[sd.nRov].SatPos).norm() < 1e6 ||
            Vector3d(EpkBas->SatPVT[sd.nBas].SatPos).norm() < 1e6) continue;
        if (EpkRov->SatPVT[sd.nRov].Elevation < 10.0 * Rad) continue;

        if (sd.System == GPS && DDObs->RefPrn[0] != -1 && sd.Prn != DDObs->RefPrn[0]) {
            nGps++;
        } else if (sd.System == BDS && DDObs->RefPrn[1] != -1 && sd.Prn != DDObs->RefPrn[1]) {
            nBds++;
        }
    }

    DDObs->DDSatNum[0] = nGps;
    DDObs->DDSatNum[1] = nBds;
    DDObs->Sats = nGps + nBds;
}


// 构建单频段 nDD 颗卫星的双差权阵子块 (闭合解析解 W_sub = D_DD^{-1})
static MatrixXd BuildSubWeightMatrix(int nDD, double sigma) {
    MatrixXd W_sub(nDD, nDD);
    double sigma2 = sigma * sigma;
    double factor = 1.0 / (2.0 * sigma2 * (nDD + 1));
    for (int r = 0; r < nDD; ++r) {
        for (int c = 0; c < nDD; ++c) {
            W_sub(r, c) = (r == c ? (double)nDD : -1.0) * factor;
        }
    }
    return W_sub;
}

bool RTKFloat(RAWDAT* Raw, PPRESULT* Base, PPRESULT* Rov, const ROVERCFGINFO& cfg){
    if (!Raw || !Base || !Rov) return false;

    // 1. 设置基站真实位置并计算基站卫星轨道与高度角
    static Vector3d s_fixed_XB = Vector3d::Zero();
    static int s_base_src = 0; // 1: BESTPOS in stream, 2: SPP
    static bool s_base_notified = false;//打印标记，确保基站坐标的来源提醒只打印一次

    if (Raw->BasEpk.Pos[0] != 0.0 && Raw->BasEpk.Pos[1] != 0.0 && Raw->BasEpk.Pos[2] != 0.0) {
        // 当前历元数据流中包含基站已知坐标
        // Base->Position != 0 时，IsSuccess == false 表示坐标来自数据流
        Base->Position[0] = Raw->BasEpk.Pos[0];
        Base->Position[1] = Raw->BasEpk.Pos[1];
        Base->Position[2] = Raw->BasEpk.Pos[2];
        Base->IsSuccess = false;
        s_fixed_XB << Base->Position[0], Base->Position[1], Base->Position[2];
        s_base_src = 1;
        if (!s_base_notified) {
            printf("[Notice] Base station coordinate obtained from data stream (BESTPOS): [%.4f, %.4f, %.4f]\n",
                   s_fixed_XB[0], s_fixed_XB[1], s_fixed_XB[2]);
            s_base_notified = true;
        }
    } else if (s_fixed_XB.norm() > 1000.0) {
        // 数据流中无基站坐标，直接沿用历史已确定的坐标
        Base->Position[0] = s_fixed_XB[0];
        Base->Position[1] = s_fixed_XB[1];
        Base->Position[2] = s_fixed_XB[2];
        Base->IsSuccess = (s_base_src == 2); // 保持历史来源属性：SPP为true，数据流已知点为false
    } else {
        // 从未获取过基站坐标，使用 SPP 解算基站初始位置（SPP 内部已将 Base->IsSuccess 设为 true）
        DetectOutlier(&Raw->BasEpk);
        if (!SPP(&Raw->BasEpk, Raw, Base)) return false;
        s_fixed_XB << Base->Position[0], Base->Position[1], Base->Position[2];
        s_base_src = 2;
        if (!s_base_notified) {
            printf("[Notice] Base station coordinate not provided in data stream, initialized via SPP: [%.4f, %.4f, %.4f]\n",
                   s_fixed_XB[0], s_fixed_XB[1], s_fixed_XB[2]);
            s_base_notified = true;
        }
    }
    Vector3d X_B = s_fixed_XB;
    DetectOutlier(&Raw->BasEpk); // 必须每历元计算基站 PIF 组合，供 ComputeGPSSatOrbitAtSignalTrans 计算轨道
    ComputeGPSSatOrbitAtSignalTrans(&Raw->BasEpk, Raw->GpsEph, Raw->BdsEph, Base->Position, Raw->BasEpk.SatPVT);
    for (int i = 0; i < Raw->BasEpk.SatNum; ++i) {
        if (Raw->BasEpk.SatPVT[i].Valid) {
            double el = 0, az = 0;
            CompSatElAz(Base->Position, Raw->BasEpk.SatPVT[i].SatPos, &el, &az);
            Raw->BasEpk.SatPVT[i].Elevation = el;
            Raw->BasEpk.SatPVT[i].Azimuth = az;
        }
    }

    // 2. 流动站执行 SPP 单点定位 (以基站坐标作为初始迭代位置，防止发散)
    if (Raw->RovEpk.Pos[0] == 0.0 || (Vector3d(Raw->RovEpk.Pos) - X_B).norm() > 20000.0) {
        Raw->RovEpk.Pos[0] = Base->Position[0];
        Raw->RovEpk.Pos[1] = Base->Position[1];
        Raw->RovEpk.Pos[2] = Base->Position[2];
    }
    DetectOutlier(&Raw->RovEpk);
    SPP(&Raw->RovEpk, Raw, Rov);

    // 流动站初值：优先使用 SPP 解，若偏差过大则直接使用基站坐标
    Vector3d X_R_star = X_B;
    if (Rov->IsSuccess && (Vector3d(Rov->Position) - X_B).norm() < 50000.0) {
        X_R_star << Rov->Position[0], Rov->Position[1], Rov->Position[2];
    }

    // 确保流动站卫星轨道和高度角基于可靠几何位置稳定精化 (避免 SPP 发散污染)
    ComputeGPSSatOrbitAtSignalTrans(&Raw->RovEpk, Raw->GpsEph, Raw->BdsEph, X_R_star.data(), Raw->RovEpk.SatPVT);
    for (int i = 0; i < Raw->RovEpk.SatNum; ++i) {
        if (Raw->RovEpk.SatPVT[i].Valid) {
            double el = 0, az = 0;
            CompSatElAz(X_R_star.data(), Raw->RovEpk.SatPVT[i].SatPos, &el, &az);
            Raw->RovEpk.SatPVT[i].Elevation = el;
            Raw->RovEpk.SatPVT[i].Azimuth = az;
        }
    }

    // 3. 选取基准星并确定各系统双差卫星数
    DetRefSat(&Raw->BasEpk, &Raw->RovEpk, &Raw->SdObs, &Raw->DDObs);
    int m_DD = Raw->DDObs.DDSatNum[0]; // GPS 双差星数
    int n_DD = Raw->DDObs.DDSatNum[1]; // BDS 双差星数
    int S    = m_DD + n_DD;           // 双差总星数
    int n    = 3 + 2 * S;             // 待估参数总数 (3坐标改正数 + 2S双差模糊度)
    int l    = 4 * S;                 // 双差观测方程总数 (双系统双频双差)

    if (l <= n || S < 3) return false; // 多余观测数 r = l - n 必须大于 0
    if (m_DD > 0 && Raw->DDObs.RefPos[0] < 0) return false;//没有找到参考星
    if (n_DD > 0 && Raw->DDObs.RefPos[1] < 0) return false;
    int ref_idx[2] = { Raw->DDObs.RefPos[0], Raw->DDObs.RefPos[1] };

    // 4. 计算基站到各可用卫星的固定几何距离 d_B
    std::vector<double> d_B(Raw->SdObs.SatNum, 0.0);
    for (int i = 0; i < Raw->SdObs.SatNum; ++i) {
        int nBas = Raw->SdObs.SdSatObs[i].nBas;
        const double* sat_pos = Raw->BasEpk.SatPVT[nBas].SatPos;
        double dx = sat_pos[0] - X_B[0], dy = sat_pos[1] - X_B[1], dz = sat_pos[2] - X_B[2];
        d_B[i] = sqrt(dx * dx + dy * dy + dz * dz);
    }

    // 5. 组装全局双差权矩阵 W (准对角分块矩阵)
    MatrixXd W = MatrixXd::Zero(l, l);
    double sigP = (cfg.CodeNoise > 0.0) ? cfg.CodeNoise : 0.3;
    double sigL = (cfg.CPNoise > 0.0)   ? cfg.CPNoise   : 0.003;
    if (m_DD > 0) {
        MatrixXd W_PG1 = BuildSubWeightMatrix(m_DD, sigP);
        MatrixXd W_LG1 = BuildSubWeightMatrix(m_DD, sigL);
        MatrixXd W_PG2 = BuildSubWeightMatrix(m_DD, sigP);
        MatrixXd W_LG2 = BuildSubWeightMatrix(m_DD, sigL);
        W.block(0 * m_DD, 0 * m_DD, m_DD, m_DD) = W_PG1;
        W.block(1 * m_DD, 1 * m_DD, m_DD, m_DD) = W_LG1;
        W.block(2 * m_DD, 2 * m_DD, m_DD, m_DD) = W_PG2;
        W.block(3 * m_DD, 3 * m_DD, m_DD, m_DD) = W_LG2;
    }
    if (n_DD > 0) {
        int offset = 4 * m_DD;
        MatrixXd W_PB1 = BuildSubWeightMatrix(n_DD, sigP);
        MatrixXd W_LB1 = BuildSubWeightMatrix(n_DD, sigL);
        MatrixXd W_PB3 = BuildSubWeightMatrix(n_DD, sigP);
        MatrixXd W_LB3 = BuildSubWeightMatrix(n_DD, sigL);
        W.block(offset + 0 * n_DD, offset + 0 * n_DD, n_DD, n_DD) = W_PB1;
        W.block(offset + 1 * n_DD, offset + 1 * n_DD, n_DD, n_DD) = W_LB1;
        W.block(offset + 2 * n_DD, offset + 2 * n_DD, n_DD, n_DD) = W_PB3;
        W.block(offset + 3 * n_DD, offset + 3 * n_DD, n_DD, n_DD) = W_LB3;
    }

    // 6. 高斯-牛顿最小二乘迭代循环
    bool converged = false;
    MatrixXd Qxx(n, n);
    VectorXd x_hat(n);
    MatrixXd H = MatrixXd::Zero(l, n);
    VectorXd z = VectorXd::Zero(l);

    for (int iter = 0; iter < 10; ++iter) {
        // A. 计算流动站到各系统参考星的几何距离与方向余弦 (l, m, n)
        double d_R_ref[2] = {0}, l_ref[2] = {0}, m_ref[2] = {0}, n_ref[2] = {0};
        for (int sys = 0; sys < 2; ++sys) {
            if (ref_idx[sys] == -1) continue;
            int nRov = Raw->SdObs.SdSatObs[ref_idx[sys]].nRov;
            const double* spos = Raw->RovEpk.SatPVT[nRov].SatPos;
            double dx = spos[0] - X_R_star[0], dy = spos[1] - X_R_star[1], dz = spos[2] - X_R_star[2];
            d_R_ref[sys] = sqrt(dx * dx + dy * dy + dz * dz);
            if (d_R_ref[sys] > 0.0) {
                l_ref[sys] = dx / d_R_ref[sys];
                m_ref[sys] = dy / d_R_ref[sys];
                n_ref[sys] = dz / d_R_ref[sys];
            }
        }

        // B. 初始化设计矩阵 H 与 OMC 观测向量 z (对应模型 z = H * x + Delta)
        H.setZero();
        z.setZero();
        int k_G = 0, k_C = 0;

        for (int i = 0; i < Raw->SdObs.SatNum; ++i) {
            const SDSATOBS& sd = Raw->SdObs.SdSatObs[i];
            if (sd.Valid == 0) continue;
            if (!Raw->RovEpk.SatPVT[sd.nRov].Valid || !Raw->BasEpk.SatPVT[sd.nBas].Valid) continue;
            if (Vector3d(Raw->BasEpk.SatPVT[sd.nBas].SatPos).norm() < 1e6 ||
                Vector3d(Raw->RovEpk.SatPVT[sd.nRov].SatPos).norm() < 1e6) continue;
            double elev_mask = (cfg.ElevThreshold > 0.0) ? cfg.ElevThreshold * Rad : 10.0 * Rad;
            if (Raw->RovEpk.SatPVT[sd.nRov].Elevation < elev_mask) continue;

            int sys = (sd.System == GPS) ? 0 : 1;
            if (ref_idx[sys] == -1) continue;
            if (i == ref_idx[sys]) continue; // 跳过参考星

            int nRov = sd.nRov;
            const double* spos = Raw->RovEpk.SatPVT[nRov].SatPos;
            double dx = spos[0] - X_R_star[0], dy = spos[1] - X_R_star[1], dz = spos[2] - X_R_star[2];
            double d_R_j = sqrt(dx * dx + dy * dy + dz * dz);
            if (d_R_j <= 0.0) continue;
            double l_j = dx / d_R_j, m_j = dy / d_R_j, n_j = dz / d_R_j;

            // 方向余弦差 (待估星 - 参考星)
            double dl = l_j - l_ref[sys];
            double dm = m_j - m_ref[sys];
            double dn = n_j - n_ref[sys];

            // 双差几何距离近似计算量 rho_DD_star
            double rho_DD_star = (d_R_j - d_R_ref[sys]) - (d_B[i] - d_B[ref_idx[sys]]);

            const SDSATOBS& sd_ref = Raw->SdObs.SdSatObs[ref_idx[sys]];
            if (Raw->RovEpk.Time.SecOfWeek >= 116706.5 && Raw->RovEpk.Time.SecOfWeek <= 116708.5 && iter == 0) {
                if (std::isnan(rho_DD_star) || std::isnan(sd.dP[0]) || std::isnan(sd_ref.dP[0]) ||
                    std::isnan(sd.dL[0]) || std::isnan(sd_ref.dL[0])) {
                    printf("  DEBUG NAN Sat %c%02d: d_R_j=%.2f, d_R_ref=%.2f, d_B[i]=%.2f, d_B[ref]=%.2f, dP[0]=%.2f, ref_dP[0]=%.2f\n",
                           sd.System == GPS ? 'G' : 'C', sd.Prn, d_R_j, d_R_ref[sys], d_B[i], d_B[ref_idx[sys]],
                           sd.dP[0], sd_ref.dP[0]);
                }
            }
            if (sd.System == GPS && m_DD > 0 && k_G < m_DD) {
                // GPS L1 伪距与载波
                z(0 * m_DD + k_G) = (sd.dP[0] - sd_ref.dP[0]) - rho_DD_star;
                H.block<1, 3>(0 * m_DD + k_G, 0) << -dl, -dm, -dn;

                z(1 * m_DD + k_G) = (sd.dL[0] - sd_ref.dL[0]) - rho_DD_star;
                H.block<1, 3>(1 * m_DD + k_G, 0) << -dl, -dm, -dn;
                H(1 * m_DD + k_G, 3 + k_G) = WL1_GPS;

                // GPS L2 伪距与载波
                z(2 * m_DD + k_G) = (sd.dP[1] - sd_ref.dP[1]) - rho_DD_star;
                H.block<1, 3>(2 * m_DD + k_G, 0) << -dl, -dm, -dn;

                z(3 * m_DD + k_G) = (sd.dL[1] - sd_ref.dL[1]) - rho_DD_star;
                H.block<1, 3>(3 * m_DD + k_G, 0) << -dl, -dm, -dn;
                H(3 * m_DD + k_G, 3 + m_DD + k_G) = WL2_GPS;

                k_G++;
            } else if (sd.System == BDS && n_DD > 0 && k_C < n_DD) {
                int off = 4 * m_DD;
                // BDS B1 伪距与载波
                z(off + 0 * n_DD + k_C) = (sd.dP[0] - sd_ref.dP[0]) - rho_DD_star;
                H.block<1, 3>(off + 0 * n_DD + k_C, 0) << -dl, -dm, -dn;

                z(off + 1 * n_DD + k_C) = (sd.dL[0] - sd_ref.dL[0]) - rho_DD_star;
                H.block<1, 3>(off + 1 * n_DD + k_C, 0) << -dl, -dm, -dn;
                H(off + 1 * n_DD + k_C, 3 + 2 * m_DD + k_C) = WL1_BDS;

                // BDS B3 伪距与载波
                z(off + 2 * n_DD + k_C) = (sd.dP[1] - sd_ref.dP[1]) - rho_DD_star;
                H.block<1, 3>(off + 2 * n_DD + k_C, 0) << -dl, -dm, -dn;

                z(off + 3 * n_DD + k_C) = (sd.dL[1] - sd_ref.dL[1]) - rho_DD_star;
                H.block<1, 3>(off + 3 * n_DD + k_C, 0) << -dl, -dm, -dn;
                H(off + 3 * n_DD + k_C, 3 + 2 * m_DD + n_DD + k_C) = WL3_BDS;

                k_C++;
            }
        }

        // C. 最小二乘求解法方程: (H^T * W * H) * x_hat = H^T * W * z
        MatrixXd HtWH = H.transpose() * W * H;
        VectorXd HtWz = H.transpose() * W * z;

        Eigen::LDLT<MatrixXd> ldlt(HtWH);
        if (ldlt.info() == Eigen::Success) {
            x_hat = ldlt.solve(HtWz);
            Qxx = HtWH.inverse();
        } else {
            Qxx = HtWH.inverse();
            x_hat = Qxx * HtWz;
        }

        // D. 更新流动站近似坐标: X_R_star = X_R_star + x_hat[0..2]
        X_R_star[0] += x_hat(0);
        X_R_star[1] += x_hat(1);
        X_R_star[2] += x_hat(2);

        // E. 收敛判定 (坐标改正数模长 < 0.1 mm)
        double dx_norm = x_hat.head<3>().norm();
        if (Raw->RovEpk.Time.SecOfWeek >= 116706.5 && Raw->RovEpk.Time.SecOfWeek <= 116708.5) {
            printf("  [TOW %.1f Iter %d] dx_norm: %e, x_hat[0..2]: [%.4f, %.4f, %.4f], k_G=%d/%d, k_C=%d/%d\n",
                   Raw->RovEpk.Time.SecOfWeek, iter, dx_norm, x_hat(0), x_hat(1), x_hat(2),
                   k_G, m_DD, k_C, n_DD);
            printf("    X_R_star: [%.3f, %.3f, %.3f], d_R_ref: [%.3f, %.3f]\n",
                   X_R_star[0], X_R_star[1], X_R_star[2], d_R_ref[0], d_R_ref[1]);
            if (iter == 0) {
                for (int k = 0; k < l; ++k) {
                    if (std::isnan(z(k)) || std::isinf(z(k))) {
                        printf("    z(%d) is NaN/Inf! l=%d\n", k, l);
                    }
                }
            }
        }
        if (dx_norm < 1.0e-4) {
            converged = true;
            break;
        }
    }

    if (!converged) {
        if (Raw->RovEpk.Time.SecOfWeek >= 116706.5 && Raw->RovEpk.Time.SecOfWeek <= 116708.5) {
            printf("  [TOW %.1f FAIL] GN did not converge!\n", Raw->RovEpk.Time.SecOfWeek);
        }
        return false;
    }

    // 7. 解算成功，更新流动站与基线结果
    Rov->Position[0] = X_R_star[0];
    Rov->Position[1] = X_R_star[1];
    Rov->Position[2] = X_R_star[2];
    Rov->IsSuccess = true;
    Rov->Time = Raw->RovEpk.Time;

    Raw->DDObs.dPos[0] = X_R_star[0] - X_B[0];
    Raw->DDObs.dPos[1] = X_R_star[1] - X_B[1];
    Raw->DDObs.dPos[2] = X_R_star[2] - X_B[2];

    // 残差向量与验后单位权中误差
    VectorXd v = z - H * x_hat;
    double vTWv = (v.transpose() * W * v)(0);
    double sigma0 = (l > n) ? sqrt(std::max(0.0, vTWv / (l - n))) : 0.0;
    Rov->SigmaPos = sigma0;

    // 8. 提取浮点参数并调用 LAMBDA 执行模糊度搜索与固定基线恢复
    Vector3d b_hat = X_R_star - X_B;
    VectorXd a_hat = x_hat.tail(2 * S);
    Matrix3d Qbb = Qxx.block<3, 3>(0, 0);
    MatrixXd Qaa = Qxx.block(3, 3, 2 * S, 2 * S);
    MatrixXd Qba = Qxx.block(0, 3, 3, 2 * S);

    FixAmbiguity(Raw, Base, Rov, b_hat, a_hat, Qbb, Qaa, Qba, cfg.RatioThres);

    return true;
}