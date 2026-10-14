#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "RTK_Structs.h"
#include "lambda.h"

using namespace Eigen;

/*------------------------------------------------------------------------------
* RTKFilter.cpp : 模块化扩展卡尔曼滤波 (EKF) 相对定位算法
*
* 架构规范：
*   1. 状态与通道直接利用老师预留的 RTKEKF 结构体：
*      - Ekf->Index[k] : 记录第 k 个双差通道绑定的卫星 PRN 号
*      - Ekf->nPos[k]  : 记录该通道在当前历元单差数组 SdObs 中的真实下标 i
*      - Ekf->nSats    : 当前活跃的双差卫星总数 (m_DD + n_DD)
*   2. 状态向量 X 与协方差 P 严格对齐：
*      - X[0..2]       : 流动站三维坐标 (X_R, Y_R, Z_R) [m]
*      - X[3 + 2*k]    : 通道 k 卫星频段 1 双差整周模糊度 [cycles]
*      - X[3 + 2*k + 1]: 通道 k 卫星频段 2 双差整周模糊度 [cycles]
*   3. 模块化三大独立子函数：
*      - InitEKF()       : 滤波器首历元状态初始化
*      - TimeUpdateEKF() : 时间预测 (状态与协方差外推、通道匹配、周跳/新星重置)
*      - MeasUpdateEKF() : 量测更新 (构建 H, y, R，计算增益 K，Joseph 稳定方差更新)
*-----------------------------------------------------------------------------*/

#define TOTAL_STATE_DIM (3 + MAXCHANNUM * 2)

// 构建单系统双差完全相关量测协方差矩阵 R_sub = W_sub^{-1} (维数: nDD x nDD)
static MatrixXd BuildRSubMatrix(int nDD, double sigma)
{
    MatrixXd R_sub(nDD, nDD);
    double sig2 = sigma * sigma;
    for (int r = 0; r < nDD; ++r) {
        for (int c = 0; c < nDD; ++c) {
            R_sub(r, c) = (r == c) ? (2.0 * sig2) : (sig2);
        }
    }
    return R_sub;
}

/* =============================================================================
 * 子模块 1: 滤波器首历元状态初始化 (InitEKF)
 * ============================================================================= */
static void InitEKF(RTKEKF* Ekf, const Vector3d& X_R_init)
{
    if (!Ekf) return;
    Ekf->IsInit = false;
    Ekf->nSats = 0;
    for (int i = 0; i < MAXCHANNUM; ++i) {
        Ekf->Index[i] = -1;
        Ekf->nPos[i] = -1;
        Ekf->FixAmb[i] = -1;
    }
    for (int i = 0; i < TOTAL_STATE_DIM; ++i) {
        Ekf->X[i] = Ekf->X0[i] = 0.0;
        for (int j = 0; j < TOTAL_STATE_DIM; ++j) {
            Ekf->P[i * TOTAL_STATE_DIM + j] = Ekf->P0[i * TOTAL_STATE_DIM + j] = 0.0;
        }
    }

    // 初始位置设定 (以 SPP 解算坐标初始化)
    Ekf->X[0] = X_R_init[0];
    Ekf->X[1] = X_R_init[1];
    Ekf->X[2] = X_R_init[2];

    // 坐标初始方差 (10.0 m)^2 = 100 m^2
    double var_pos_init = 100.0;
    Ekf->P[0 * TOTAL_STATE_DIM + 0] = var_pos_init;
    Ekf->P[1 * TOTAL_STATE_DIM + 1] = var_pos_init;
    Ekf->P[2 * TOTAL_STATE_DIM + 2] = var_pos_init;
}

/* =============================================================================
 * 子模块 2: 时间预测 (TimeUpdateEKF)
 * 功能：
 *   1. 坐标状态外推，累加过程噪声；
 *   2. 通道直接对齐：GPS 占 [0 .. m_DD-1]，BDS 占 [m_DD .. m_DD+n_DD-1]；
 *   3. 连续无周跳星继承历史模糊度；新升星/周跳星/参考星切换时重置初值与方差；
 *   4. 清理闲置通道。
 * ============================================================================= */
static void TimeUpdateEKF(RTKEKF* Ekf, const RAWDAT* Raw, const int ref_idx[2],
                          int m_DD, int n_DD, const bool ref_changed[2], double q_pos)
{
    // 1. 坐标状态时间预测 (静态短基线模式注入极小过程噪声防止无约束漂移)
    Ekf->P[0 * TOTAL_STATE_DIM + 0] += q_pos;
    Ekf->P[1 * TOTAL_STATE_DIM + 1] += q_pos;
    Ekf->P[2 * TOTAL_STATE_DIM + 2] += q_pos;

    double init_var = 25.0; // 模糊度先验方差 [cycles^2]
    double q_amb = 1.0e-10;  // 连续跟踪无周跳时的微小过程噪声 [cycles^2]
    int S = m_DD + n_DD;

    // 2. 遍历 GPS 通道 [0 .. m_DD - 1]
    for (int k = 0; k < m_DD; ++k) {
        int i = Ekf->nPos[k]; // 单差数组下标
        const SDSATOBS& sd = Raw->SdObs.SdSatObs[i];
        const SDSATOBS& sd_ref = Raw->SdObs.SdSatObs[ref_idx[0]];
        int prn = sd.Prn;
        bool reset = ref_changed[0] || (sd.Valid == 0) || (Ekf->Index[k] != prn);

        int idx1 = 3 + 2 * k;
        int idx2 = 3 + 2 * k + 1;

        if (!reset) {
            // 连续跟踪无周跳：保持模糊度不变，微弱累加过程噪声
            Ekf->P[idx1 * TOTAL_STATE_DIM + idx1] += q_amb;
            Ekf->P[idx2 * TOTAL_STATE_DIM + idx2] += q_amb;
        } else {
            // 新星、周跳或参考星改变：重新赋初值并重置方差
            Ekf->Index[k] = prn;
            double dd_P1 = sd.dP[0] - sd_ref.dP[0];
            double dd_L1 = sd.dL[0] - sd_ref.dL[0];
            double dd_P2 = sd.dP[1] - sd_ref.dP[1];
            double dd_L2 = sd.dL[1] - sd_ref.dL[1];
            Ekf->X[idx1] = (dd_P1 - dd_L1) / WL1_GPS;
            Ekf->X[idx2] = (dd_P2 - dd_L2) / WL2_GPS;

            for (int m = 0; m < TOTAL_STATE_DIM; ++m) {
                Ekf->P[idx1 * TOTAL_STATE_DIM + m] = Ekf->P[m * TOTAL_STATE_DIM + idx1] = 0.0;
                Ekf->P[idx2 * TOTAL_STATE_DIM + m] = Ekf->P[m * TOTAL_STATE_DIM + idx2] = 0.0;
            }
            Ekf->P[idx1 * TOTAL_STATE_DIM + idx1] = init_var;
            Ekf->P[idx2 * TOTAL_STATE_DIM + idx2] = init_var;
        }
    }

    // 3. 遍历 BDS 通道 [m_DD .. S - 1]
    for (int k_C = 0; k_C < n_DD; ++k_C) {
        int k = m_DD + k_C;
        int i = Ekf->nPos[k];
        const SDSATOBS& sd = Raw->SdObs.SdSatObs[i];
        const SDSATOBS& sd_ref = Raw->SdObs.SdSatObs[ref_idx[1]];
        int prn = sd.Prn;
        bool reset = ref_changed[1] || (sd.Valid == 0) || (Ekf->Index[k] != prn);

        int idx1 = 3 + 2 * k;
        int idx2 = 3 + 2 * k + 1;

        if (!reset) {
            Ekf->P[idx1 * TOTAL_STATE_DIM + idx1] += q_amb;
            Ekf->P[idx2 * TOTAL_STATE_DIM + idx2] += q_amb;
        } else {
            Ekf->Index[k] = prn;
            double dd_P1 = sd.dP[0] - sd_ref.dP[0];
            double dd_L1 = sd.dL[0] - sd_ref.dL[0];
            double dd_P2 = sd.dP[1] - sd_ref.dP[1];
            double dd_L2 = sd.dL[1] - sd_ref.dL[1];
            Ekf->X[idx1] = (dd_P1 - dd_L1) / WL1_BDS;
            Ekf->X[idx2] = (dd_P2 - dd_L2) / WL3_BDS;

            for (int m = 0; m < TOTAL_STATE_DIM; ++m) {
                Ekf->P[idx1 * TOTAL_STATE_DIM + m] = Ekf->P[m * TOTAL_STATE_DIM + idx1] = 0.0;
                Ekf->P[idx2 * TOTAL_STATE_DIM + m] = Ekf->P[m * TOTAL_STATE_DIM + idx2] = 0.0;
            }
            Ekf->P[idx1 * TOTAL_STATE_DIM + idx1] = init_var;
            Ekf->P[idx2 * TOTAL_STATE_DIM + idx2] = init_var;
        }
    }

    // 4. 清理当前历元闲置的槽位 [S .. MAXCHANNUM - 1]
    for (int k = S; k < MAXCHANNUM; ++k) {
        Ekf->Index[k] = -1;
        Ekf->nPos[k] = -1;
        int idx1 = 3 + 2 * k;
        int idx2 = 3 + 2 * k + 1;
        Ekf->X[idx1] = Ekf->X[idx2] = 0.0;
        for (int m = 0; m < TOTAL_STATE_DIM; ++m) {
            Ekf->P[idx1 * TOTAL_STATE_DIM + m] = Ekf->P[m * TOTAL_STATE_DIM + idx1] = 0.0;
            Ekf->P[idx2 * TOTAL_STATE_DIM + m] = Ekf->P[m * TOTAL_STATE_DIM + idx2] = 0.0;
        }
    }
}

/* =============================================================================
 * 子模块 3: 量测更新 (MeasUpdateEKF)
 * 功能：
 *   1. 直接提取活跃状态 X[0..3+2S-1] 与子协方差 P[0..3+2S-1, 0..3+2S-1]；
 *   2. 构建设计矩阵 H 与 OMC 残差向量 y；
 *   3. 构建完全相关的量测噪声矩阵 R；
 *   4. 计算增益 K 并更新状态增量；
 *   5. Joseph 稳定形式更新协方差矩阵；
 *   6. 写回结构体 Ekf 并输出基线结果与精度指标。
 * ============================================================================= */
static bool MeasUpdateEKF(RTKEKF* Ekf, RAWDAT* Raw, PPRESULT* Base, PPRESULT* Rov,
                          const std::vector<double>& d_B,
                          const int ref_idx[2],
                          int m_DD, int n_DD, double sigP, double sigL, double ratioThres)
{
    int S = m_DD + n_DD;
    int l = 4 * S;
    int n_act = 3 + 2 * S;

    // 1. 活跃状态向量与协方差矩阵直接从 Ekf 头部切片提取 (零中转、零额外映射)
    VectorXd X_act(n_act);
    MatrixXd P_act(n_act, n_act);
    for (int r = 0; r < n_act; ++r) {
        X_act(r) = Ekf->X[r];
        for (int c = 0; c < n_act; ++c) {
            P_act(r, c) = Ekf->P[r * TOTAL_STATE_DIM + c];
        }
    }

    // 2. 构造设计矩阵 H 与 OMC 残差向量 y
    Vector3d X_R_est(X_act(0), X_act(1), X_act(2));

    double d_R_ref[2] = {0}, l_ref[2] = {0}, m_ref[2] = {0}, n_ref[2] = {0};
    for (int sys = 0; sys < 2; ++sys) {
        if (ref_idx[sys] == -1) continue;
        int nRov = Raw->SdObs.SdSatObs[ref_idx[sys]].nRov;
        const double* spos = Raw->RovEpk.SatPVT[nRov].SatPos;
        double dx = spos[0] - X_R_est[0], dy = spos[1] - X_R_est[1], dz = spos[2] - X_R_est[2];
        d_R_ref[sys] = sqrt(dx * dx + dy * dy + dz * dz);
        l_ref[sys] = dx / d_R_ref[sys];
        m_ref[sys] = dy / d_R_ref[sys];
        n_ref[sys] = dz / d_R_ref[sys];
    }

    MatrixXd H = MatrixXd::Zero(l, n_act);
    VectorXd y = VectorXd::Zero(l);

    // 填充 GPS 观测方程
    if (m_DD > 0 && ref_idx[0] >= 0) {
        const SDSATOBS& sd_ref = Raw->SdObs.SdSatObs[ref_idx[0]];
        for (int k = 0; k < m_DD; ++k) {
            int i = Ekf->nPos[k];
            const SDSATOBS& sd = Raw->SdObs.SdSatObs[i];
            int nRov = sd.nRov;
            const double* spos = Raw->RovEpk.SatPVT[nRov].SatPos;
            double dx = spos[0] - X_R_est[0], dy = spos[1] - X_R_est[1], dz = spos[2] - X_R_est[2];
            double d_R_j = sqrt(dx * dx + dy * dy + dz * dz);
            double dl = (dx / d_R_j) - l_ref[0];
            double dm = (dy / d_R_j) - m_ref[0];
            double dn = (dz / d_R_j) - n_ref[0];
            double rho_DD = (d_R_j - d_R_ref[0]) - (d_B[i] - d_B[ref_idx[0]]);

            int col1 = 3 + 2 * k;
            int col2 = 3 + 2 * k + 1;
            double amb1 = X_act(col1);
            double amb2 = X_act(col2);

            // L1 伪距与载波
            y(0 * m_DD + k) = (sd.dP[0] - sd_ref.dP[0]) - rho_DD;
            H.block<1, 3>(0 * m_DD + k, 0) << -dl, -dm, -dn;

            y(1 * m_DD + k) = (sd.dL[0] - sd_ref.dL[0]) - (rho_DD + WL1_GPS * amb1);
            H.block<1, 3>(1 * m_DD + k, 0) << -dl, -dm, -dn;
            H(1 * m_DD + k, col1) = WL1_GPS;

            // L2 伪距与载波
            y(2 * m_DD + k) = (sd.dP[1] - sd_ref.dP[1]) - rho_DD;
            H.block<1, 3>(2 * m_DD + k, 0) << -dl, -dm, -dn;

            y(3 * m_DD + k) = (sd.dL[1] - sd_ref.dL[1]) - (rho_DD + WL2_GPS * amb2);
            H.block<1, 3>(3 * m_DD + k, 0) << -dl, -dm, -dn;
            H(3 * m_DD + k, col2) = WL2_GPS;
        }
    }

    // 填充 BDS 观测方程
    if (n_DD > 0 && ref_idx[1] >= 0) {
        int off = 4 * m_DD;
        const SDSATOBS& sd_ref = Raw->SdObs.SdSatObs[ref_idx[1]];
        for (int k_C = 0; k_C < n_DD; ++k_C) {
            int k = m_DD + k_C;
            int i = Ekf->nPos[k];
            const SDSATOBS& sd = Raw->SdObs.SdSatObs[i];
            int nRov = sd.nRov;
            const double* spos = Raw->RovEpk.SatPVT[nRov].SatPos;
            double dx = spos[0] - X_R_est[0], dy = spos[1] - X_R_est[1], dz = spos[2] - X_R_est[2];
            double d_R_j = sqrt(dx * dx + dy * dy + dz * dz);
            double dl = (dx / d_R_j) - l_ref[1];
            double dm = (dy / d_R_j) - m_ref[1];
            double dn = (dz / d_R_j) - n_ref[1];
            double rho_DD = (d_R_j - d_R_ref[1]) - (d_B[i] - d_B[ref_idx[1]]);

            int col1 = 3 + 2 * k;
            int col2 = 3 + 2 * k + 1;
            double amb1 = X_act(col1);
            double amb2 = X_act(col2);

            // B1 伪距与载波
            y(off + 0 * n_DD + k_C) = (sd.dP[0] - sd_ref.dP[0]) - rho_DD;
            H.block<1, 3>(off + 0 * n_DD + k_C, 0) << -dl, -dm, -dn;

            y(off + 1 * n_DD + k_C) = (sd.dL[0] - sd_ref.dL[0]) - (rho_DD + WL1_BDS * amb1);
            H.block<1, 3>(off + 1 * n_DD + k_C, 0) << -dl, -dm, -dn;
            H(off + 1 * n_DD + k_C, col1) = WL1_BDS;

            // B3 伪距与载波
            y(off + 2 * n_DD + k_C) = (sd.dP[1] - sd_ref.dP[1]) - rho_DD;
            H.block<1, 3>(off + 2 * n_DD + k_C, 0) << -dl, -dm, -dn;

            y(off + 3 * n_DD + k_C) = (sd.dL[1] - sd_ref.dL[1]) - (rho_DD + WL3_BDS * amb2);
            H.block<1, 3>(off + 3 * n_DD + k_C, 0) << -dl, -dm, -dn;
            H(off + 3 * n_DD + k_C, col2) = WL3_BDS;
        }
    }

    // 3. 构建量测噪声协方差矩阵 R
    MatrixXd R = MatrixXd::Zero(l, l);
    if (m_DD > 0) {
        R.block(0 * m_DD, 0 * m_DD, m_DD, m_DD) = BuildRSubMatrix(m_DD, sigP);
        R.block(1 * m_DD, 1 * m_DD, m_DD, m_DD) = BuildRSubMatrix(m_DD, sigL);
        R.block(2 * m_DD, 2 * m_DD, m_DD, m_DD) = BuildRSubMatrix(m_DD, sigP);
        R.block(3 * m_DD, 3 * m_DD, m_DD, m_DD) = BuildRSubMatrix(m_DD, sigL);
    }
    if (n_DD > 0) {
        int off = 4 * m_DD;
        R.block(off + 0 * n_DD, off + 0 * n_DD, n_DD, n_DD) = BuildRSubMatrix(n_DD, sigP);
        R.block(off + 1 * n_DD, off + 1 * n_DD, n_DD, n_DD) = BuildRSubMatrix(n_DD, sigL);
        R.block(off + 2 * n_DD, off + 2 * n_DD, n_DD, n_DD) = BuildRSubMatrix(n_DD, sigP);
        R.block(off + 3 * n_DD, off + 3 * n_DD, n_DD, n_DD) = BuildRSubMatrix(n_DD, sigL);
    }

    // 4. 计算 Kalman 增益与状态更新
    MatrixXd Ht = H.transpose();
    MatrixXd S_mat = H * P_act * Ht + R;
    MatrixXd S_inv = S_mat.inverse();
    MatrixXd K = P_act * Ht * S_inv;

    VectorXd delta_X = K * y;
    X_act += delta_X;

    // 5. Joseph 稳定形式更新协方差矩阵: P = (I - K*H)*P*(I - K*H)^T + K*R*K^T
    MatrixXd I_KH = MatrixXd::Identity(n_act, n_act) - K * H;
    P_act = I_KH * P_act * I_KH.transpose() + K * R * K.transpose();

    // 6. 写回结构体 Ekf
    for (int r = 0; r < n_act; ++r) {
        Ekf->X[r] = X_act(r);
        for (int c = 0; c < n_act; ++c) {
            Ekf->P[r * TOTAL_STATE_DIM + c] = P_act(r, c);
        }
    }
    Ekf->Time = Raw->RovEpk.Time;
    Ekf->CurDDObs = Raw->DDObs;
    Ekf->IsInit = true;

    // 7. 更新流动站基线向量与精度指标
    Vector3d X_B(Raw->BasEpk.Pos[0], Raw->BasEpk.Pos[1], Raw->BasEpk.Pos[2]);
    Rov->Position[0] = X_act(0);
    Rov->Position[1] = X_act(1);
    Rov->Position[2] = X_act(2);
    Rov->IsSuccess = true;
    Rov->Time = Raw->RovEpk.Time;

    Raw->DDObs.dPos[0] = X_act(0) - X_B[0];
    Raw->DDObs.dPos[1] = X_act(1) - X_B[1];
    Raw->DDObs.dPos[2] = X_act(2) - X_B[2];

    // 验后残差平方和与验后单位权中误差 sigma0 [m]
    VectorXd v = y - H * delta_X;
    double vTRinvV = (v.transpose() * R.inverse() * v)(0);
    double sigma0 = (l > n_act) ? sqrt(std::max(0.0, vTRinvV / (l - n_act))) : 0.0;
    Rov->SigmaPos = sigma0;

    // 8. 提取浮点参数并调用 LAMBDA 执行模糊度搜索与固定基线恢复
    Vector3d b_hat = X_act.head<3>() - X_B;
    VectorXd a_hat = X_act.tail(2 * S);
    Matrix3d Qbb = P_act.block<3, 3>(0, 0);
    MatrixXd Qaa = P_act.block(3, 3, 2 * S, 2 * S);
    MatrixXd Qba = P_act.block(0, 3, 3, 2 * S);

    FixAmbiguity(Raw, Base, Rov, b_hat, a_hat, Qbb, Qaa, Qba, ratioThres);

    return true;
}

/* =============================================================================
 * 顶层调度入口: 单历元卡尔曼滤波函数 (RTKFilter)
 * 说明：单次调用严格推进单个历元，状态持久化保存在传入的 Ekf 结构体指针中。
 * ============================================================================= */
bool RTKFilter(RAWDAT* Raw, PPRESULT* Base, PPRESULT* Rov, RTKEKF* Ekf, const ROVERCFGINFO& cfg)
{
    if (!Raw || !Base || !Rov || !Ekf) return false;

    // 1. 基准站已知坐标维护与卫星轨道/高度角计算
    static Vector3d s_fixed_XB = Vector3d::Zero();
    static int s_base_src = 0;

    if (Raw->BasEpk.Pos[0] != 0.0 && Raw->BasEpk.Pos[1] != 0.0 && Raw->BasEpk.Pos[2] != 0.0) {
        Base->Position[0] = Raw->BasEpk.Pos[0];
        Base->Position[1] = Raw->BasEpk.Pos[1];
        Base->Position[2] = Raw->BasEpk.Pos[2];
        Base->IsSuccess = false;
        s_fixed_XB << Base->Position[0], Base->Position[1], Base->Position[2];
        s_base_src = 1;
    } else if (s_fixed_XB.norm() > 1000.0) {
        Base->Position[0] = s_fixed_XB[0];
        Base->Position[1] = s_fixed_XB[1];
        Base->Position[2] = s_fixed_XB[2];
        Base->IsSuccess = (s_base_src == 2);
    } else {
        DetectOutlier(&Raw->BasEpk);
        if (!SPP(&Raw->BasEpk, Raw, Base)) return false;
        s_fixed_XB << Base->Position[0], Base->Position[1], Base->Position[2];
        s_base_src = 2;
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

    // 2. 流动站执行 SPP 单点定位获取初值
    if (Raw->RovEpk.Pos[0] == 0.0 || (Vector3d(Raw->RovEpk.Pos) - X_B).norm() > 20000.0) {
        Raw->RovEpk.Pos[0] = Base->Position[0];
        Raw->RovEpk.Pos[1] = Base->Position[1];
        Raw->RovEpk.Pos[2] = Base->Position[2];
    }
    DetectOutlier(&Raw->RovEpk);
    SPP(&Raw->RovEpk, Raw, Rov);

    Vector3d X_R_star = X_B;
    if (Rov->IsSuccess && (Vector3d(Rov->Position) - X_B).norm() < 50000.0) {//是否真的需要有待考虑
        X_R_star << Rov->Position[0], Rov->Position[1], Rov->Position[2];
    }

    ComputeGPSSatOrbitAtSignalTrans(&Raw->RovEpk, Raw->GpsEph, Raw->BdsEph, X_R_star.data(), Raw->RovEpk.SatPVT);
    for (int i = 0; i < Raw->RovEpk.SatNum; ++i) {
        if (Raw->RovEpk.SatPVT[i].Valid) {
            double el = 0, az = 0;
            CompSatElAz(X_R_star.data(), Raw->RovEpk.SatPVT[i].SatPos, &el, &az);
            Raw->RovEpk.SatPVT[i].Elevation = el;
            Raw->RovEpk.SatPVT[i].Azimuth = az;
        }
    }

    // 3. 首次进入时执行滤波器状态初始化 (InitEKF)
    if (!Ekf->IsInit) {
        InitEKF(Ekf, X_R_star);
    }

    // 4. 选取参考星并统计双差可用卫星数
    DetRefSat(&Raw->BasEpk, &Raw->RovEpk, &Raw->SdObs, &Raw->DDObs);
    int m_DD = Raw->DDObs.DDSatNum[0];
    int n_DD = Raw->DDObs.DDSatNum[1];
    int S    = m_DD + n_DD;
    int l    = 4 * S;

    if (S < 3 || l <= 3) return false;

    int ref_idx[2] = { Raw->DDObs.RefPos[0], Raw->DDObs.RefPos[1] };
    int ref_prn[2] = { Raw->DDObs.RefPrn[0], Raw->DDObs.RefPrn[1] };

    // 5. 计算基准站几何距离
    std::vector<double> d_B(Raw->SdObs.SatNum, 0.0);
    for (int i = 0; i < Raw->SdObs.SatNum; ++i) {
        int nBas = Raw->SdObs.SdSatObs[i].nBas;
        const double* spos = Raw->BasEpk.SatPVT[nBas].SatPos;
        double dx = spos[0] - X_B[0], dy = spos[1] - X_B[1], dz = spos[2] - X_B[2];
        d_B[i] = sqrt(dx * dx + dy * dy + dz * dz);
    }

    // 6. 将当前历元双差星直接填入 Ekf->nPos 数组
    double elev_mask = (cfg.ElevThreshold > 0.0) ? cfg.ElevThreshold * Rad : 10.0 * Rad;
    int k_G = 0;
    if (ref_idx[0] >= 0 && m_DD > 0) {
        for (int i = 0; i < Raw->SdObs.SatNum && k_G < m_DD; ++i) {
            const SDSATOBS& sd = Raw->SdObs.SdSatObs[i];
            if (sd.System != GPS || sd.Valid == 0) continue;
            if (i == ref_idx[0]) continue;
            if (!Raw->RovEpk.SatPVT[sd.nRov].Valid || !Raw->BasEpk.SatPVT[sd.nBas].Valid) continue;
            if (Vector3d(Raw->BasEpk.SatPVT[sd.nBas].SatPos).norm() < 1e6 ||
                Vector3d(Raw->RovEpk.SatPVT[sd.nRov].SatPos).norm() < 1e6) continue;
            if (Raw->RovEpk.SatPVT[sd.nRov].Elevation < elev_mask) continue;

            Ekf->nPos[k_G] = i; // 直接记录单差数组下标！
            k_G++;
        }
    }

    int k_C = 0;
    if (ref_idx[1] >= 0 && n_DD > 0) {
        for (int i = 0; i < Raw->SdObs.SatNum && k_C < n_DD; ++i) {
            const SDSATOBS& sd = Raw->SdObs.SdSatObs[i];
            if (sd.System != BDS || sd.Valid == 0) continue;
            if (i == ref_idx[1]) continue;
            if (!Raw->RovEpk.SatPVT[sd.nRov].Valid || !Raw->BasEpk.SatPVT[sd.nBas].Valid) continue;
            if (Vector3d(Raw->BasEpk.SatPVT[sd.nBas].SatPos).norm() < 1e6 ||
                Vector3d(Raw->RovEpk.SatPVT[sd.nRov].SatPos).norm() < 1e6) continue;
            if (Raw->RovEpk.SatPVT[sd.nRov].Elevation < elev_mask) continue;

            Ekf->nPos[m_DD + k_C] = i; // 直接记录单差数组下标！
            k_C++;
        }
    }

    m_DD = k_G;
    n_DD = k_C;
    S = m_DD + n_DD;
    l = 4 * S;
    if (S < 3 || l <= 3) return false;
    Ekf->nSats = S;

    // 7. 检查参考星是否切换
    bool ref_changed[2] = { false, false };
    if (Ekf->CurDDObs.RefPrn[0] != 0 && Ekf->CurDDObs.RefPrn[0] != ref_prn[0]) ref_changed[0] = true;
    if (Ekf->CurDDObs.RefPrn[1] != 0 && Ekf->CurDDObs.RefPrn[1] != ref_prn[1]) ref_changed[1] = true;

    // 8. 依次执行滤波子模块：时间更新与量测更新
    double q_pos = 1.0e-6; // 静态短基线模式过程噪声 [m^2]
    TimeUpdateEKF(Ekf, Raw, ref_idx, m_DD, n_DD, ref_changed, q_pos);

    double sigP = (cfg.CodeNoise > 0.0) ? cfg.CodeNoise : 0.3;
    double sigL = (cfg.CPNoise > 0.0)   ? cfg.CPNoise   : 0.003;
    return MeasUpdateEKF(Ekf, Raw, Base, Rov, d_B, ref_idx, m_DD, n_DD, sigP, sigL, cfg.RatioThres);
}
