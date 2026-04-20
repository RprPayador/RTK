#include "RTK_Structs.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;


// 信号发射时刻卫星位置速度计算 
void ComputeGPSSatOrbitAtSignalTrans(const EPOCHOBS* Epk, GPSEPHREC* GpsEph, GPSEPHREC* BDSEph, double RcvPos[3], SATMIDRES* MidRes) {
    for (int i = 0; i < Epk->SatNum; i++) {
        const SATOBS& obs = Epk->SatObs[i];
        SATMIDRES& pvt = MidRes[i]; // 将计算结果写入外部传入的缓冲区
        pvt.Valid = false;

        // 直接获取已在 detect.cpp 中计算好的无电离层组合伪距 (PIF)
        double rho = Epk->ComObs[i].PIF;
        if (rho <= 0) continue; 

        // 1. 初始化卫星钟差
        double dts = 0.0;
        
        // 2. 迭代计算信号发射时刻 (针对该历元观测时刻进行倒推)
        GPSTIME t_tr;
        for (int iter = 0; iter < 2; iter++) {
            t_tr = Epk->Time;
            t_tr.SecOfWeek -= (rho / C_Light + dts);
            
            // 处理跨周
            if (t_tr.SecOfWeek < 0) { t_tr.SecOfWeek += 604800.0; t_tr.Week -= 1; }
            if (t_tr.SecOfWeek >= 604800.0) { t_tr.SecOfWeek -= 604800.0; t_tr.Week += 1; }
            
            // 计算卫星钟差
            if (!CompSatClkOff(obs.Prn, obs.System, &t_tr, GpsEph, BDSEph, &pvt)) break;

            // 计算卫星轨迹（包含相对论效应改正）
            bool stat = false;
            if (obs.System == GPS) stat = CompGPSSatPVT(obs.Prn, &t_tr, &GpsEph[obs.Prn - 1], &pvt);
            else stat = CompBDSSatPVT(obs.Prn, &t_tr, &BDSEph[obs.Prn - 1], &pvt);
            
            if (!stat) break;
            dts = pvt.SatClkOft;
        }
        
        if (!pvt.Valid) continue;

        // 3. 地球自转改正
        double dt1 = rho / C_Light;
        if (RcvPos[0] != 0 || RcvPos[1] != 0 || RcvPos[2] != 0) {
            double dx = pvt.SatPos[0] - RcvPos[0];
            double dy = pvt.SatPos[1] - RcvPos[1];
            double dz = pvt.SatPos[2] - RcvPos[2];
            dt1 = sqrt(dx*dx + dy*dy + dz*dz) / C_Light;
        }

        double omega_e = (obs.System == BDS) ? Omega_BDS : Omega_WGS;
        double theta = omega_e * dt1;

        double cos_t = cos(theta);
        double sin_t = sin(theta);

        double Xs = pvt.SatPos[0];
        double Ys = pvt.SatPos[1];
        double Vxs = pvt.SatVel[0];
        double Vys = pvt.SatVel[1];
        
        pvt.SatPos[0] = Xs * cos_t + Ys * sin_t;
        pvt.SatPos[1] = -Xs * sin_t + Ys * cos_t;
        pvt.SatVel[0] = Vxs * cos_t + Vys * sin_t;
        pvt.SatVel[1] = -Vxs * sin_t + Vys * cos_t;
        
        pvt.Valid = true;
    }
}


// 单点定位解算主函数 (严格对齐 RTK_Structs.h 第 397 行定义)
bool SPP(EPOCHOBS* Epoch, RAWDAT* Raw, PPRESULT* Result) {
    
    // 1. 设置初略位置与钟差
    Vector3d X_R(0, 0, 0);
    double dt_G = 0.0; 
    double dt_C = 0.0; 

    // 初值继承
    if (Epoch->Pos[0] != 0 && Epoch->Pos[1] != 0 && Epoch->Pos[2] != 0) {
        X_R << Epoch->Pos[0], Epoch->Pos[1], Epoch->Pos[2];
    }

    // 2. 迭代求解过程
    bool is_converged = false;
    double PDOP = 999.9;
    
    for (int iter = 0; iter < 25; ++iter) {
        
        // A. 轨道精化计算
        ComputeGPSSatOrbitAtSignalTrans(Epoch, Raw->GpsEph, Raw->BdsEph, X_R.data(), Epoch->SatPVT);

        // B. 卫星可用性统计
        std::vector<int> valid_sats;
        int gps_count = 0, bds_count = 0;
        
        for (int i = 0; i < Epoch->SatNum; i++) {
            SATMIDRES& pvt = Epoch->SatPVT[i];
            if (!pvt.Valid) continue;
            
            // 计算高度角（用于对流层改正和信息输出）
            double el = 0, az = 0;
            if (X_R.norm() > 1000.0) {
                double sat_pos[3] = {pvt.SatPos[0], pvt.SatPos[1], pvt.SatPos[2]};
                CompSatElAz(X_R.data(), sat_pos, &el, &az);
                pvt.Elevation = el; pvt.Azimuth = az;
            }
            valid_sats.push_back(i);
            if (Epoch->SatObs[i].System == GPS) gps_count++;
            if (Epoch->SatObs[i].System == BDS) bds_count++;
        }
        
        int total_sats = valid_sats.size();
        int state_size = 3 + (gps_count > 0 ? 1 : 0) + (bds_count > 0 ? 1 : 0);//判断最小二乘矩阵列数是4还是5
        
        if (total_sats < state_size) return false;//说明卫星数小于4
        
        // C. 建立线性化方程
        int idx_G = (gps_count > 0) ? 3 : -1;
        int idx_C = (bds_count > 0) ? (3 + (gps_count > 0 ? 1 : 0)) : -1;

        MatrixXd B = MatrixXd::Zero(total_sats, state_size);
        VectorXd w(total_sats);
        
        for (int k = 0; k < total_sats; k++) {
            int idx = valid_sats[k];
            SATOBS& obs = Epoch->SatObs[idx];
            SATMIDRES& pvt = Epoch->SatPVT[idx];
            double rho_IF = Epoch->ComObs[idx].PIF;
            
            double dx = pvt.SatPos[0] - X_R[0];
            double dy = pvt.SatPos[1] - X_R[1];
            double dz = pvt.SatPos[2] - X_R[2];
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            if (r < 1e-3) r = 1e-3;
            
            B(k, 0) = -dx / r;
            B(k, 1) = -dy / r;
            B(k, 2) = -dz / r;
            
            if (obs.System == GPS && idx_G != -1) B(k, idx_G) = 1.0;
            if (obs.System == BDS && idx_C != -1) B(k, idx_C) = 1.0;

            // 模型改正 (对流层 + BDS TGD)
            double T_trop = 0.0;
            if (X_R.norm() > 1000.0) {
                double blh[3];
                XYZToBLH(X_R.data(), blh, R_WGS84, F_WGS84);
                T_trop = Hopfield(blh[2], pvt.Elevation);
                pvt.TropCorr = T_trop; // 保存以便输出
            }
            
            double Tgd_corr = 0.0;
            if (obs.System == BDS) {
                // BDS B1/B3 IF TGD: (f1^2 / (f1^2 - f3^2)) * TGD1
                double f1 = 1561.098e6, f3 = 1268.52e6;
                Tgd_corr = (f1*f1) / (f1*f1 - f3*f3) * pvt.Tgd1 * C_Light;
            }
            double current_dt_r = (obs.System == GPS) ? dt_G : dt_C;
            w(k) = rho_IF - (r + current_dt_r - C_Light * pvt.SatClkOft + T_trop + Tgd_corr);
        }
        
        // D. 最小二乘求解 (等权)
        MatrixXd BTB = B.transpose() * B;
        VectorXd dx_hat = BTB.ldlt().solve(B.transpose() * w);
        
        X_R[0] += dx_hat(0);
        X_R[1] += dx_hat(1);
        X_R[2] += dx_hat(2);
        
        if (idx_G != -1) dt_G += dx_hat(idx_G);
        if (idx_C != -1) dt_C += dx_hat(idx_C);
        
        // E. 收敛检查与精度评定
        if (dx_hat.head<3>().norm() < 0.001) {
            MatrixXd Q = BTB.inverse();
            PDOP = sqrt(Q(0, 0) + Q(1, 1) + Q(2, 2));
            
            VectorXd v = B * dx_hat - w;
            Result->SigmaPos = sqrt(v.dot(v) / (total_sats - state_size));
            Result->Time = Epoch->Time;
            is_converged = true;
            break;
        }
    }
    
    // 3. 输出平差结果
    if (is_converged) {
        Result->Position[0] = X_R[0]; Result->Position[1] = X_R[1]; Result->Position[2] = X_R[2];
        Result->RcvClkOft[0] = dt_G / C_Light; Result->RcvClkOft[1] = dt_C / C_Light;
        Result->PDOP = PDOP;
        
        int gps_num = 0, bds_num = 0;
        for (int i = 0; i < Epoch->SatNum; i++) {
            if (!Epoch->SatPVT[i].Valid) continue;
            if (Epoch->SatObs[i].System == GPS) gps_num++;
            if (Epoch->SatObs[i].System == BDS) bds_num++;
        }
        Result->GPSSatNum = gps_num;
        Result->BDSSatNum = bds_num;
        Result->AllSatNum = gps_num + bds_num;
        Result->IsSuccess = true;
        
        // 坐标状态传递
        Epoch->Pos[0] = X_R[0]; Epoch->Pos[1] = X_R[1]; Epoch->Pos[2] = X_R[2];
        return true;
    }
    
    Result->IsSuccess = false;
    return false;
}
