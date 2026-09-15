#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>          // fabs() 所需
#include "RTK_Structs.h"

/* 静态缓冲区，跨调用保持未处理的字节
 * DecodeNovOem7Dat 解不完的半帧数据留在 Buff 头部、Len 记录长度,
 * 下次调用从 Len 处继续读、继续解。若不加 static,函数返回后
 * 半帧字节丢失,字节流会错位,导致随机丢历元/CRC 校验失败 */

// 返回值:  1=同步成功  0=文件结束  -1=未同步(继续读)
int GetSynObs(FILE* FBas, FILE* FRov, SOCKET& BasSock, SOCKET& RovSock, RAWDAT* Raw, int IsFileData, double SynThres){
    static unsigned char RovBuff[MAXRAWLEN], BasBuff[MAXRAWLEN];
    static int RovLen = 0, BasLen = 0;
    int RovReadLen, BasReadLen;
    int val;

    // 读流动站，解码出一个完整历元
    bool RovDone = false;
    while (!RovDone) {
        if (IsFileData == 1) {
            RovReadLen = fread(RovBuff + RovLen, 1, MAXRAWLEN - RovLen, FRov);
            if (RovReadLen == 0) {//文件读完
                while(true){
                    val=DecodeNovOem7Dat(RovBuff,RovLen,&Raw->RovEpk,Raw->GpsEph,Raw->BdsEph);
                    if(val==-1) return 0;//数据解码完
                    if(val==1) {RovDone=true;break;}
                }
                break;
            }  
        } else {
            Sleep(980);
            RovReadLen = recv(RovSock, (char*)(RovBuff + RovLen), MAXRAWLEN - RovLen, 0);
            if (RovReadLen == 0) return 0;
        }
        RovLen += RovReadLen;

        while (true) {
            val = DecodeNovOem7Dat(RovBuff, RovLen, &Raw->RovEpk, Raw->GpsEph, Raw->BdsEph);
            if (val == -1) break;           // 帧不完整，回去继续读
            if (val == 1) { RovDone = true; break; }  // 解出观测历元
        }
    }

    // 比较流动站时刻与上次基站时刻
    double dt = GetDiffTime(&Raw->RovEpk.Time, &Raw->BasEpk.Time);

    // 流动站比基站早，说明基站还没跟上，暂时无法同步
    if (dt < -SynThres) return -1;
    if (std::fabs(dt)<=SynThres) return 1;

    // 如果基站时刻已经同步（或基站未初始化 dt 很大），进入读基站
    // 读基站，直到基站时刻追上流动站
    while (dt > SynThres) {
        if (IsFileData == 1) {
            BasReadLen = fread(BasBuff + BasLen, 1, MAXRAWLEN - BasLen, FBas);
            if (BasReadLen == 0) {//文件读完
                while(true){
                    val=DecodeNovOem7Dat(BasBuff,BasLen,&Raw->BasEpk,Raw->GpsEph,Raw->BdsEph);
                    if(val==-1) return 0;//数据解码完
                    if (val == 1) {                                 
                        dt = GetDiffTime(&Raw->RovEpk.Time,&Raw->BasEpk.Time);                                       
                        if (dt < -SynThres) return -1;              
                        if (std::fabs(dt) <= SynThres) return 1;
                    }
                }  
            } 
        }
        else {
            BasReadLen = recv(BasSock, (char*)(BasBuff + BasLen), MAXRAWLEN - BasLen, 0);
            if (BasReadLen == 0) {//文件读完
                while(true){
                    val=DecodeNovOem7Dat(BasBuff,BasLen,&Raw->BasEpk,Raw->GpsEph,Raw->BdsEph);
                    if(val==-1) return 0;//数据解码完
                    if (val == 1) {                                 
                        dt = GetDiffTime(&Raw->RovEpk.Time,&Raw->BasEpk.Time);                                       
                        if (dt < -SynThres) return -1;              
                        if (std::fabs(dt) <= SynThres) return 1;   // 追平即返回                                                  
                    } 
                }
            }
        }
        BasLen += BasReadLen;

        while (true) {
            val = DecodeNovOem7Dat(BasBuff, BasLen, &Raw->BasEpk, Raw->GpsEph, Raw->BdsEph);
            if (val == -1) break;
            if (val == 1) {
                dt = GetDiffTime(&Raw->RovEpk.Time, &Raw->BasEpk.Time);
                if (dt < -SynThres) return -1; 
                else if(std::fabs(dt)<=SynThres) return 1;
            }
        }
    }
    return 1;

}