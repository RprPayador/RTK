#include "RTK_Structs.h"
#include <string.h>
#include <math.h>


#define POLYCRC32 0xEDB88320u

unsigned int crc32(const unsigned char* buff, int len) {
	int i, j;
	unsigned int crc = 0;
	for (i = 0; i < len; i++) {
		crc ^= buff[i];
		for (j = 0; j < 8; j++) {
			if (crc & 1) crc = (crc >> 1) ^ POLYCRC32;
			else crc >>= 1;
		}
	}
	return crc;
}

int DecodeNovOem7Dat(unsigned char Buff[], int& Len, EPOCHOBS* obs, GPSEPHREC geph[], GPSEPHREC beph[]) {
	int i = 0;
	int ret = -1;

	while (i + 3 < Len) {
		if (Buff[i] == 0xAA && Buff[i + 1] == 0x44 && Buff[i + 2] == 0x12) {//查找同步头
			if (i + 28 > Len)
				break;

			unsigned short msgID = *(unsigned short*)(Buff + i + 4);//读取消息ID
			unsigned short msgLen = *(unsigned short*)(Buff + i + 8);//读取消息长度

			if (i + 28 + msgLen + 4 > Len) // 总长 = 头28 + 体msgLen + 校验4
				break;


			unsigned long crcInCalculate = crc32(Buff + i, 28 + msgLen);//计算CRC
			unsigned long crcInMsg = *(unsigned long*)(Buff + i + 28 + msgLen);//读取CRC

			if (crcInCalculate != crcInMsg) { // CRC校验失败，跳过同步头继续搜索
				i += 3;
				continue;
			}

			switch (msgID) {
			case 1696: // BDSEPHEM: 北斗星历
				ret = decode_bdsephem(Buff + i, beph); 
				break;
			case 7:    // GPSEPHEM: GPS星历
				ret = decode_gpsephem(Buff + i, geph); 
				break;
			case 43:   // RANGEB: 原始观测值
				ret = decode_rangeb_oem7(Buff + i, obs); 
				break;
			case 47:   // PSRPOS: 伪距定位解
				PPRESULT pos;
				ret = decode_psrpos(Buff + i, &pos); 
				if (ret == 2) {
					obs->Pos[0] = pos.Position[0];
					obs->Pos[1] = pos.Position[1];
					obs->Pos[2] = pos.Position[2];
				}
				break;
			}

			i = i + 28 + msgLen + 4; // 处理完毕，移动到下一条消息起始位置

			if (ret == 1)//如果是观测值消息，直接返回
				break;
		}

		else i++;//查找下一个同步头

	}


	int remain = Len - i;
	memmove(Buff, Buff + i, remain);//剩余消息前移
	Len = remain;
	return ret;
}

int decode_rangeb_oem7(unsigned char* buff, EPOCHOBS* obs) {
	unsigned char* body = buff + 28;

	// 1. 从 Header 提取 GPS 时间
	unsigned short week = *(unsigned short*)(buff + 14); 
	unsigned long  ms = *(unsigned long*)(buff + 16);    

	int i = 0;
	unsigned long Nobs = *(unsigned long*)body;//观测值总数
	i = 4;

	int satNum = 0;
	obs->SatNum = 0;
	for (int k = 0; k < MAXCHANNUM; k++) {
		obs->SatObs[k] = SATOBS();    // 重置观测值
		obs->SatPVT[k] = SATMIDRES(); // 重置卫星位置等计算结果
	}

	obs->Time.Week = (short)week;
	obs->Time.SecOfWeek = ms / 1000.0;

	// 2. 遍历每条观测记录 (每条记录长 44 字节)
	for (int j = 0; j < (int)Nobs; j++) {
		// 简单的缓冲区安全检查
		if (i + 44 > MAXRAWLEN) break; 

		unsigned short prn = *(unsigned short*)(body + i);
		double psr = *(double*)(body + 4 + i);
		float psr_sigma = *(float*)(body + 12 + i);
		double adr = *(double*)(body + 16 + i);
		float adr_sigma = *(float*)(body + 24 + i);
		float dopp = *(float*)(body + 28 + i);
		float CNo = *(float*)(body + 32 + i);
		float locktime = *(float*)(body + 36 + i);
		unsigned long status = *(unsigned long*)(body + 40 + i);

		int SystemID = (status >> 16) & 0x07;
		int SignalTypeID = (status >> 21) & 0x1F;
		int half = status & 0x01;


		GNSSSys sys = UNKS;
		int s = -1;

		if (SystemID == 0) { // GPS 系统
			sys = GPS;
			if (SignalTypeID == 0) s = 0; // L1 C/A
			else if (SignalTypeID == 9) s = 1; // L2 P(Y)
		}
		else if (SystemID == 4) { // BDS 北斗系统
			sys = BDS;
			if (SignalTypeID == 0 || SignalTypeID == 4)s = 0; // B1I (D1=0, D2=4)
			else if (SignalTypeID == 2 || SignalTypeID == 6)s = 1; // B3I (D1=2, D2=6)
		}

		// 如果不是我们关注的系统或频点，跳过
		if (sys == UNKS || s == -1) {
			i += 44;
			continue;
		}

		// 3. 查找或加入卫星列表
		int idx = -1;
		for (int k = 0; k < satNum; k++) {
			if (obs->SatObs[k].System == sys && obs->SatObs[k].Prn == prn) {
				idx = k;
				break;
			}
		}

		if (idx == -1 && satNum < MAXCHANNUM) {
			idx = satNum++;
			obs->SatObs[idx].Prn = prn;
			obs->SatObs[idx].System = sys;
			obs->SatObs[idx].Valid = true;
		}

		// 4. 填充观测值数据
		if (idx != -1) {
			if (s == 0) {
				obs->SatObs[idx].c1 = psr;
				obs->SatObs[idx].l1 = adr;
				obs->SatObs[idx].d1 = dopp;
				obs->SatObs[idx].cn0[0] = CNo;
				obs->SatObs[idx].LockTime[0] = locktime;
				obs->SatObs[idx].half[0] = (unsigned char)half;
			}
			else if (s == 1) {
				obs->SatObs[idx].p2 = psr;
				obs->SatObs[idx].l2 = adr;
				obs->SatObs[idx].d2 = dopp;
				obs->SatObs[idx].cn0[1] = CNo;
				obs->SatObs[idx].LockTime[1] = locktime;
				obs->SatObs[idx].half[1] = (unsigned char)half;
			}
		}
		i += 44;
	}
    obs->SatNum = (short)satNum;

	return 1;
}

int decode_gpsephem(unsigned char* buff, GPSEPHREC* eph) {
	unsigned char* body = buff + 28;
	unsigned long prn = *(unsigned long*)body;

	if (prn < 1 || prn > MAXGPSNUM) return -1;

	GPSEPHREC* p = &eph[prn - 1];
	p->PRN = (unsigned short)prn;
	p->Sys = GPS;
	p->SVHealth = short(*(unsigned long*)(body + 12));
	p->IODE = double(*(unsigned long*)(body + 16));
	p->TOC.Week = p->TOE.Week = short(*(unsigned long*)(body + 24));
	p->TOE.SecOfWeek = *(double*)(body + 32);
	p->SqrtA = sqrt(*(double*)(body + 40));
	p->DeltaN = *(double*)(body + 48);
	p->M0 = *(double*)(body + 56);
	p->e = *(double*)(body + 64);
	p->omega = *(double*)(body + 72);
	p->Cuc = *(double*)(body + 80);
	p->Cus = *(double*)(body + 88);
	p->Crc = *(double*)(body + 96);
	p->Crs = *(double*)(body + 104);
	p->Cic = *(double*)(body + 112);
	p->Cis = *(double*)(body + 120);
	p->i0 = *(double*)(body + 128);
	p->iDot = *(double*)(body + 136);
	p->OMEGA = *(double*)(body + 144);
	p->OMEGADot = *(double*)(body + 152);
	p->IODC = double(*(unsigned long*)(body + 160));
	p->TOC.SecOfWeek = *(double*)(body + 164);
	p->TGD1 = *(double*)(body + 172);
	p->ClkBias = *(double*)(body + 180);
	p->ClkDrift = *(double*)(body + 188);
	p->ClkDriftRate = *(double*)(body + 196);
	p->SVAccuracy = *(double*)(body + 216);

	return 0;

}





int decode_bdsephem(unsigned char* buff, GPSEPHREC* eph) {
	unsigned char* body = buff + 28;
	unsigned long prn = *(unsigned long*)body;

	if (prn < 1 || prn > MAXBDSNUM) return -1;

	GPSEPHREC* p = &eph[prn - 1];

	p->Sys = BDS;
	p->PRN = (unsigned short)prn;
	p->TOC.Week = p->TOE.Week = short(*(unsigned long*)(body + 4));
    p->SVAccuracy = *(double*)(body + 8);
	p->SVHealth = short(*(unsigned long*)(body + 16));
	p->TGD1 = *(double*)(body + 20);
	p->TGD2 = *(double*)(body + 28);
    p->IODC = double(*(unsigned long*)(body + 36));
	p->TOC.SecOfWeek = double(*(unsigned long*)(body + 40));
	p->ClkBias = *(double*)(body + 44);
	p->ClkDrift = *(double*)(body + 52);
	p->ClkDriftRate = *(double*)(body + 60);
    p->IODE = double(*(unsigned long*)(body + 68));
	p->TOE.SecOfWeek = double(*(unsigned long*)(body + 72));
	p->SqrtA = *(double*)(body + 76);
    p->e=*(double*)(body+84);
	p->omega = *(double*)(body + 92);
	p->DeltaN = *(double*)(body + 100);
	p->M0 = *(double*)(body + 108);
	p->OMEGA = *(double*)(body + 116);
	p->OMEGADot = *(double*)(body + 124);
	p->i0 = *(double*)(body + 132);
    p->iDot = *(double*)(body + 140);
	p->Cuc = *(double*)(body + 148);
	p->Cus = *(double*)(body + 156);
	p->Crc = *(double*)(body + 164);
	p->Crs = *(double*)(body + 172);
	p->Cic = *(double*)(body + 180);
	p->Cis = *(double*)(body + 188);

	return 0;
}

int decode_psrpos(unsigned char* buff, PPRESULT* pos) { 
	unsigned char* body = buff + 28;

	unsigned long sol_status = *(unsigned long*)(body + 0);
	double lat = *(double*)(body + 8);
	double lon = *(double*)(body + 16);
	double hgt = *(double*)(body + 24);

	double BLH[3], XYZ[3];
	BLH[0] = lat * Rad;
	BLH[1] = lon * Rad;
	BLH[2] = hgt;
	BLHToXYZ(BLH, XYZ, R_WGS84, F_WGS84);

	pos->Position[0] = XYZ[0];
	pos->Position[1] = XYZ[1];
	pos->Position[2] = XYZ[2];
	pos->IsSuccess = (sol_status == 0);

	return 2; 
}
