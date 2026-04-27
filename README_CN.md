# GNSS 单点定位与测速 (SPP/SPV)

[![Language](https://img.shields.io/badge/language-C%2B%2B-blue.svg)](https://isocpp.org/)
[![Platform](https://img.shields.io/badge/platform-Windows-lightgrey.svg)](https://www.microsoft.com/windows)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

> 基于 C++ 的 GNSS 单点定位与测速程序，支持 GPS 和北斗双系统。

## 项目简介

本项目实现了一套完整的 GNSS 数据处理流程，支持单点定位和单点测速功能。可通过 TCP 套接字实时接收 NovAtel OEM7 接收机数据，或离线处理二进制数据文件。

### 主要特性

- **双系统支持**：GPS (L1/L2) 和北斗 (B1/B3) 卫星导航系统
- **SPP & SPV**：单点定位与单点测速
- **实时与离线**：支持 TCP 网络数据流和本地文件两种模式
- **粗差探测**：MW (墨尔本-乌比纳) 组合与 GF (几何无关) 组合周跳探测
- **误差改正**：Hopfield 对流层模型、TGD 改正、相对论效应改正
- **热启动**：利用上一历元解算结果作为当前历元初值

## 项目架构

```
RTK_Raw/
├── RTK.cpp              # 主程序入口与数据处理循环
├── RTK_Structs.h        # 核心数据结构体与常数定义
├── DecodeNovOem7.cpp    # NovAtel OEM7 二进制协议解码器
├── SPP.cpp              # SPP/SPV 算法实现（伪距定位与多普勒测速）
├── SatPVT.cpp           # 卫星轨道与钟差计算
├── CoordinateConvert.cpp # 坐标转换（地心地固/大地/站心坐标系）
├── detect.cpp           # 粗差探测与周跳检测
├── TimeConvert.cpp      # 时间系统转换
├── sockets.cpp          # TCP 套接字封装（实时数据接收）
└── data/                # 示例数据文件
```

## 编译环境

- **编译器**：MinGW-w64 GCC 或 MSVC，支持 C++11 标准
- **依赖库**：
  - [Eigen3](https://eigen.tuxfamily.org/)（线性代数库）
  - Windows Socket API (`ws2_32.lib`)
- **操作系统**：Windows（因使用 WinSock 网络库）

## 编译说明

### 使用 VS Code

1. 在 VS Code 中打开本项目
2. 按 `Ctrl+Shift+B` 执行构建任务
3. 构建配置位于 `.vscode/tasks.json`

### 命令行编译

```bash
cd RTK_Raw
g++ -fdiagnostics-color=always -g \
    RTK.cpp DecodeNovOem7.cpp CoordinateConvert.cpp \
    TimeConvert.cpp SatPVT.cpp SPP.cpp detect.cpp \
    -I "D:\GNSS Algorithm\RTK\RTK\eigen-5.0.0\eigen-5.0.0" \
    -o RTK.exe \
    -lws2_32
```

## 使用说明

### 配置模式

在 `RTK_Structs.h` 中设置运行模式：

```cpp
#define FILEMODE 1   // 1 = 文件模式（离线处理），0 = 实时 TCP 网络模式
```

### 文件模式

将 NovAtel OEM7 二进制数据文件放置于：
```
D:\GNSS Algorithm\RTK\RTK\oem719-202603111200.bin
```

### 实时网络模式

在 `RTK.cpp` 中配置 TCP 服务器地址：
```cpp
if(OpenSocket(NetGps, "47.114.134.129", 7190) == false)
```

### 运行程序

```bash
./RTK.exe
```

解算结果将输出至 `result.txt`。

## 输出格式

### 单颗卫星信息（每个历元）
```
G01 X= -12694065.000 Y=   8930477.000 Z=  21578501.000 Clk=-7.123456e-04 Vx= -945.6000 Vy=  2905.7000 Vz=  1373.1000 Clkd= 3.45678e-11 PIF=20200356.1000 Trop=  2.100 E= 45.123deg
```

### 定位结果（每个历元）
```
SPP: 2294 273600.000 X:-1269406.5000 Y:893047.7000 Z:2157850.1000 B:   30.12345678 L:  120.98765432 H:  123.456 Vx:   0.0000 Vy:   0.0000 Vz:   0.0000 GPS Clk:     12.345 BDS Clk:      0.000 PDOP:   1.234 Sigma:   2.345 GPSSats:  8 BDSSats:  6 Sats: 14
```

## 算法说明

### 1. 数据解码
- **协议**：NovAtel OEM7 二进制数据格式
- **消息类型**：RANGEB（原始观测值）、GPSEPHEM/BDSEPHEM（广播星历）、PSRPOS（参考位置）
- **观测数据**：双频伪距、载波相位、多普勒观测值

### 2. 粗差探测
- **MW 组合**：宽巷模糊度检测（阈值：3.0 周）
- **GF 组合**：电离层延迟变化检测（阈值：0.05 米）
- **双频 PIF**：无电离层伪距组合用于定位解算

### 3. SPP 定位算法
- **估计方法**：最小二乘平差，迭代线性化
- **状态向量**：[X, Y, Z, dt_GPS, dt_BDS]（三维坐标 + GPS/北斗双系统接收机钟差）
- **收敛条件**：最大迭代 25 次，状态变化量小于 0.1 毫米
- **质量指标**：PDOP（位置精度因子）、Sigma（验后单位权中误差）

### 4. SPV 测速算法
- **观测数据**：多普勒观测值转换为距离变化率
- **状态向量**：[Vx, Vy, Vz, dtr_dot]（三维速度 + 接收机钟漂）
- **模型**：利用广播星历计算卫星速度，结合站星几何关系

### 5. 误差改正
- **对流层延迟**：Hopfield 模型，标准气象参数
- **TGD 改正**：北斗 B1/B3 无电离层组合的 TGD 改正
- **相对论效应**：卫星钟差的 Schwarzschild 项改正
- **地球自转**：信号传播期间的 Sagnac 效应改正

## 坐标系统

| 坐标系统 | 说明 | 参考标准 |
|---------|------|---------|
| WGS-84 | GPS 默认坐标系统 | GPS ICD |
| CGCS2000 | 北斗坐标系统 | 北斗 ICD |
| ENU | 站心东北天坐标系 | 测站中心 |

## 物理常数

| 参数 | 数值 | 来源 |
|-----|------|------|
| 光速 | 299,792,458 m/s | IAU 1976 |
| WGS-84 长半轴 a | 6,378,137.0 m | WGS-84 |
| WGS-84 扁率 f | 1/298.257223563 | WGS-84 |
| 地心引力常数 GM | 398600.5×10⁹ m³/s² | WGS-84 |
| 地球自转速率 Ω (WGS) | 7.2921151467×10⁻⁵ rad/s | WGS-84 |
| 地球自转速率 Ω (BDS) | 7.2921150×10⁻⁵ rad/s | CGCS2000 |

## 项目结构

```
RTK_Raw/
├── .vscode/                 # VS Code 配置
│   ├── c_cpp_properties.json
│   ├── launch.json
│   ├── settings.json
│   └── tasks.json
├── OEM719原始数据/          # 原始数据样本
├── data/                    # 处理后数据输出
├── .gitignore              # Git 忽略规则
├── CoordinateConvert.cpp   # 坐标转换函数
├── DecodeNovOem7.cpp       # OEM7 协议解码
├── detect.cpp              # 周跳探测算法
├── RTK.cpp                 # 主程序
├── RTK_Structs.h           # 数据结构体与常数
├── SatPVT.cpp              # 卫星位置速度计算
├── SPP.cpp                 # SPP/SPV 核心算法
├── sockets.cpp             # 网络套接字封装
├── test_CompSatElAz.cpp    # 高度角方位角计算测试
└── TimeConvert.cpp         # 时间系统转换
```

## 参与贡献

欢迎提交 Issue 或 Pull Request！

## 致谢

- 本项目由武汉大学测绘学院开发
- 广播星历算法参考 GPS ICD 200 与北斗 ICD 规范

## 参考资料

1. IS-GPS-200N. *Navstar GPS Space Segment/Navigation User Segment Interface*
2. 北斗卫星导航系统空间信号接口控制文件 - 公开服务信号B1I (Version 3.0)
3. NovAtel OEM7 Firmware Reference Manual
4. Hofmann-Wellenhof et al., *GNSS – Global Navigation Satellite Systems*, Springer

## 许可证

MIT 许可证 - 详见 [LICENSE](LICENSE) 文件。

---

**作者**：RprPayador  
**单位**：  
**邮箱**: 17707234049lijy@gmail.com  
**版本**：1.0
