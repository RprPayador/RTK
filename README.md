# GNSS Single Point Positioning (SPP)

[![Language](https://img.shields.io/badge/language-C%2B%2B-blue.svg)](https://isocpp.org/)
[![Platform](https://img.shields.io/badge/platform-Windows-lightgrey.svg)](https://www.microsoft.com/windows)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

> A C++ implementation of GNSS Single Point Positioning (SPP) and Single Point Velocity (SPV) supporting both GPS and BeiDou (BDS) systems.

## Overview

This project implements a complete GNSS data processing pipeline for single point positioning and velocity estimation. It supports real-time data streaming from NovAtel OEM7 receivers via TCP socket or offline processing of binary data files.

### Key Features

- **Dual-System Support**: GPS (L1/L2) and BeiDou (B1/B3) satellite systems
- **SPP & SPV**: Single Point Positioning and Single Point Velocity estimation
- **Real-time & Offline**: Support for TCP socket streaming and file-based processing
- **Outlier Detection**: MW (Melbourne-Wübbena) and GF (Geometry-Free) combination for cycle slip detection
- **Error Corrections**: Troposphere delay (Hopfield model), TGD corrections, and relativistic effects
- **Hot Start**: Position warm-start capability using previous epoch solutions

## Architecture

```
RTK_Raw/
├── RTK.cpp              # Main entry point and processing loop
├── RTK_Structs.h        # Core data structures and constants
├── DecodeNovOem7.cpp    # NovAtel OEM7 binary protocol decoder
├── SPP.cpp              # SPP/SPV algorithms (pseudorange positioning)
├── SatPVT.cpp           # Satellite orbit and clock computation
├── CoordinateConvert.cpp # Coordinate transformations (XYZ/BLH/ENU)
├── detect.cpp           # Outlier detection and cycle slip detection
├── TimeConvert.cpp      # Time system conversions
├── sockets.cpp          # TCP socket wrapper for real-time data
└── data/                # Sample data files
```

## Build Requirements

- **Compiler**: MinGW-w64 GCC or MSVC with C++11 support
- **Dependencies**: 
  - [Eigen3](https://eigen.tuxfamily.org/) (linear algebra library)
  - Windows Socket API (`ws2_32.lib`)
- **OS**: Windows (due to WinSock dependency)

## Build Instructions

### Using VS Code

1. Open the project in VS Code
2. Press `Ctrl+Shift+B` to build
3. The build task is configured in `.vscode/tasks.json`

### Manual Build

```bash
cd RTK_Raw
g++ -fdiagnostics-color=always -g \
    RTK.cpp DecodeNovOem7.cpp CoordinateConvert.cpp \
    TimeConvert.cpp SatPVT.cpp SPP.cpp detect.cpp \
    -I "D:\GNSS Algorithm\RTK\RTK\eigen-5.0.0\eigen-5.0.0" \
    -o RTK.exe \
    -lws2_32
```

## Usage

### Configuration

Edit `RTK_Structs.h` to set operation mode:

```cpp
#define FILEMODE 1   // 1 = File mode, 0 = Real-time TCP mode
```

### File Mode

Place the NovAtel OEM7 binary data file at:
```
D:\GNSS Algorithm\RTK\RTK\oem719-202603111200.bin
```

### Real-time Mode

Configure the TCP server address in `RTK.cpp`:
```cpp
if(OpenSocket(NetGps, "47.114.134.129", 7190) == false)
```

### Run

```bash
./RTK.exe
```

Results will be saved to `result.txt`.

## Output Format

### Satellite Information (per epoch)
```
G01 X= -12694065.000 Y=   8930477.000 Z=  21578501.000 Clk=-7.123456e-04 Vx= -945.6000 Vy=  2905.7000 Vz=  1373.1000 Clkd= 3.45678e-11 PIF=20200356.1000 Trop=  2.100 E= 45.123deg
```

### Position Solution (per epoch)
```
SPP: 2294 273600.000 X:-1269406.5000 Y:893047.7000 Z:2157850.1000 B:   30.12345678 L:  120.98765432 H:  123.456 Vx:   0.0000 Vy:   0.0000 Vz:   0.0000 GPS Clk:     12.345 BDS Clk:      0.000 PDOP:   1.234 Sigma:   2.345 GPSSats:  8 BDSSats:  6 Sats: 14
```

## Algorithm Details

### 1. Data Decoding
- **Protocol**: NovAtel OEM7 binary format
- **Messages**: RANGEB (observations), GPSEPHEM/BDSEPHEM (ephemeris), PSRPOS (reference position)
- **Observations**: Dual-frequency pseudorange, carrier phase, Doppler measurements

### 2. Outlier Detection
- **MW Combination**: Wide-lane ambiguity detection (threshold: 3.0 cycles)
- **GF Combination**: Ionospheric delay variation detection (threshold: 0.05m)
- **Dual-frequency PIF**: Ionosphere-free pseudorange for positioning

### 3. SPP Algorithm
- **Estimator**: Least Squares with iterative linearization
- **State Vector**: [X, Y, Z, dt_GPS, dt_BDS] (3D position + dual-system clock biases)
- **Convergence**: 25 iterations max, 0.1mm state change threshold
- **Quality Metrics**: PDOP, Sigma (post-fit RMS)

### 4. SPV Algorithm
- **Measurements**: Doppler observations converted to range-rate
- **State Vector**: [Vx, Vy, Vz, dtr_dot] (3D velocity + receiver clock drift)
- **Models**: Satellite velocity from broadcast ephemeris, LOS geometry

### 5. Error Corrections
- **Troposphere**: Hopfield model with standard meteorological parameters
- **TGD**: BDS B1/B3 ionosphere-free TGD correction
- **Relativity**: Schwarzschild term in satellite clock correction
- **Earth Rotation**: Sagnac effect during signal propagation

## Coordinate Systems

| System | Definition | Reference |
|--------|-----------|-----------|
| WGS-84 | GPS default coordinate system | GPS ICD |
| CGCS2000 | BeiDou coordinate system | BDS ICD |
| ENU | Local North-East-Up frame | Site-centered |

## Physical Constants

| Parameter | Value | Source |
|-----------|-------|--------|
| Speed of Light | 299,792,458 m/s | IAU 1976 |
| WGS-84 a | 6,378,137.0 m | WGS-84 |
| WGS-84 f | 1/298.257223563 | WGS-84 |
| GM (Earth) | 398600.5×10⁹ m³/s² | WGS-84 |
| Ω (WGS) | 7.2921151467×10⁻⁵ rad/s | WGS-84 |
| Ω (BDS) | 7.2921150×10⁻⁵ rad/s | CGCS2000 |

## Project Structure

```
RTK_Raw/
├── .vscode/                 # VS Code configuration
│   ├── c_cpp_properties.json
│   ├── launch.json
│   ├── settings.json
│   └── tasks.json
├── OEM719原始数据/          # Raw data samples
├── data/                    # Processed data outputs
├── .gitignore              # Git ignore rules
├── CoordinateConvert.cpp   # Coordinate transformations
├── DecodeNovOem7.cpp       # OEM7 protocol decoder
├── detect.cpp              # Cycle slip detection
├── RTK.cpp                 # Main program
├── RTK_Structs.h           # Data structures & constants
├── SatPVT.cpp              # Satellite PVT computation
├── SPP.cpp                 # SPP/SPV algorithms
├── sockets.cpp             # Network socket wrapper
├── test_CompSatElAz.cpp    # Unit test for elevation/azimuth
└── TimeConvert.cpp         # Time system conversions
```

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Acknowledgments

- This project is developed at the School of Geodesy and Geomatics, Wuhan University
- Broadcast ephemeris algorithms based on GPS ICD 200 and BeiDou ICD

## References

1. IS-GPS-200N. *Navstar GPS Space Segment/Navigation User Segment Interface*
2. 北斗卫星导航系统空间信号接口控制文件 - 公开服务信号B1I (Version 3.0)
3. NovAtel OEM7 Firmware Reference Manual
4. Hofmann-Wellenhof et al., *GNSS – Global Navigation Satellite Systems*, Springer

## License

MIT License - See [LICENSE](LICENSE) for details.

---

**Author**: RprPayador
**Institution**:
**Email**: 17707234049lijy@gmail.com  
**Version**: 1.0
