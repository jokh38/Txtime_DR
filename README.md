# Treatment Time Calculator with Dynamic Range Analysis

A GUI-based tool for calculating proton therapy treatment times and analyzing Dynamic Range effects in RayStation v2025.

## 📋 Overview

This program is a comprehensive GUI application designed for RayStation v2025 to calculate total treatment time for proton therapy and analyze treatment time variations according to Dynamic Range parameters. It includes respiratory gating functionality to accurately predict treatment times in real clinical environments.

## ✨ Key Features

### Core Functionality
- **Dynamic Range Analysis**: Calculate and compare treatment times across various DR values
- **Respiratory Gating Calculation**: Real-time respiratory signal-based gating efficiency analysis
- **Treatment Time Optimization**: Analysis of beam-on-time and total treatment time
- **Results Visualization**: Real-time canvas-based graph display

### Respiratory Signal Patterns
- **Sine**: Basic sinusoidal pattern
- **Tail**: Asymmetric respiratory curve
- **Long-tail**: Extended expiration phase pattern

### Gating Options
- **Gating Window**: 5%, 10%, 15%, 20% options
- **Phase/Amplitude** based gating
- **Real-time efficiency calculation**

## 🚀 Installation & Setup

### System Requirements
- **RayStation** v2025 (mandatory)
- **Python 3.7** or higher
- **Windows** environment recommended

## 📖 Usage Guide

### Basic Settings
| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| **Plan Selection** | Select treatment plan to analyze | First plan |
| **Resp. Period** | Respiratory period (seconds) | 4 sec |
| **Layer Switch Time** | Layer switching time (seconds) | 2 sec |

### DR Analysis Range
| Parameter | Description | Recommended Value |
|-----------|-------------|-------------------|
| **DR Start** | Starting DR value | 10 |
| **DR End** | Ending DR value | 200 |
| **DR Step** | Increment step | 20 |

### Analysis Procedure
1. **Plan Selection**: Choose treatment plan from dropdown menu
2. **Parameter Setup**: Configure respiratory period, gating window, etc.
3. **DR Range Setting**: Input Dynamic Range analysis range
4. **Run Analysis**: Click "Analyze DR Range" button
5. **Review Results**: Examine graphs and numerical results
6. **Export Results**: Save results to CSV file using "Export Results" button

## 🏗️ Code Architecture

### Class Structure
```
GUIProgram
├── Initialization & UI Creation
├── Respiratory Signal Processing
├── Gating Calculations
├── DR Analysis Engine
└── Results Visualization
```

### Key Methods
| Method | Functionality |
|--------|---------------|
| `resp_signal()` | Generate respiratory signal patterns |
| `update_gating_calculation()` | Calculate gating parameters |
| `f_bot_computation_with_dr()` | Calculate BOT with DR considerations |
| `calculate_total_treatment_time()` | Calculate total treatment time |
| `analyze_dr_range()` | Execute DR range analysis |

## 📊 Output Results

### Real-time Calculation Results
- **T_on**: Beam-on time (seconds)
- **T_off**: Beam-off time (seconds)
- **Efficiency**: Gating efficiency ratio

### Analysis Results (CSV Export)
```csv
DR,Total_Time,Beam_On_Time,Efficiency
10,1245.67,856.34,0.687
30,1089.23,856.34,0.786
50,978.45,856.34,0.875
...
```

### Visualization Features
- **Respiratory Signal Display**: Real-time respiratory pattern visualization
- **DR vs Treatment Time Graph**: Interactive analysis results plotting
- **Progress Monitoring**: Real-time analysis progress indicator

## ⚙️ Technical Details

### Calculation Algorithm
1. **Line-Scanning Analysis**: Process each energy layer's line segment weights
2. **Dose Rate Optimization**: Apply machine-specific dose rate limits
3. **Dynamic Range Application**: Adjust scanning parameters based on DR values
4. **Gating Integration**: Incorporate respiratory gating timing into calculations
5. **Time Resolution**: Apply 0.1ms time resolution for accurate results

### Machine Parameters
- **Maximum Scanning Speed**: 2000 cm/sec
- **Minimum Dose Rate**: 1.4 Gy/min
- **Energy Range**: 70-230 MeV (0.4 MeV steps)
- **Time Resolution**: 0.0001 seconds

## 🔧 Configuration Options

### Respiratory Gating Settings
| Setting | Options | Impact |
|---------|---------|---------|
| **Pattern Type** | Sine, Tail, Long-tail | Affects gating efficiency |
| **Gating Window** | 5%-20% | Determines beam-on duration |
| **Phase/Amplitude** | Toggle option | Changes gating reference |



## 📞 Support
Please contact me at kwanghyun.jo@samsung.com

For technical support or clinical questions:
- Review troubleshooting section
- Check RayStation documentation
- Validate input parameters
- Ensure proper system configuration
