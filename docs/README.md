# 红外干涉法测量碳化硅外延层厚度 - 项目文档

## 项目概述

本项目基于红外干涉法原理，开发了一套完整的碳化硅(SiC)外延层厚度测量系统。通过数学建模和算法实现，提供了从理论分析到实际测量的完整解决方案。

## 文档结构

```
docs/
├── README.md                    # 项目文档总览（本文件）
├── user_guide/                  # 用户指南
│   ├── installation.md         # 安装指南
│   ├── quick_start.md          # 快速开始
│   ├── user_manual.md          # 用户手册
│   └── troubleshooting.md      # 故障排除
├── technical/                   # 技术文档
│   ├── theory_background.md    # 理论背景
│   ├── algorithm_design.md     # 算法设计
│   ├── implementation.md       # 实现细节
│   └── performance_analysis.md # 性能分析
├── api/                        # API文档
│   ├── functions_reference.md  # 函数参考
│   ├── data_structures.md      # 数据结构
│   └── examples.md             # 示例代码
├── validation/                 # 验证文档
│   ├── test_cases.md          # 测试用例
│   ├── validation_results.md  # 验证结果
│   └── benchmark.md           # 基准测试
└── maintenance/               # 维护文档
    ├── development_guide.md   # 开发指南
    ├── coding_standards.md   # 编码规范
    └── version_history.md    # 版本历史
```

## 核心功能模块

### 1. 理论分析模块 (Problem 1)
- **菲涅尔公式计算**: 实现s偏振和p偏振的反射、透射系数计算
- **相位差分析**: 基于光程差的相位差计算和干涉级数确定
- **厚度关系推导**: 建立厚度与干涉条件的数学关系

### 2. 算法实现模块 (Problem 2)
- **厚度测量算法**: 多种算法实现厚度确定
  - 相邻极值法
  - 干涉级数拟合法
  - 傅里叶变换法
- **可靠性分析**: 统计分析和精度评估
- **SiC数据处理**: 专门针对SiC材料的数据处理流程

### 3. 多光束干涉分析模块 (Problem 3)
- **多光束条件判断**: 自动识别多光束干涉条件
- **硅片数据分析**: 硅基底的光学特性分析
- **适用性评估**: 评估测量方法的适用范围

## 技术特点

### 算法优势
- **多算法融合**: 结合多种厚度测量算法，提高测量精度
- **自适应处理**: 根据数据特征自动选择最优处理方法
- **鲁棒性设计**: 具备良好的抗噪声和异常值处理能力

### 数据处理
- **多格式支持**: 支持Excel、CSV、文本等多种数据格式
- **智能预处理**: 自动数据清洗、平滑和异常值检测
- **实时分析**: 支持实时数据处理和结果反馈

### 可视化功能
- **专业绘图**: 提供光谱分析专用的可视化工具
- **交互式分析**: 支持交互式数据探索和参数调整
- **报告生成**: 自动生成专业的分析报告

## 系统要求

### 软件环境
- **MATLAB**: R2018b或更高版本
- **工具箱**: Signal Processing Toolbox, Curve Fitting Toolbox
- **操作系统**: Windows 10/11, macOS 10.14+, Linux Ubuntu 18.04+

### 硬件要求
- **内存**: 最小4GB，推荐8GB以上
- **存储**: 至少1GB可用空间
- **处理器**: Intel i5或同等性能处理器

## 快速开始

### 1. 环境准备
```matlab
% 添加项目路径
addpath(genpath('path/to/project'));

% 检查依赖
check_dependencies();
```

### 2. 运行主程序
```matlab
% 启动主程序
main();
```

### 3. 加载数据
```matlab
% 加载光谱数据
[data, info] = data_loader('your_data_file.xlsx');
```

### 4. 执行分析
```matlab
% 运行厚度测量
result = determine_thickness_algorithm(data, incident_angle);
```

## 主要函数接口

### 核心计算函数
```matlab
% 菲涅尔公式计算
[r_s, r_p, t_s, t_p] = fresnel_formula(n1, n2, angle);

% 相位差计算
phase_diff = calculate_phase_difference(wavenumber, thickness, n_film, angle);

% 厚度算法
result = determine_thickness_algorithm(data, angle, options);

% 可靠性分析
reliability = analyze_reliability(thickness_estimates, options);
```

### 数据处理函数
```matlab
% 数据加载
[data, info] = data_loader(filename, options);

% 信号处理
processed_data = signal_processor(data, options);

% 绘图工具
fig_handle = plot_spectral_data(data, options);
```

## 配置说明

### 物理常数配置
```matlab
% 加载物理常数
constants = load_constants();

% 主要常数
% constants.c - 光速
% constants.n_air - 空气折射率
% constants.n_si - 硅折射率
% constants.n_sic - SiC折射率
```

### 系统参数配置
```matlab
% 加载系统参数
params = load_parameters();

% 主要参数组
% params.measurement - 测量参数
% params.algorithm - 算法参数
% params.processing - 处理参数
```

## 数据格式

### 输入数据格式
```
波数(cm^-1)    反射率
800.0         0.234
801.0         0.236
...
```

### 输出结果格式
```matlab
result = struct(
    'thickness_estimates', [1.45, 1.48, 1.46],  % μm
    'final_thickness', 1.46,                     % μm
    'uncertainty', 0.02,                         % μm
    'reliability_score', 0.95,                   % 0-1
    'algorithm_used', 'combined',
    'processing_info', struct(...)
);
```

## 验证与测试

### 运行测试套件
```matlab
% 运行所有测试
test_results = test_all();

% 查看测试结果
disp(test_results);
```

### 性能基准测试
```matlab
% 运行性能测试
benchmark_results = run_benchmark();
```

## 常见问题

### Q1: 如何处理噪声数据？
A: 使用signal_processor函数进行预处理，包括平滑、去噪等操作。

### Q2: 厚度测量精度如何？
A: 在理想条件下，测量精度可达到±0.01μm，实际精度取决于数据质量。

### Q3: 支持哪些数据格式？
A: 支持Excel(.xlsx)、CSV(.csv)、文本(.txt)和MATLAB(.mat)格式。

## 技术支持

### 联系方式
- **项目主页**: [项目仓库链接]
- **问题反馈**: [Issue链接]
- **技术讨论**: [讨论区链接]

### 贡献指南
欢迎提交bug报告、功能建议和代码贡献。请参考开发指南了解详细信息。

## 许可证

本项目采用MIT许可证，详见LICENSE文件。

## 引用

如果您在研究中使用了本项目，请引用：

```
[项目引用信息]
```

---

**版本**: 1.0.0  
**最后更新**: 2024年12月  
**维护者**: 红外干