# 红外干涉法测量碳化硅外延层厚度

基于红外干涉法原理的碳化硅(SiC)外延层厚度测量系统，提供完整的数学建模、算法实现和数据分析解决方案。

## 项目概述

本项目针对碳化硅外延层厚度测量的实际需求，基于红外干涉法的物理原理，开发了一套完整的测量分析系统。系统包含理论建模、算法实现、数据处理和结果分析等完整功能模块。

### 主要特点

- 🔬 **理论完备**: 基于菲涅尔公式和干涉理论的完整数学建模
- 🧮 **算法多样**: 集成多种厚度测量算法，提高测量精度和可靠性
- 📊 **数据处理**: 支持多种数据格式，提供专业的信号处理功能
- 📈 **可视化**: 专业的光谱数据可视化和分析工具
- 🔍 **质量控制**: 完善的可靠性分析和误差评估机制
- 🚀 **易于使用**: 清晰的模块化设计，简单易用的接口

## 项目结构

```
B/
├── main.m                    # 主程序入口
├── README.md                 # 项目说明（本文件）
├── config/                   # 配置文件
│   ├── constants.m          # 物理常数定义
│   └── parameters.m         # 系统参数配置
├── models/                   # 核心算法模块
│   ├── problem1/            # 理论分析模块
│   │   ├── fresnel_formula.m
│   │   ├── phase_difference.m
│   │   └── thickness_relation.m
│   ├── problem2/            # 算法实现模块
│   │   ├── thickness_algorithm.m
│   │   ├── reliability_analysis.m
│   │   └── sic_data_processor.m
│   └── problem3/            # 多光束干涉分析
│       ├── multi_beam_conditions.m
│       └── si_data_analyzer.m
├── utils/                    # 工具函数
│   ├── data_loader.m        # 数据加载器
│   ├── signal_processor.m   # 信号处理工具
│   ├── plot_utils.m         # 绘图工具
│   └── math_utils.m         # 数学工具
├── data/                     # 数据目录
│   ├── 附件1.xlsx           # SiC晶圆片10°入射角数据
│   ├── 附件2.xlsx           # SiC晶圆片15°入射角数据
│   ├── 附件3.xlsx           # Si晶圆片10°入射角数据
│   ├── 附件4.xlsx           # Si晶圆片15°入射角数据
│   └── README.md            # 数据说明
├── results/                  # 计算结果
│   ├── problem1/            # 问题一结果
│   ├── problem2/            # 问题二结果
│   ├── problem3/            # 问题三结果
│   └── README.md            # 结果说明
├── tests/                    # 测试文件
│   └── test_all.m           # 综合测试
└── docs/                     # 项目文档
    └── README.md            # 文档说明
```

## 快速开始

### 1. 环境要求

- **MATLAB**: R2018b或更高版本
- **推荐工具箱**: Signal Processing Toolbox, Curve Fitting Toolbox
- **系统要求**: Windows/macOS/Linux，4GB+内存

### 2. 运行项目

```matlab
% 1. 启动MATLAB并切换到项目目录
cd('path/to/project/B');

% 2. 运行主程序
main();
```

### 3. 基本使用

```matlab
% 加载数据
[data, info] = data_loader('data/附件1.xlsx');

% 执行厚度测量
angle = 45; % 入射角度（度）
result = determine_thickness_algorithm(data, angle);

% 查看结果
fprintf('测量厚度: %.3f μm\n', result.final_thickness);
fprintf('测量不确定度: %.3f μm\n', result.uncertainty);
```

## 核心功能

### 问题一：理论分析

基于菲涅尔公式和干涉理论，建立厚度测量的数学模型：

- **菲涅尔公式计算**: 计算界面反射和透射系数
- **相位差分析**: 分析光程差和干涉条件
- **厚度关系推导**: 建立厚度与光谱特征的关系

```matlab
% 运行理论分析
run_problem1();
```

### 问题二：算法实现

实现多种厚度测量算法，提供可靠性分析：

- **相邻极值法**: 基于干涉极值的厚度计算
- **干涉级数拟合法**: 通过拟合干涉级数确定厚度
- **傅里叶变换法**: 频域分析方法
- **可靠性评估**: 统计分析和精度评估

```matlab
% 运行算法实现
run_problem2();
```

### 问题三：多光束干涉分析

分析多光束干涉条件，评估测量方法适用性：

- **多光束条件判断**: 自动识别多光束干涉
- **硅基底分析**: 分析硅基底对测量的影响
- **适用性评估**: 评估方法的适用范围

```matlab
% 运行多光束分析
run_problem3();
```

## 数据格式

### 输入数据

支持多种格式的光谱数据：

- **Excel格式** (.xlsx): 推荐格式，支持多工作表
- **CSV格式** (.csv): 通用文本格式
- **文本格式** (.txt): 空格或制表符分隔
- **MATLAB格式** (.mat): MATLAB原生格式

数据应包含两列：
```
波数(cm^-1)    反射率
800.0         0.234
801.0         0.236
...
```

### 输出结果

```matlab
result = struct(
    'thickness_estimates', [1.45, 1.48, 1.46],  % 各算法估计值 (μm)
    'final_thickness', 1.46,                     % 最终厚度 (μm)
    'uncertainty', 0.02,                         % 不确定度 (μm)
    'reliability_score', 0.95,                   % 可靠性评分 (0-1)
    'algorithm_used', 'combined',                % 使用的算法
    'processing_info', struct(...)               % 处理信息
);
```

## 配置说明

### 物理常数

在 `config/constants.m` 中定义了所需的物理常数：

```matlab
constants = load_constants();
% constants.c      - 光速
% constants.n_air  - 空气折射率
% constants.n_si   - 硅折射率
% constants.n_sic  - SiC折射率
```

### 系统参数

在 `config/parameters.m` 中配置系统参数：

```matlab
params = load_parameters();
% params.measurement - 测量参数
% params.algorithm   - 算法参数
% params.processing  - 处理参数
```

## 测试与验证

### 运行测试

```matlab
% 运行完整测试套件
test_results = test_all();

% 查看测试结果
fprintf('测试通过率: %.1f%%\n', ...
    test_results.passed_tests / test_results.total_tests * 100);
```

### 性能指标

- **测量精度**: ±0.01-0.05 μm（取决于数据质量）
- **处理速度**: 典型光谱数据 < 1秒
- **可靠性**: 95%以上的测量可靠性评分

## 示例应用

### 单次测量

```matlab
% 加载SiC外延层光谱数据
[data, ~] = data_loader('data/附件1.xlsx');

% 设置测量参数
angle = 45;  % 入射角
options = struct('method', 'combined', 'smooth', true);

% 执行测量
result = determine_thickness_algorithm(data, angle, options);

% 显示结果
fprintf('SiC外延层厚度: %.3f ± %.3f μm\n', ...
    result.final_thickness, result.uncertainty);
```

### 批量处理

```matlab
% 批量处理多个样品
data_files = dir('data/*.xlsx');
results = cell(length(data_files), 1);

for i = 1:length(data_files)
    filename = fullfile(data_files(i).folder, data_files(i).name);
    [data, ~] = data_loader(filename);
    results{i} = determine_thickness_algorithm(data, 45);
    
    fprintf('样品 %d: %.3f μm\n', i, results{i}.final_thickness);
end
```

## 常见问题

### Q: 如何提高测量精度？
A: 
1. 确保光谱数据质量良好，信噪比高
2. 使用适当的数据预处理（平滑、去噪）
3. 选择合适的入射角度
4. 使用多算法融合结果

### Q: 如何处理噪声数据？
A: 使用 `signal_processor` 函数进行预处理：
```matlab
options = struct('denoise', true, 'smooth', true, 'smooth_window', 5);
processed_data = signal_processor(raw_data, options);
```

### Q: 支持哪些SiC材料类型？
A: 系统主要针对4H-SiC和6H-SiC多型，可通过调整折射率参数适应其他类型。

## 技术支持

- 📧 **技术咨询**: [联系邮箱]
- 📖 **详细文档**: 查看 `docs/` 目录
- 🐛 **问题反馈**: [Issue链接]
- 💬 **技术讨论**: [讨论区链接]

## 贡献

欢迎提交问题报告、功能建议和代码贡献！

1. Fork 项目
2. 创建功能分支
3. 提交更改
4. 发起 Pull Request

## 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件。

## 引用

如果您在研究中使用了本项目，请引用：

```bibtex
@software{sic_thickness_measurement,
  title={红外干涉法测量碳化硅外延层厚度系统},
  author={项目团队},
  year={2024},
  url={项目链接}
}
```

---

**版本**: 1.0.0  
**最后更新**: 2024年12月  
**开发团队**: 红外干涉法测量项目组