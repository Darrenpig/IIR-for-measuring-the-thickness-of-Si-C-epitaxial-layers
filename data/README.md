# 数据目录说明

本目录用于存储红外干涉法测量碳化硅外延层厚度项目的数据文件。

## 目录结构

```
data/
├── raw/                    # 原始数据
│   ├── sic_samples/       # SiC样品数据
│   ├── si_reference/      # 硅参考样品数据
│   └── calibration/       # 校准数据
├── processed/             # 处理后数据
│   ├── filtered/          # 滤波后数据
│   ├── smoothed/          # 平滑后数据
│   └── normalized/        # 归一化数据
└── README.md              # 本说明文件
```

## 数据格式说明

### 原始数据格式
- **Excel文件** (.xlsx, .xls): 包含波数和反射率两列数据
- **CSV文件** (.csv): 逗号分隔的文本文件
- **文本文件** (.txt, .dat): 空格或制表符分隔的数据
- **MAT文件** (.mat): MATLAB数据文件

### 数据列说明
1. **第一列**: 波数 (cm⁻¹)
2. **第二列**: 反射率 (无量纲，0-1之间)

### 文件命名规范
- SiC样品: `sic_sample_[编号]_[角度]deg_[日期].xlsx`
- 硅参考: `si_reference_[角度]deg_[日期].xlsx`
- 校准数据: `calibration_[类型]_[日期].xlsx`

## 使用说明

1. 将原始测量数据放入 `raw/` 对应子目录
2. 使用 `data_loader.m` 函数加载数据
3. 处理后的数据自动保存到 `processed/` 对应子目录

## 注意事项

- 确保数据文件格式正确
- 波数应为递增或递减排列
- 反射率值应在合理范围内(0-1.2)
- 避免数据文件包含缺失值或异常值

## 示例数据加载

```matlab
% 加载SiC样品数据
[data, info] = data_loader('data/raw/sic_samples/sic_sample_01_45deg_20231201.xlsx');

% 加载硅参考数据
[si_data, si_info] = data_loader('data/raw/si_reference/si_reference_45deg_20231201.xlsx');
```