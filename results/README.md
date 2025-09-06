# 计算结果目录

本目录用于存储红外干涉法测量碳化硅外延层厚度项目的计算结果。

## 目录结构

```
results/
├── README.md              # 本说明文件
├── problem1/              # 问题一结果
│   ├── fresnel_results/   # 菲涅尔公式计算结果
│   ├── phase_results/     # 相位差计算结果
│   └── theory_results/    # 理论推导结果
├── problem2/              # 问题二结果
│   ├── thickness_results/ # 厚度测量结果
│   ├── reliability_results/ # 可靠性分析结果
│   └── sic_analysis/      # SiC数据分析结果
├── problem3/              # 问题三结果
│   ├── multi_beam_analysis/ # 多光束干涉分析
│   ├── si_analysis/       # 硅片数据分析
│   └── comparison_results/ # 对比分析结果
└── reports/               # 综合报告
    ├── summary_reports/   # 总结报告
    ├── figures/          # 图表文件
    └── data_exports/     # 导出数据
```

## 结果文件格式

### MATLAB结果文件 (.mat)
- 包含完整的计算结果结构体
- 便于后续MATLAB程序读取和处理
- 保留所有数值精度

### Excel结果文件 (.xlsx)
- 适合数据查看和进一步分析
- 包含多个工作表，分别存储不同类型的结果
- 便于与其他软件交互

### 图像文件 (.png, .fig)
- .png格式：用于报告和展示
- .fig格式：MATLAB图形文件，可重新编辑

### 文本报告 (.txt, .md)
- 包含分析总结和关键结论
- Markdown格式便于文档管理

## 文件命名规范

### 基本格式
```
[项目类型]_[分析方法]_[日期时间].[扩展名]
```

### 示例
- `problem1_fresnel_20241201_143022.mat`
- `problem2_thickness_algorithm_20241201_143022.xlsx`
- `problem3_multi_beam_analysis_20241201_143022.png`
- `reliability_analysis_summary_20241201.md`

### 命名说明
- **项目类型**: problem1, problem2, problem3, summary
- **分析方法**: fresnel, phase_diff, thickness, reliability, multi_beam等
- **日期时间**: YYYYMMDD_HHMMSS格式
- **扩展名**: 根据文件类型确定

## 结果数据结构

### 问题一结果结构
```matlab
result_problem1 = struct(
    'fresnel_coefficients', struct(...),  % 菲涅尔系数
    'phase_differences', struct(...),     % 相位差计算
    'thickness_relation', struct(...),    % 厚度关系
    'theoretical_curves', struct(...),    % 理论曲线
    'analysis_parameters', struct(...),   % 分析参数
    'computation_info', struct(...)       % 计算信息
);
```

### 问题二结果结构
```matlab
result_problem2 = struct(
    'thickness_estimates', [],            % 厚度估计值
    'algorithm_results', struct(...),     % 算法结果
    'reliability_analysis', struct(...),  % 可靠性分析
    'sic_data_analysis', struct(...),     % SiC数据分析
    'processing_info', struct(...),       % 处理信息
    'quality_metrics', struct(...)        % 质量指标
);
```

### 问题三结果结构
```matlab
result_problem3 = struct(
    'multi_beam_judgment', logical,       % 多光束判断
    'si_data_analysis', struct(...),      % 硅片分析
    'comparison_results', struct(...),    % 对比结果
    'applicability_assessment', struct(...), % 适用性评估
    'recommendations', struct(...)        % 建议
);
```

## 使用说明

### 保存结果
```matlab
% 保存MATLAB结果
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('problem1_analysis_%s.mat', timestamp);
save(fullfile('results', 'problem1', filename), 'result_problem1');

% 导出Excel结果
excel_filename = sprintf('problem2_thickness_%s.xlsx', timestamp);
export_to_excel(result_problem2, fullfile('results', 'problem2', excel_filename));
```

### 加载结果
```matlab
% 加载最新结果
result_files = dir(fullfile('results', 'problem1', '*.mat'));
if ~isempty(result_files)
    [~, idx] = max([result_files.datenum]);
    latest_file = result_files(idx).name;
    load(fullfile('results', 'problem1', latest_file));
end
```

### 生成报告
```matlab
% 生成综合报告
generate_summary_report(result_problem1, result_problem2, result_problem3);
```

## 数据备份

### 自动备份
- 重要结果会自动创建备份副本
- 备份文件添加 `_backup` 后缀
- 定期清理过期备份文件

### 手动备份
```matlab
% 备份重要结果
backup_results('results', 'backup_folder');
```

## 注意事项

1. **文件大小**: 注意大型结果文件的存储空间占用
2. **版本管理**: 重要结果建议保留多个版本
3. **数据完整性**: 定期检查结果文件的完整性
4. **访问权限**: 确保结果文件的适当访问权限
5. **清理策略**: 定期清理临时和过期的结果文件

## 质量控制

### 结果验证
- 自动检查结果的合理性
- 与理论值进行对比验证
- 记录异常结果和处理方法

### 追溯性
- 记录完整的计算参数
- 保存输入数据的引用
- 记录计算环境和版本信息

---

**最后更新**: 2024年12月
**维护者**: 红外干涉法测量项目组