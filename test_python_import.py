#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
测试Python数据导入功能
"""

import sys
import os
sys.path.append('code')

from data_process import DataProcessor

def test_python_import():
    """测试Python版本的数据导入"""
    try:
        dp = DataProcessor()
        
        # 测试数据加载
        print("正在测试Python数据导入...")
        success = dp.load_excel_data('附件1.xlsx')
        
        if not success:
            raise Exception("数据加载失败")
            
        # 获取加载的数据
        data = dp.raw_data
        
        print(f"Python数据读取成功")
        print(f"波数数据点数: {len(data['wavenumber'])}")
        print(f"反射率数据点数: {len(data['reflectance'])}")
        print(f"波数范围: {min(data['wavenumber']):.2f} - {max(data['wavenumber']):.2f} cm^-1")
        print(f"反射率范围: {min(data['reflectance']):.4f} - {max(data['reflectance']):.4f}")
        
        # 运行质量检查
        print("\n运行数据质量检查...")
        quality_report = dp.quality_check()
        print(f"数据质量检查完成")
        
        print(f"\n数据统计信息:")
        print(f"- 数据点数: {quality_report['total_points']}")
        print(f"- 缺失值: {quality_report['missing_values']}")
        print(f"- 异常值: {quality_report['outliers']} ({quality_report['outlier_ratio']:.1f}%)")
        print(f"- 信噪比: {quality_report['snr_db']:.1f} dB")
        print(f"- 质量等级: {quality_report['quality_grade']}")
        print(f"- 噪声水平: {quality_report['noise_level']:.4f}")
        print(f"- 不规则步长: {quality_report['irregular_steps']}")
        
        return True, quality_report
        
    except Exception as e:
        print(f"Python导入错误: {e}")
        return False, str(e)

if __name__ == "__main__":
    success, result = test_python_import()
    if success:
        print("\nPython数据导入测试通过！")
    else:
        print(f"\nPython数据导入测试失败: {result}")